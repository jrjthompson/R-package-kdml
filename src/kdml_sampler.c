#include "kdml.h"

#include <R_ext/Utils.h>

#include <limits.h>
#include <math.h>
#include <string.h>

#define KDML_ADAPT_GAIN 0.25
#define KDML_ADAPT_POWER 0.6
#define KDML_MIN_PROPOSAL_SD 1e-3
#define KDML_MAX_PROPOSAL_SD 10.0
#define KDML_TWO_PI 6.283185307179586476925286766559

static uint64_t kdml_rng_next(kdml_rng *rng)
{
    uint64_t value = rng->state;
    value ^= value >> 12;
    value ^= value << 25;
    value ^= value >> 27;
    rng->state = value;
    return value * UINT64_C(2685821657736338717);
}

static int kdml_cancel_requested(volatile int *cancel)
{
    int value = 0;
    if (cancel == NULL) {
        return 0;
    }
#ifdef _OPENMP
# pragma omp atomic read
#endif
    value = *cancel;
    return value != 0;
}

static void kdml_publish_completed(volatile int *completed_sweeps, int value)
{
    if (completed_sweeps == NULL) {
        return;
    }
#ifdef _OPENMP
# pragma omp atomic write
#endif
    *completed_sweeps = value;
}

static double kdml_rng_uniform(kdml_rng *rng)
{
    return ((double) (kdml_rng_next(rng) >> 11) + 0.5) /
           9007199254740992.0;
}

static double kdml_rng_normal(kdml_rng *rng)
{
    double radius;
    double angle;

    if (rng->has_spare) {
        rng->has_spare = 0;
        return rng->spare;
    }
    radius = sqrt(-2.0 * log(kdml_rng_uniform(rng)));
    angle = KDML_TWO_PI * kdml_rng_uniform(rng);
    rng->spare = radius * sin(angle);
    rng->has_spare = 1;
    return radius * cos(angle);
}

static void kdml_count_invalid(int *count)
{
    if (*count < INT_MAX) {
        ++(*count);
    }
}

static int kdml_metropolis_accept(double proposed, double current,
                                  kdml_rng *rng)
{
    const double log_ratio = proposed - current;
    if (!R_FINITE(proposed)) {
        return 0;
    }
    if (log_ratio >= 0.0) {
        return 1;
    }
    return log(kdml_rng_uniform(rng)) < log_ratio;
}

static int kdml_backend_build(kdml_score_backend *backend,
                              const int *kernel_code,
                              const double *bandwidth, double *score)
{
    backend->error[0] = '\0';
    return backend->ops->build(
        backend, kernel_code, bandwidth, score
    );
}

static int kdml_backend_propose(kdml_score_backend *backend,
                                const int *current_kernel,
                                const double *current_bandwidth,
                                int feature, int proposed_kernel,
                                double proposed_bandwidth, double *score)
{
    backend->error[0] = '\0';
    return backend->ops->propose(
        backend, current_kernel, current_bandwidth,
        feature, proposed_kernel, proposed_bandwidth, score
    );
}

static int kdml_backend_accept(kdml_score_backend *backend,
                               int feature, int proposed_kernel,
                               double proposed_bandwidth)
{
    backend->error[0] = '\0';
    return backend->ops->accept(
        backend, feature, proposed_kernel, proposed_bandwidth
    );
}

static int kdml_propose_other_kernel(const kdml_mcmc_config *config,
                                     int feature, int current,
                                     kdml_rng *rng)
{
    const int begin = config->candidate_offsets[feature];
    const int end = config->candidate_offsets[feature + 1];
    const int alternatives = end - begin - 1;
    int target;
    int position;

    if (alternatives <= 0) {
        return current;
    }
    target = (int) floor(
        kdml_rng_uniform(rng) * (double) alternatives
    );
    if (target >= alternatives) {
        target = alternatives - 1;
    }
    for (position = begin; position < end; ++position) {
        const int candidate = config->candidate_codes[position];
        if (candidate == current) {
            continue;
        }
        if (target == 0) {
            return candidate;
        }
        --target;
    }
    return current;
}

static void kdml_adapt_scale(double *proposal_sd, int count,
                             int accepted, double target)
{
    double log_scale = log(*proposal_sd);
    const double gain = KDML_ADAPT_GAIN *
                        pow((double) count, -KDML_ADAPT_POWER);
    log_scale += gain * ((accepted ? 1.0 : 0.0) - target);
    if (log_scale < log(KDML_MIN_PROPOSAL_SD)) {
        log_scale = log(KDML_MIN_PROPOSAL_SD);
    } else if (log_scale > log(KDML_MAX_PROPOSAL_SD)) {
        log_scale = log(KDML_MAX_PROPOSAL_SD);
    }
    *proposal_sd = exp(log_scale);
}

static double kdml_model_log_prior(const kdml_data *data,
                                   const kdml_mcmc_config *config)
{
    int feature;
    double answer = 0.0;
    for (feature = 0; feature < data->p; ++feature) {
        const int count = config->candidate_offsets[feature + 1] -
                          config->candidate_offsets[feature];
        answer -= log((double) count);
    }
    return answer;
}

static void kdml_store_draw(const kdml_data *data,
                            const kdml_mcmc_config *config,
                            kdml_mcmc_output *output,
                            int chain, int draw,
                            const int *kernel, const double *theta,
                            const double *bandwidth,
                            double score, double log_target)
{
    int feature;
    const R_xlen_t score_index =
        (R_xlen_t) draw + (R_xlen_t) config->draws * chain;

    output->score[score_index] = score;
    output->log_target[score_index] = log_target;
    for (feature = 0; feature < data->p; ++feature) {
        const R_xlen_t index =
            (R_xlen_t) draw +
            (R_xlen_t) config->draws *
                ((R_xlen_t) feature + (R_xlen_t) data->p * chain);
        output->kernel_code[index] = kernel[feature];
        output->theta[index] = theta[feature];
        output->bandwidth[index] = bandwidth[feature];
    }
}

void kdml_initialize_mcmc_output(const kdml_data *data,
                                 const kdml_mcmc_config *config,
                                 kdml_mcmc_output *output)
{
    const size_t diagnostics_count =
        (size_t) data->p * (size_t) config->chains;

    memset(output->kernel_attempts, 0,
           diagnostics_count * sizeof(double));
    memset(output->kernel_accepts, 0,
           diagnostics_count * sizeof(double));
    memset(output->bandwidth_attempts, 0,
           diagnostics_count * sizeof(double));
    memset(output->bandwidth_accepts, 0,
           diagnostics_count * sizeof(double));
    memset(output->sampling_kernel_attempts, 0,
           diagnostics_count * sizeof(double));
    memset(output->sampling_kernel_accepts, 0,
           diagnostics_count * sizeof(double));
    memset(output->sampling_bandwidth_attempts, 0,
           diagnostics_count * sizeof(double));
    memset(output->sampling_bandwidth_accepts, 0,
           diagnostics_count * sizeof(double));
    memset(output->invalid_proposals, 0,
           (size_t) config->chains * sizeof(int));
}

int kdml_run_chain(const kdml_data *data,
                   const kdml_mcmc_config *config,
                   kdml_mcmc_output *output,
                   kdml_score_backend *backend, kdml_progress *progress,
                   int chain, int *state_kernel, double *state_theta,
                   double *state_bandwidth, double *state_proposal_sd,
                   int *adapt_count, kdml_rng *rng,
                   volatile int *completed_sweeps, volatile int *cancel,
                   int check_interrupt)
{
    int feature;
    int sweep;
    int retained = 0;
    const int total_sweeps =
        config->warmup + config->draws * config->thin;
    const double log_model_prior = kdml_model_log_prior(data, config);
    double current_score;
    double current_log_prior;
    double current_log_target;

    if (progress != NULL) {
        kdml_progress_begin(progress, chain + 1);
    }
    memcpy(state_kernel, config->initial_kernel,
           (size_t) data->p * sizeof(int));
    memcpy(state_theta, config->initial_theta,
           (size_t) data->p * sizeof(double));
    memcpy(state_proposal_sd, config->proposal_sd,
           (size_t) data->p * sizeof(double));
    memset(adapt_count, 0, (size_t) data->p * sizeof(int));

    for (feature = 0; feature < data->p; ++feature) {
        if (!kdml_theta_to_bandwidth(
                data, feature, state_kernel[feature],
                state_theta[feature], &state_bandwidth[feature])) {
            return 1;
        }
    }
    if (!kdml_backend_build(
            backend, state_kernel, state_bandwidth, &current_score)) {
        return backend->error[0] == '\0' ? 1 : 4;
    }
    current_log_prior = kdml_log_prior(
        data, state_theta, config->prior_mean, config->prior_sd
    );
    current_log_target = kdml_log_target(
        current_score, config->beta, current_log_prior, log_model_prior
    );
    if (!R_FINITE(current_log_target)) {
        return 1;
    }

    for (sweep = 0; sweep < total_sweeps; ++sweep) {
        const int in_sampling = sweep >= config->warmup;

        if (kdml_cancel_requested(cancel)) {
            return 5;
        }
        if (check_interrupt && (sweep & 15) == 0) {
            R_CheckUserInterrupt();
        }

        for (feature = 0; feature < data->p; ++feature) {
            const int candidate_count =
                config->candidate_offsets[feature + 1] -
                config->candidate_offsets[feature];
            const size_t diagnostic =
                (size_t) feature + (size_t) data->p * (size_t) chain;

            if (candidate_count > 1) {
                const int proposed_kernel = kdml_propose_other_kernel(
                    config, feature, state_kernel[feature], rng
                );
                double proposed_bandwidth;
                int accepted = 0;

                output->kernel_attempts[diagnostic] += 1.0;
                if (in_sampling) {
                    output->sampling_kernel_attempts[diagnostic] += 1.0;
                }
                if (!kdml_theta_to_bandwidth(
                        data, feature, proposed_kernel,
                        state_theta[feature], &proposed_bandwidth)) {
                    kdml_count_invalid(&output->invalid_proposals[chain]);
                } else {
                    double proposed_score;
                    if (!kdml_backend_propose(
                            backend, state_kernel, state_bandwidth,
                            feature, proposed_kernel,
                            proposed_bandwidth, &proposed_score)) {
                        return 4;
                    }
                    {
                        const double proposed_target = kdml_log_target(
                            proposed_score, config->beta,
                            current_log_prior, log_model_prior
                        );
                        if (!R_FINITE(proposed_target)) {
                            kdml_count_invalid(
                                &output->invalid_proposals[chain]
                            );
                        } else if (kdml_metropolis_accept(
                                       proposed_target,
                                       current_log_target, rng)) {
                            accepted = 1;
                            if (!kdml_backend_accept(
                                    backend, feature, proposed_kernel,
                                    proposed_bandwidth)) {
                                return 4;
                            }
                            state_kernel[feature] = proposed_kernel;
                            state_bandwidth[feature] = proposed_bandwidth;
                            current_score = proposed_score;
                            current_log_target = proposed_target;
                        }
                    }
                }
                if (accepted) {
                    output->kernel_accepts[diagnostic] += 1.0;
                    if (in_sampling) {
                        output->sampling_kernel_accepts[diagnostic] += 1.0;
                    }
                }
            }

            {
                const double old_theta = state_theta[feature];
                const double proposed_theta = old_theta +
                    state_proposal_sd[feature] * kdml_rng_normal(rng);
                double proposed_bandwidth;
                double proposed_score = R_NegInf;
                double proposed_log_prior = R_NegInf;
                double proposed_target = R_NegInf;
                int accepted = 0;

                output->bandwidth_attempts[diagnostic] += 1.0;
                if (in_sampling) {
                    output->sampling_bandwidth_attempts[diagnostic] += 1.0;
                }
                if (kdml_theta_to_bandwidth(
                        data, feature, state_kernel[feature],
                        proposed_theta, &proposed_bandwidth)) {
                    if (!kdml_backend_propose(
                            backend, state_kernel, state_bandwidth,
                            feature, state_kernel[feature],
                            proposed_bandwidth, &proposed_score)) {
                        return 4;
                    }
                    if (R_FINITE(proposed_score)) {
                        state_theta[feature] = proposed_theta;
                        proposed_log_prior = kdml_log_prior(
                            data, state_theta,
                            config->prior_mean, config->prior_sd
                        );
                        state_theta[feature] = old_theta;
                        proposed_target = kdml_log_target(
                            proposed_score, config->beta,
                            proposed_log_prior, log_model_prior
                        );
                    }
                }

                if (!R_FINITE(proposed_target)) {
                    kdml_count_invalid(&output->invalid_proposals[chain]);
                } else if (kdml_metropolis_accept(
                               proposed_target, current_log_target, rng)) {
                    accepted = 1;
                    if (!kdml_backend_accept(
                            backend, feature, state_kernel[feature],
                            proposed_bandwidth)) {
                        return 4;
                    }
                    state_theta[feature] = proposed_theta;
                    state_bandwidth[feature] = proposed_bandwidth;
                    current_score = proposed_score;
                    current_log_prior = proposed_log_prior;
                    current_log_target = proposed_target;
                }

                if (accepted) {
                    output->bandwidth_accepts[diagnostic] += 1.0;
                    if (in_sampling) {
                        output->sampling_bandwidth_accepts[diagnostic] += 1.0;
                    }
                }
                if (!in_sampling && config->adapt) {
                    ++adapt_count[feature];
                    kdml_adapt_scale(
                        &state_proposal_sd[feature],
                        adapt_count[feature], accepted,
                        config->target_accept
                    );
                }
            }
        }

        if (((sweep + 1) & 127) == 0) {
            if (!kdml_backend_build(
                    backend, state_kernel, state_bandwidth,
                    &current_score)) {
                return backend->error[0] == '\0' ? 2 : 4;
            }
            current_log_prior = kdml_log_prior(
                data, state_theta, config->prior_mean, config->prior_sd
            );
            current_log_target = kdml_log_target(
                current_score, config->beta,
                current_log_prior, log_model_prior
            );
            if (!R_FINITE(current_log_target)) {
                return 2;
            }
        }

        if (in_sampling &&
            ((sweep + 1 - config->warmup) % config->thin == 0)) {
            kdml_store_draw(
                data, config, output, chain, retained,
                state_kernel, state_theta, state_bandwidth,
                current_score, current_log_target
            );
            ++retained;
        }
        kdml_publish_completed(completed_sweeps, sweep + 1);
        if (progress != NULL) {
            kdml_progress_tick(progress, sweep + 1);
        }
    }

    if (retained != config->draws) {
        return 3;
    }
    for (feature = 0; feature < data->p; ++feature) {
        const size_t diagnostic =
            (size_t) feature + (size_t) data->p * (size_t) chain;
        output->final_proposal_sd[diagnostic] =
            state_proposal_sd[feature];
    }
    if (progress != NULL) {
        kdml_progress_end(progress, "done");
    }
    return 0;
}

int kdml_run_mcmc(const kdml_data *data, const kdml_mcmc_config *config,
                  kdml_mcmc_output *output,
                  kdml_score_backend *backend, kdml_progress *progress,
                  int *state_kernel, double *state_theta,
                  double *state_bandwidth, double *state_proposal_sd,
                  int *adapt_count, kdml_rng *rng, int *error_chain)
{
    int chain;

    kdml_initialize_mcmc_output(data, config, output);
    for (chain = 0; chain < config->chains; ++chain) {
        const int status = kdml_run_chain(
            data, config, output, backend, progress, chain,
            state_kernel, state_theta, state_bandwidth,
            state_proposal_sd, adapt_count, &rng[chain],
            NULL, NULL, 1
        );
        if (status != 0) {
            *error_chain = chain;
            return status;
        }
    }
    return 0;
}
