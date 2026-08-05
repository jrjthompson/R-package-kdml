#include "kdml.h"

#include <R_ext/Utils.h>

#include <float.h>
#include <math.h>
#include <string.h>

#define KDML_PI 3.141592653589793238462643383279502884
#define KDML_SQRT_TWO 1.414213562373095048801688724209698079
#define KDML_LOG_TWO 0.693147180559945309417232121458176568
#define KDML_LOG_TWO_PI 1.837877066409345483560659472811235279

static double kdml_continuous_kernel(int code, double z)
{
    const double absolute = fabs(z);
    double base;
    double tail;

    switch (code) {
    case 0:
        return exp(-0.5 * z * z) / sqrt(2.0 * KDML_PI);
    case 1:
        return absolute <= 1.0 ? 0.75 * (1.0 - z * z) : 0.0;
    case 2:
        return absolute <= 1.0 ? 0.5 : 0.0;
    case 3:
        return absolute <= 1.0 ? 1.0 - absolute : 0.0;
    case 4:
        if (absolute > 1.0) {
            return 0.0;
        }
        base = 1.0 - z * z;
        return (15.0 / 16.0) * base * base;
    case 5:
        if (absolute > 1.0) {
            return 0.0;
        }
        base = 1.0 - z * z;
        return (35.0 / 32.0) * base * base * base;
    case 6:
        if (absolute > 1.0) {
            return 0.0;
        }
        base = 1.0 - absolute * absolute * absolute;
        return (70.0 / 81.0) * base * base * base;
    case 7:
        return absolute <= 1.0
                   ? (KDML_PI / 4.0) * cos(KDML_PI * z / 2.0)
                   : 0.0;
    case 8:
        tail = exp(-absolute);
        return tail / ((1.0 + tail) * (1.0 + tail));
    case 9:
        tail = exp(-absolute);
        return (2.0 * tail) / (KDML_PI * (1.0 + tail * tail));
    case 10:
        base = absolute / KDML_SQRT_TWO;
        return 0.5 * exp(-base) * sin(base + KDML_PI / 4.0);
    default:
        return R_NaReal;
    }
}

int kdml_kernel_is_valid(int type, int code)
{
    if (type == KDML_CONTINUOUS) {
        return code >= 0 && code <= 10;
    }
    if (type == KDML_NOMINAL) {
        return code >= 0 && code <= 1;
    }
    if (type == KDML_ORDINAL) {
        return code >= 0 && code <= 4;
    }
    return 0;
}

double kdml_bandwidth_upper(int type, int code, int level_count)
{
    if (type == KDML_CONTINUOUS) {
        return R_PosInf;
    }
    if (type == KDML_NOMINAL && code == 1) {
        return ((double) level_count - 1.0) / (double) level_count;
    }
    return 1.0;
}

int kdml_theta_to_bandwidth(const kdml_data *data, int feature,
                            int kernel_code, double theta,
                            double *bandwidth)
{
    const int type = data->type[feature];
    double value;
    double logistic;
    double upper;

    if (!R_FINITE(theta)) {
        return 0;
    }
    if (type == KDML_CONTINUOUS) {
        value = data->scale[feature] * exp(theta);
        if (!R_FINITE(value) || value <= 0.0) {
            return 0;
        }
        *bandwidth = value;
        return 1;
    }

    if (theta >= 0.0) {
        logistic = 1.0 / (1.0 + exp(-theta));
    } else {
        value = exp(theta);
        logistic = value / (1.0 + value);
    }
    upper = kdml_bandwidth_upper(type, kernel_code,
                                 data->level_count[feature]);
    value = upper * logistic;
    if (!R_FINITE(value) || value <= 0.0 || value >= upper) {
        return 0;
    }
    *bandwidth = value;
    return 1;
}

double kdml_feature_value(const kdml_data *data, int feature,
                          double first, double second,
                          int kernel_code, double bandwidth)
{
    const int type = data->type[feature];
    const double difference = first - second;
    const double distance = fabs(difference);
    double value;

    if (!R_FINITE(bandwidth) || bandwidth <= 0.0 ||
        !kdml_kernel_is_valid(type, kernel_code)) {
        return R_NaReal;
    }

    if (type == KDML_CONTINUOUS) {
        value = kdml_continuous_kernel(kernel_code,
                                       difference / bandwidth);
        if (data->metric == KDML_DKPS) {
            value /= bandwidth;
        }
        return value;
    }

    if (bandwidth >= kdml_bandwidth_upper(type, kernel_code,
                                          data->level_count[feature])) {
        return R_NaReal;
    }

    if (type == KDML_NOMINAL) {
        if (kernel_code == 0) {
            return difference == 0.0 ? 1.0 : bandwidth;
        }
        if (difference == 0.0) {
            return 1.0 - bandwidth;
        }
        return bandwidth / ((double) data->level_count[feature] - 1.0);
    }

    if (type == KDML_ORDINAL) {
        const int ordinal_distance = (int) floor(distance + 0.5);
        const int ordinal_range = data->level_count[feature] - 1;

        switch (kernel_code) {
        case 0:
            return ordinal_distance == 0
                       ? 1.0 - bandwidth
                       : 0.5 * (1.0 - bandwidth) *
                             pow(bandwidth, (double) ordinal_distance);
        case 1:
            return pow(bandwidth,
                       (double) ordinal_distance * (double) ordinal_distance);
        case 2:
            return ordinal_distance == 0
                       ? bandwidth
                       : (1.0 - bandwidth) *
                             exp(-(double) ordinal_distance * KDML_LOG_TWO);
        case 3:
            if (ordinal_distance > ordinal_range) {
                return 0.0;
            }
            value = lgamma((double) ordinal_range + 1.0) -
                    lgamma((double) ordinal_distance + 1.0) -
                    lgamma((double) (ordinal_range - ordinal_distance) + 1.0) +
                    (double) ordinal_distance * log(bandwidth) +
                    (double) (ordinal_range - ordinal_distance) *
                        log1p(-bandwidth);
            return exp(value);
        case 4:
            return pow(bandwidth, (double) ordinal_distance);
        default:
            return R_NaReal;
        }
    }

    return R_NaReal;
}

double kdml_pair_similarity(const kdml_data *data, int first, int second,
                            const int *kernel_code,
                            const double *bandwidth)
{
    int feature;
    int has_continuous = 0;
    double continuous = data->metric == KDML_DKPS ? 1.0 : 0.0;
    double categorical = 0.0;

    for (feature = 0; feature < data->p; ++feature) {
        const double value = kdml_feature_value(
            data, feature,
            data->x[(size_t) first + (size_t) data->n * (size_t) feature],
            data->x[(size_t) second + (size_t) data->n * (size_t) feature],
            kernel_code[feature], bandwidth[feature]
        );
        if (!R_FINITE(value)) {
            return R_NaReal;
        }
        if (data->type[feature] == KDML_CONTINUOUS) {
            has_continuous = 1;
            if (data->metric == KDML_DKPS) {
                continuous *= value;
            } else {
                continuous += value;
            }
        } else {
            categorical += value;
        }
    }

    if (data->metric == KDML_DKPS && !has_continuous) {
        continuous = 0.0;
    }
    return continuous + categorical;
}

static int kdml_proposed_self_is_finite(const kdml_data *data,
                                        const int *current_kernel,
                                        const double *current_bandwidth,
                                        int proposed_feature,
                                        int proposed_kernel,
                                        double proposed_bandwidth)
{
    int feature;
    int has_continuous = 0;
    double continuous = data->metric == KDML_DKPS ? 1.0 : 0.0;
    double categorical = 0.0;

    for (feature = 0; feature < data->p; ++feature) {
        const int code = feature == proposed_feature
                             ? proposed_kernel
                             : current_kernel[feature];
        const double bandwidth = feature == proposed_feature
                                     ? proposed_bandwidth
                                     : current_bandwidth[feature];
        const double observation = data->x[
            (size_t) data->n * (size_t) feature
        ];
        const double value = kdml_feature_value(
            data, feature, observation, observation, code, bandwidth
        );

        if (!R_FINITE(value)) {
            return 0;
        }
        if (data->type[feature] == KDML_CONTINUOUS) {
            has_continuous = 1;
            if (data->metric == KDML_DKPS) {
                continuous *= value;
                if (!R_FINITE(continuous)) {
                    return 0;
                }
            } else {
                continuous += value;
            }
        } else {
            categorical += value;
        }
    }
    if (data->metric == KDML_DKPS && !has_continuous) {
        continuous = 0.0;
    }
    return R_FINITE(continuous + categorical);
}

double kdml_score_from_rows(const kdml_data *data, const double *row_sums)
{
    int row;
    double total = 0.0;
    const double denominator = (double) data->n - 1.0;

    for (row = 0; row < data->n; ++row) {
        const double leave_one_out = row_sums[row] / denominator;
        if (!R_FINITE(leave_one_out) || leave_one_out <= 0.0) {
            return R_NegInf;
        }
        total += log(leave_one_out);
    }
    total /= (double) data->n;
    return R_FINITE(total) ? total : R_NegInf;
}

double kdml_score_full(const kdml_data *data, const int *kernel_code,
                       const double *bandwidth, double *similarity)
{
    int first;
    int second;
    double *row_sums = (double *) R_alloc((size_t) data->n, sizeof(double));

    memset(row_sums, 0, (size_t) data->n * sizeof(double));
    if (similarity != NULL) {
        const size_t cells = (size_t) data->n * (size_t) data->n;
        memset(similarity, 0, cells * sizeof(double));
    }
    if (!R_FINITE(kdml_pair_similarity(
            data, 0, 0, kernel_code, bandwidth))) {
        return R_NegInf;
    }
    if (similarity != NULL) {
        for (first = 0; first < data->n; ++first) {
            similarity[(size_t) first + (size_t) data->n * (size_t) first] = kdml_pair_similarity(
                data, first, first, kernel_code, bandwidth
            );
        }
    }

    for (first = 0; first < data->n - 1; ++first) {
        for (second = first + 1; second < data->n; ++second) {
            const double value = kdml_pair_similarity(
                data, first, second, kernel_code, bandwidth
            );
            if (!R_FINITE(value)) {
                return R_NegInf;
            }
            row_sums[first] += value;
            row_sums[second] += value;
            if (similarity != NULL) {
                similarity[(size_t) first +
                           (size_t) data->n * (size_t) second] = value;
                similarity[(size_t) second +
                           (size_t) data->n * (size_t) first] = value;
            }
        }
    }
    return kdml_score_from_rows(data, row_sums);
}

int kdml_cache_build(const kdml_data *data, const int *kernel_code,
                     const double *bandwidth, kdml_cache *cache,
                     double *score, int check_interrupt)
{
    int first;
    int second;
    size_t pair = 0;

    if (!R_FINITE(kdml_pair_similarity(
            data, 0, 0, kernel_code, bandwidth))) {
        return 0;
    }
    memset(cache->row_sums, 0, (size_t) data->n * sizeof(double));
    for (first = 0; first < data->n - 1; ++first) {
        if (check_interrupt && (first & 127) == 0) {
            R_CheckUserInterrupt();
        }
        for (second = first + 1; second < data->n; ++second, ++pair) {
            int feature;
            int zero_count = 0;
            double nonzero_product = 1.0;
            double continuous_sum = 0.0;
            double categorical_sum = 0.0;
            double similarity;

            for (feature = 0; feature < data->p; ++feature) {
                const double value = kdml_feature_value(
                    data, feature,
                    data->x[(size_t) first +
                            (size_t) data->n * (size_t) feature],
                    data->x[(size_t) second +
                            (size_t) data->n * (size_t) feature],
                    kernel_code[feature], bandwidth[feature]
                );
                if (!R_FINITE(value)) {
                    return 0;
                }
                if (data->type[feature] == KDML_CONTINUOUS) {
                    if (data->metric == KDML_DKPS) {
                        if (value == 0.0) {
                            ++zero_count;
                        } else {
                            nonzero_product *= value;
                        }
                    } else {
                        continuous_sum += value;
                    }
                } else {
                    categorical_sum += value;
                }
            }

            if (data->metric == KDML_DKPS) {
                cache->continuous_nonzero_product[pair] = nonzero_product;
                cache->continuous_zero_count[pair] = zero_count;
                continuous_sum = data->continuous_count == 0 || zero_count > 0
                                     ? 0.0
                                     : nonzero_product;
            }
            similarity = continuous_sum + categorical_sum;
            if (!R_FINITE(similarity)) {
                return 0;
            }
            cache->row_sums[first] += similarity;
            cache->row_sums[second] += similarity;
        }
    }

    *score = kdml_score_from_rows(data, cache->row_sums);
    return R_FINITE(*score);
}

double kdml_cache_proposal(const kdml_data *data,
                           const int *current_kernel,
                           const double *current_bandwidth,
                           const kdml_cache *cache,
                           int feature, int proposed_kernel,
                           double proposed_bandwidth,
                           kdml_proposal_workspace *workspace,
                           int check_interrupt)
{
    int first;
    int second;
    size_t pair = 0;

    if (!kdml_proposed_self_is_finite(
            data, current_kernel, current_bandwidth,
            feature, proposed_kernel, proposed_bandwidth)) {
        return R_NegInf;
    }
    memcpy(workspace->row_sums, cache->row_sums,
           (size_t) data->n * sizeof(double));
    for (first = 0; first < data->n - 1; ++first) {
        if (check_interrupt && (first & 127) == 0) {
            R_CheckUserInterrupt();
        }
        for (second = first + 1; second < data->n; ++second, ++pair) {
            const double first_value = data->x[
                (size_t) first + (size_t) data->n * (size_t) feature
            ];
            const double second_value = data->x[
                (size_t) second + (size_t) data->n * (size_t) feature
            ];
            const double old_value = kdml_feature_value(
                data, feature, first_value, second_value,
                current_kernel[feature], current_bandwidth[feature]
            );
            const double new_value = kdml_feature_value(
                data, feature, first_value, second_value,
                proposed_kernel, proposed_bandwidth
            );
            double difference;

            if (!R_FINITE(old_value) || !R_FINITE(new_value)) {
                return R_NegInf;
            }

            if (data->metric == KDML_DKPS &&
                data->type[feature] == KDML_CONTINUOUS) {
                int other;
                int new_zero_count = 0;
                double new_product = 1.0;
                double old_continuous;
                double new_continuous;

                /* Recompute the proposed product in feature order.  Updating
                 * it by repeated division and multiplication accumulates
                 * enough rounding error over a long chain to affect MSCV and
                 * Metropolis decisions.  This matches a fresh cache build. */
                for (other = 0; other < data->p; ++other) {
                    double other_value;
                    if (data->type[other] != KDML_CONTINUOUS) {
                        continue;
                    }
                    if (other == feature) {
                        other_value = new_value;
                    } else {
                        other_value = kdml_feature_value(
                            data, other,
                            data->x[(size_t) first +
                                    (size_t) data->n * (size_t) other],
                            data->x[(size_t) second +
                                    (size_t) data->n * (size_t) other],
                            current_kernel[other], current_bandwidth[other]
                        );
                    }
                    if (!R_FINITE(other_value)) {
                        return R_NegInf;
                    }
                    if (other_value == 0.0) {
                        ++new_zero_count;
                    } else {
                        new_product *= other_value;
                    }
                }
                old_continuous =
                    data->continuous_count == 0 ||
                    cache->continuous_zero_count[pair] > 0
                        ? 0.0
                        : cache->continuous_nonzero_product[pair];
                new_continuous =
                    data->continuous_count == 0 || new_zero_count > 0
                        ? 0.0
                        : new_product;
                workspace->continuous_nonzero_product[pair] = new_product;
                workspace->continuous_zero_count[pair] = new_zero_count;
                difference = new_continuous - old_continuous;
            } else {
                difference = new_value - old_value;
            }

            workspace->row_sums[first] += difference;
            workspace->row_sums[second] += difference;
        }
    }

    return kdml_score_from_rows(data, workspace->row_sums);
}

void kdml_cache_accept(const kdml_data *data, int feature,
                       kdml_cache *cache,
                       const kdml_proposal_workspace *workspace)
{
    memcpy(cache->row_sums, workspace->row_sums,
           (size_t) data->n * sizeof(double));
    if (data->metric == KDML_DKPS &&
        data->type[feature] == KDML_CONTINUOUS) {
        memcpy(cache->continuous_nonzero_product,
               workspace->continuous_nonzero_product,
               data->pair_count * sizeof(double));
        memcpy(cache->continuous_zero_count,
               workspace->continuous_zero_count,
               data->pair_count * sizeof(int));
    }
}

double kdml_log_prior(const kdml_data *data, const double *theta,
                      const double *prior_mean, const double *prior_sd)
{
    int feature;
    double answer = 0.0;

    for (feature = 0; feature < data->p; ++feature) {
        const double standardized =
            (theta[feature] - prior_mean[feature]) / prior_sd[feature];
        const double contribution = -log(prior_sd[feature]) -
                                    0.5 * KDML_LOG_TWO_PI -
                                    0.5 * standardized * standardized;
        if (!R_FINITE(contribution)) {
            return R_NegInf;
        }
        answer += contribution;
    }
    return R_FINITE(answer) ? answer : R_NegInf;
}

double kdml_log_target(double score, double beta, double log_prior,
                       double log_model_prior)
{
    double answer;
    if (!R_FINITE(score) || !R_FINITE(log_prior) ||
        !R_FINITE(log_model_prior)) {
        return R_NegInf;
    }
    answer = beta * score + log_prior + log_model_prior;
    return R_FINITE(answer) ? answer : R_NegInf;
}
