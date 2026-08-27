#ifndef KDML_H
#define KDML_H

/*
 * Use only the prefixed R API.  Besides avoiding namespace pollution, this
 * keeps R's legacy `match` macro from rewriting OpenMP declare-variant
 * pragmas in recent omp.h headers.
 */
#ifndef R_NO_REMAP
# define R_NO_REMAP
#endif

#include <R.h>
#include <Rinternals.h>

#include <stddef.h>
#include <stdint.h>

enum {
    KDML_CONTINUOUS = 0,
    KDML_NOMINAL = 1,
    KDML_ORDINAL = 2
};

enum {
    KDML_DKPS = 0,
    KDML_DKSS = 1
};

enum {
    KDML_BACKEND_CPU = 0,
    KDML_BACKEND_CUDA = 1,
    KDML_BACKEND_AUTO = 2
};

#define KDML_BACKEND_ERROR_SIZE 512

typedef struct {
    int n;
    int p;
    int metric;
    const double *x;
    const int *type;
    const int *level_count;
    const double *scale;
    int continuous_count;
    size_t pair_count;
} kdml_data;

typedef struct {
    double *row_sums;
    double *continuous_nonzero_product;
    int *continuous_zero_count;
} kdml_cache;

typedef struct {
    double *row_sums;
    double *continuous_nonzero_product;
    int *continuous_zero_count;
} kdml_proposal_workspace;

typedef struct kdml_score_backend kdml_score_backend;

typedef struct {
    int (*available)(void);
    const char *(*info)(void);
    int (*build)(kdml_score_backend *backend,
                 const int *kernel_code, const double *bandwidth,
                 double *score);
    int (*propose)(kdml_score_backend *backend,
                   const int *current_kernel,
                   const double *current_bandwidth,
                   int feature, int proposed_kernel,
                   double proposed_bandwidth, double *score);
    int (*accept)(kdml_score_backend *backend,
                  int feature, int proposed_kernel,
                  double proposed_bandwidth);
    void (*destroy)(kdml_score_backend *backend);
} kdml_score_backend_ops;

struct kdml_score_backend {
    const char *name;
    const kdml_score_backend_ops *ops;
    const kdml_data *data;
    void *context;
    char error[KDML_BACKEND_ERROR_SIZE];
};

typedef struct {
    int chains;
    int warmup;
    int draws;
    int thin;
    int progress;
    int progress_in_place;
    int chain_threads;
    double beta;
    int adapt;
    double target_accept;
    const double *proposal_sd;
    const double *prior_mean;
    const double *prior_sd;
    const int *candidate_codes;
    const int *candidate_offsets;
    const int *initial_kernel;
    const double *initial_theta;
} kdml_mcmc_config;

typedef struct {
    double *theta;
    double *bandwidth;
    int *kernel_code;
    double *score;
    double *log_target;
    double *final_proposal_sd;
    double *kernel_attempts;
    double *kernel_accepts;
    double *bandwidth_attempts;
    double *bandwidth_accepts;
    double *sampling_kernel_attempts;
    double *sampling_kernel_accepts;
    double *sampling_bandwidth_attempts;
    double *sampling_bandwidth_accepts;
    int *invalid_proposals;
} kdml_mcmc_output;

typedef struct {
    uint64_t state;
    int has_spare;
    double spare;
} kdml_rng;

typedef struct {
    int enabled;
    int active;
    int in_place;
    int aggregate;
    int chain;
    int chains;
    long long total;
    long long completed;
    long long stride;
    double start_seconds;
    double last_emit_seconds;
    int last_width;
} kdml_progress;

int kdml_kernel_is_valid(int type, int code);
double kdml_bandwidth_upper(int type, int code, int level_count);
int kdml_theta_to_bandwidth(const kdml_data *data, int feature,
                            int kernel_code, double theta,
                            double *bandwidth);
double kdml_feature_value(const kdml_data *data, int feature,
                          double first, double second,
                          int kernel_code, double bandwidth);
double kdml_pair_similarity(const kdml_data *data, int first, int second,
                            const int *kernel_code,
                            const double *bandwidth);
double kdml_score_full(const kdml_data *data, const int *kernel_code,
                       const double *bandwidth, double *similarity);
double kdml_score_from_rows(const kdml_data *data, const double *row_sums);

int kdml_cache_build(const kdml_data *data, const int *kernel_code,
                     const double *bandwidth, kdml_cache *cache,
                     double *score, int check_interrupt);
double kdml_cache_proposal(const kdml_data *data,
                           const int *current_kernel,
                           const double *current_bandwidth,
                           const kdml_cache *cache,
                           int feature, int proposed_kernel,
                           double proposed_bandwidth,
                           kdml_proposal_workspace *workspace,
                           int check_interrupt);
void kdml_cache_accept(const kdml_data *data, int feature,
                       kdml_cache *cache,
                       const kdml_proposal_workspace *workspace);

double kdml_log_prior(const kdml_data *data, const double *theta,
                      const double *prior_mean, const double *prior_sd);
double kdml_log_target(double score, double beta, double log_prior,
                       double log_model_prior);

void kdml_progress_initialize(kdml_progress *progress, int enabled,
                              int in_place, int chains, long long total);
void kdml_progress_begin(kdml_progress *progress, int chain);
void kdml_progress_tick(kdml_progress *progress, long long completed);
void kdml_progress_end(kdml_progress *progress, const char *status);
void kdml_progress_group_begin(kdml_progress *progress,
                               const int *completed);
void kdml_progress_group_tick(kdml_progress *progress,
                              const int *completed);
void kdml_progress_group_end(kdml_progress *progress,
                             const int *completed, const char *status);

int kdml_run_mcmc(const kdml_data *data, const kdml_mcmc_config *config,
                  kdml_mcmc_output *output,
                  kdml_score_backend *backend, kdml_progress *progress,
                  int *state_kernel, double *state_theta,
                  double *state_bandwidth, double *state_proposal_sd,
                  int *adapt_count, kdml_rng *rng, int *error_chain);

void kdml_initialize_mcmc_output(const kdml_data *data,
                                 const kdml_mcmc_config *config,
                                 kdml_mcmc_output *output);
int kdml_run_chain(const kdml_data *data,
                   const kdml_mcmc_config *config,
                   kdml_mcmc_output *output,
                   kdml_score_backend *backend, kdml_progress *progress,
                   int chain, int *state_kernel, double *state_theta,
                   double *state_bandwidth, double *state_proposal_sd,
                   int *adapt_count, kdml_rng *rng,
                   volatile int *completed_sweeps, volatile int *cancel,
                   int check_interrupt);

/* Registered .Call entry points. */
SEXP kdml_mcmc_call(SEXP spec);
SEXP kdml_score_call(SEXP spec);
SEXP kdml_distance_call(SEXP spec);
SEXP kdml_cuda_info_call(void);

#endif
