#include "kdml.h"
#include "kdml_cuda_api.h"

#include <R_ext/Random.h>
#include <R_ext/Utils.h>

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
# include <omp.h>
#endif

#ifdef _WIN32
# include <windows.h>
#else
# include <sys/select.h>
#endif

static SEXP kdml_list_get(SEXP list, const char *name)
{
    SEXP names;
    R_xlen_t index;

    if (TYPEOF(list) != VECSXP) {
        Rf_error("Internal KDML specification must be a list.");
    }
    names = Rf_getAttrib(list, R_NamesSymbol);
    if (TYPEOF(names) != STRSXP) {
        Rf_error("Internal KDML specification must be named.");
    }
    for (index = 0; index < XLENGTH(list); ++index) {
        if (strcmp(CHAR(STRING_ELT(names, index)), name) == 0) {
            return VECTOR_ELT(list, index);
        }
    }
    Rf_error("Internal KDML specification is missing `%s`.", name);
    return R_NilValue;
}

static SEXP kdml_list_get_optional(SEXP list, const char *name)
{
    SEXP names;
    R_xlen_t index;

    if (TYPEOF(list) != VECSXP) {
        Rf_error("Internal KDML specification must be a list.");
    }
    names = Rf_getAttrib(list, R_NamesSymbol);
    if (TYPEOF(names) != STRSXP) {
        Rf_error("Internal KDML specification must be named.");
    }
    for (index = 0; index < XLENGTH(list); ++index) {
        if (strcmp(CHAR(STRING_ELT(names, index)), name) == 0) {
            return VECTOR_ELT(list, index);
        }
    }
    return R_NilValue;
}

static int kdml_scalar_integer(SEXP value, const char *name, int minimum)
{
    int answer;
    if (TYPEOF(value) != INTSXP || XLENGTH(value) != 1) {
        Rf_error("Internal `%s` must be a scalar integer.", name);
    }
    answer = INTEGER(value)[0];
    if (answer == NA_INTEGER || answer < minimum) {
        Rf_error("Internal `%s` is outside its supported range.", name);
    }
    return answer;
}

static double kdml_scalar_double(SEXP value, const char *name)
{
    double answer;
    if (TYPEOF(value) != REALSXP || XLENGTH(value) != 1) {
        Rf_error("Internal `%s` must be a scalar double.", name);
    }
    answer = REAL(value)[0];
    if (!R_FINITE(answer)) {
        Rf_error("Internal `%s` must be finite.", name);
    }
    return answer;
}

static int kdml_scalar_logical(SEXP value, const char *name)
{
    int answer;
    if (TYPEOF(value) != LGLSXP || XLENGTH(value) != 1) {
        Rf_error("Internal `%s` must be a scalar logical.", name);
    }
    answer = LOGICAL(value)[0];
    if (answer == NA_LOGICAL) {
        Rf_error("Internal `%s` cannot be missing.", name);
    }
    return answer != 0;
}

static void kdml_setup_data(SEXP spec, kdml_data *data,
                            int require_scale, int minimum_rows)
{
    SEXP x = kdml_list_get(spec, "x");
    SEXP dimensions;
    SEXP type = kdml_list_get(spec, "type");
    SEXP level_count = kdml_list_get(spec, "level_count");
    SEXP scale = R_NilValue;
    int feature;
    R_xlen_t cell;
    size_t n_size;

    if (TYPEOF(x) != REALSXP) {
        Rf_error("Internal `x` must be a double matrix.");
    }
    dimensions = Rf_getAttrib(x, R_DimSymbol);
    if (TYPEOF(dimensions) != INTSXP || XLENGTH(dimensions) != 2) {
        Rf_error("Internal `x` must be a matrix.");
    }
    data->n = INTEGER(dimensions)[0];
    data->p = INTEGER(dimensions)[1];
    if (data->n < minimum_rows || data->p < 1) {
        Rf_error("Internal KDML data dimensions are invalid.");
    }
    if ((R_xlen_t) data->n > R_XLEN_T_MAX / data->p ||
        XLENGTH(x) != (R_xlen_t) data->n * data->p) {
        Rf_error("Internal KDML matrix dimensions are inconsistent.");
    }
    if (TYPEOF(type) != INTSXP || XLENGTH(type) != data->p ||
        TYPEOF(level_count) != INTSXP || XLENGTH(level_count) != data->p) {
        Rf_error("Internal feature metadata have invalid dimensions.");
    }
    if (require_scale) {
        scale = kdml_list_get(spec, "scale");
        if (TYPEOF(scale) != REALSXP || XLENGTH(scale) != data->p) {
            Rf_error("Internal continuous scales have invalid dimensions.");
        }
    }

    data->x = REAL(x);
    data->type = INTEGER(type);
    data->level_count = INTEGER(level_count);
    data->scale = require_scale ? REAL(scale) : NULL;
    data->continuous_count = 0;
    data->metric = kdml_scalar_integer(
        kdml_list_get(spec, "metric"), "metric", 0
    );
    if (data->metric != KDML_DKPS && data->metric != KDML_DKSS) {
        Rf_error("Internal metric code is invalid.");
    }

    for (feature = 0; feature < data->p; ++feature) {
        const int feature_type = data->type[feature];
        if (feature_type < KDML_CONTINUOUS || feature_type > KDML_ORDINAL) {
            Rf_error("Internal feature type code is invalid.");
        }
        if (feature_type == KDML_CONTINUOUS) {
            ++data->continuous_count;
            if (require_scale &&
                (!R_FINITE(data->scale[feature]) ||
                 data->scale[feature] <= 0.0)) {
                Rf_error("Internal continuous scale is invalid.");
            }
        } else if (data->level_count[feature] < 2) {
            Rf_error("Categorical features need at least two observed levels.");
        }
    }

    for (cell = 0; cell < XLENGTH(x); ++cell) {
        if (!R_FINITE(data->x[cell])) {
            Rf_error("Internal KDML data contain a non-finite value.");
        }
    }
    for (feature = 0; feature < data->p; ++feature) {
        int row;
        if (data->type[feature] == KDML_CONTINUOUS) {
            continue;
        }
        for (row = 0; row < data->n; ++row) {
            const double value = data->x[
                (size_t) row + (size_t) data->n * (size_t) feature
            ];
            if (value < 0.0 ||
                value > (double) data->level_count[feature] - 1.0 ||
                value != floor(value)) {
                Rf_error("Internal categorical encoding is invalid.");
            }
        }
    }

    n_size = (size_t) data->n;
    if (n_size > SIZE_MAX / n_size) {
        Rf_error("The requested square matrix is too large.");
    }
    if (data->n > 1 && n_size > SIZE_MAX / (n_size - 1U)) {
        Rf_error("The number of observation pairs is too large.");
    }
    data->pair_count = data->n > 1
                           ? n_size * (n_size - 1U) / 2U
                           : 0U;
}

static void kdml_validate_state(const kdml_data *data,
                                SEXP kernel, SEXP bandwidth)
{
    int feature;
    if (TYPEOF(kernel) != INTSXP || XLENGTH(kernel) != data->p ||
        TYPEOF(bandwidth) != REALSXP || XLENGTH(bandwidth) != data->p) {
        Rf_error("Internal kernel state has invalid dimensions.");
    }
    for (feature = 0; feature < data->p; ++feature) {
        const int code = INTEGER(kernel)[feature];
        const double value = REAL(bandwidth)[feature];
        const double upper = kdml_bandwidth_upper(
            data->type[feature], code, data->level_count[feature]
        );
        if (!kdml_kernel_is_valid(data->type[feature], code)) {
            Rf_error("Internal kernel code is incompatible with its feature.");
        }
        if (!R_FINITE(value) || value <= 0.0 || value >= upper) {
            Rf_error("Internal bandwidth is outside its open MCMC support.");
        }
    }
}

static void kdml_set_names(SEXP object, const char **names, int count)
{
    int index;
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, count));
    for (index = 0; index < count; ++index) {
        SET_STRING_ELT(result_names, index, Rf_mkChar(names[index]));
    }
    Rf_setAttrib(object, R_NamesSymbol, result_names);
    UNPROTECT(1);
}

static void kdml_set_list_element(SEXP list, int index, SEXP value)
{
    SET_VECTOR_ELT(list, index, value);
}

typedef struct {
    const kdml_data *data;
    int check_interrupt;
    kdml_cache cache;
    kdml_proposal_workspace workspace;
} kdml_cpu_backend_context;

static int kdml_cpu_available(void)
{
    return 1;
}

static const char *kdml_cpu_info(void)
{
    return "portable C";
}

static int kdml_cpu_build(kdml_score_backend *backend,
                          const int *kernel_code,
                          const double *bandwidth, double *score)
{
    kdml_cpu_backend_context *context =
        (kdml_cpu_backend_context *) backend->context;
    return kdml_cache_build(
        context->data, kernel_code, bandwidth, &context->cache, score,
        context->check_interrupt
    );
}

static int kdml_cpu_propose(kdml_score_backend *backend,
                            const int *current_kernel,
                            const double *current_bandwidth,
                            int feature, int proposed_kernel,
                            double proposed_bandwidth, double *score)
{
    kdml_cpu_backend_context *context =
        (kdml_cpu_backend_context *) backend->context;
    *score = kdml_cache_proposal(
        context->data, current_kernel, current_bandwidth,
        &context->cache, feature, proposed_kernel,
        proposed_bandwidth, &context->workspace,
        context->check_interrupt
    );
    return 1;
}

static int kdml_cpu_accept(kdml_score_backend *backend,
                           int feature, int proposed_kernel,
                           double proposed_bandwidth)
{
    kdml_cpu_backend_context *context =
        (kdml_cpu_backend_context *) backend->context;
    (void) proposed_kernel;
    (void) proposed_bandwidth;
    kdml_cache_accept(
        context->data, feature, &context->cache, &context->workspace
    );
    return 1;
}

static void kdml_cpu_destroy(kdml_score_backend *backend)
{
    backend->context = NULL;
}

static int kdml_cuda_backend_available(void)
{
    return kdml_cuda_available() != 0;
}

static const char *kdml_cuda_backend_info(void)
{
    return kdml_cuda_compile_version_string();
}

static void kdml_copy_cuda_error(kdml_score_backend *backend,
                                 const char *operation, int status)
{
    const char *detail = kdml_cuda_last_error();
    const char *status_text = kdml_cuda_status_string(status);
    if (detail != NULL && detail[0] != '\0') {
        snprintf(
            backend->error, sizeof(backend->error),
            "CUDA %s failed (%s): %s", operation,
            status_text == NULL ? "unknown status" : status_text, detail
        );
    } else {
        snprintf(
            backend->error, sizeof(backend->error),
            "CUDA %s failed (%s).", operation,
            status_text == NULL ? "unknown status" : status_text
        );
    }
}

static int kdml_cuda_backend_build(kdml_score_backend *backend,
                                   const int *kernel_code,
                                   const double *bandwidth, double *score)
{
    const int status = kdml_cuda_build(
        (kdml_cuda_workspace *) backend->context,
        kernel_code, bandwidth, score
    );
    if (status != KDML_CUDA_SUCCESS) {
        kdml_copy_cuda_error(backend, "build", status);
        return 0;
    }
    return 1;
}

static int kdml_cuda_backend_propose(kdml_score_backend *backend,
                                     const int *current_kernel,
                                     const double *current_bandwidth,
                                     int feature, int proposed_kernel,
                                     double proposed_bandwidth, double *score)
{
    const int status = kdml_cuda_propose(
        (kdml_cuda_workspace *) backend->context,
        feature, proposed_kernel, proposed_bandwidth, score
    );
    (void) current_kernel;
    (void) current_bandwidth;
    if (status != KDML_CUDA_SUCCESS) {
        kdml_copy_cuda_error(backend, "proposal", status);
        return 0;
    }
    return 1;
}

static int kdml_cuda_backend_accept(kdml_score_backend *backend,
                                    int feature, int proposed_kernel,
                                    double proposed_bandwidth)
{
    const int status = kdml_cuda_accept(
        (kdml_cuda_workspace *) backend->context,
        feature, proposed_kernel, proposed_bandwidth
    );
    if (status != KDML_CUDA_SUCCESS) {
        kdml_copy_cuda_error(backend, "accept", status);
        return 0;
    }
    return 1;
}

static void kdml_cuda_backend_destroy(kdml_score_backend *backend)
{
    if (backend->context != NULL) {
        kdml_cuda_destroy((kdml_cuda_workspace *) backend->context);
        backend->context = NULL;
    }
}

static const kdml_score_backend_ops kdml_cpu_backend_ops = {
    kdml_cpu_available,
    kdml_cpu_info,
    kdml_cpu_build,
    kdml_cpu_propose,
    kdml_cpu_accept,
    kdml_cpu_destroy
};

static const kdml_score_backend_ops kdml_cuda_backend_ops = {
    kdml_cuda_backend_available,
    kdml_cuda_backend_info,
    kdml_cuda_backend_build,
    kdml_cuda_backend_propose,
    kdml_cuda_backend_accept,
    kdml_cuda_backend_destroy
};

static void kdml_backend_use_cpu(kdml_score_backend *backend,
                                 const kdml_data *data,
                                 kdml_cpu_backend_context *context,
                                 int check_interrupt)
{
    const size_t pair_allocation =
        data->pair_count > 0 ? data->pair_count : 1U;

    context->data = data;
    context->check_interrupt = check_interrupt != 0;
    context->cache.row_sums = (double *) R_alloc(
        (size_t) data->n, sizeof(double)
    );
    context->cache.continuous_nonzero_product = (double *) R_alloc(
        pair_allocation, sizeof(double)
    );
    context->cache.continuous_zero_count = (int *) R_alloc(
        pair_allocation, sizeof(int)
    );
    context->workspace.row_sums = (double *) R_alloc(
        (size_t) data->n, sizeof(double)
    );
    context->workspace.continuous_nonzero_product = (double *) R_alloc(
        pair_allocation, sizeof(double)
    );
    context->workspace.continuous_zero_count = (int *) R_alloc(
        pair_allocation, sizeof(int)
    );
    backend->name = "cpu";
    backend->ops = &kdml_cpu_backend_ops;
    backend->data = data;
    backend->context = context;
    backend->error[0] = '\0';
}

static int kdml_backend_initialize(kdml_score_backend *backend,
                                   int requested,
                                   const kdml_data *data,
                                   kdml_cpu_backend_context *cpu_context)
{
    kdml_cuda_workspace *workspace = NULL;
    int status;

    memset(backend, 0, sizeof(*backend));
    if (requested == KDML_BACKEND_CPU) {
        kdml_backend_use_cpu(backend, data, cpu_context, 1);
        return 1;
    }

    if (!kdml_cuda_backend_ops.available()) {
        const char *detail = kdml_cuda_last_error();
        if (requested == KDML_BACKEND_AUTO) {
            kdml_backend_use_cpu(backend, data, cpu_context, 1);
            return 1;
        }
        backend->name = "cuda";
        backend->ops = &kdml_cuda_backend_ops;
        backend->data = data;
        snprintf(
            backend->error, sizeof(backend->error),
            "CUDA backend is unavailable: %s",
            detail == NULL || detail[0] == '\0' ?
                "the package was built without a usable CUDA device." :
                detail
        );
        return 0;
    }

    status = kdml_cuda_create(
        data->n, data->p, data->metric, data->x, data->type,
        data->level_count, data->continuous_count, -1, &workspace
    );
    if (status == KDML_CUDA_SUCCESS) {
        backend->name = "cuda";
        backend->ops = &kdml_cuda_backend_ops;
        backend->data = data;
        backend->context = workspace;
        backend->error[0] = '\0';
        return 1;
    }

    backend->name = "cuda";
    backend->ops = &kdml_cuda_backend_ops;
    backend->data = data;
    backend->context = workspace;
    kdml_copy_cuda_error(backend, "initialization", status);
    kdml_cuda_backend_destroy(backend);
    if (requested == KDML_BACKEND_AUTO) {
        kdml_backend_use_cpu(backend, data, cpu_context, 1);
        return 1;
    }
    return 0;
}

static void kdml_backend_destroy(kdml_score_backend *backend)
{
    if (backend->ops != NULL && backend->ops->destroy != NULL) {
        backend->ops->destroy(backend);
    }
}

SEXP kdml_cuda_available_call(void)
{
    return Rf_ScalarLogical(kdml_cuda_available() != 0);
}

SEXP kdml_cuda_info_call(void)
{
    static const char *result_names[] = {
        "compiled", "available", "compile_version",
        "compile_version_string", "runtime_version",
        "device_count", "devices"
    };
    static const char *device_names[] = {
        "device", "name", "total_global_memory",
        "compute_capability_major", "compute_capability_minor",
        "multiprocessor_count", "driver_version", "runtime_version"
    };
    SEXP result;
    SEXP devices;
    int runtime_version = NA_INTEGER;
    int device_count = 0;
    int device;

    if (kdml_cuda_runtime_version(&runtime_version) != KDML_CUDA_SUCCESS) {
        runtime_version = NA_INTEGER;
    }
    if (kdml_cuda_device_count(&device_count) != KDML_CUDA_SUCCESS ||
        device_count < 0) {
        device_count = 0;
    }

    result = PROTECT(Rf_allocVector(VECSXP, 7));
    SET_VECTOR_ELT(result, 0, Rf_ScalarLogical(kdml_cuda_compiled() != 0));
    SET_VECTOR_ELT(result, 1, Rf_ScalarLogical(kdml_cuda_available() != 0));
    SET_VECTOR_ELT(result, 2, Rf_ScalarInteger(kdml_cuda_compile_version()));
    SET_VECTOR_ELT(
        result, 3,
        Rf_mkString(kdml_cuda_compile_version_string() == NULL ?
                    "" : kdml_cuda_compile_version_string())
    );
    SET_VECTOR_ELT(result, 4, Rf_ScalarInteger(runtime_version));
    SET_VECTOR_ELT(result, 5, Rf_ScalarInteger(device_count));

    devices = PROTECT(Rf_allocVector(VECSXP, device_count));
    for (device = 0; device < device_count; ++device) {
        kdml_cuda_device_properties properties;
        SEXP item;
        if (kdml_cuda_device_info(device, &properties) !=
            KDML_CUDA_SUCCESS) {
            SET_VECTOR_ELT(devices, device, R_NilValue);
            continue;
        }
        item = PROTECT(Rf_allocVector(VECSXP, 8));
        SET_VECTOR_ELT(item, 0, Rf_ScalarInteger(properties.device));
        SET_VECTOR_ELT(item, 1, Rf_mkString(properties.name));
        SET_VECTOR_ELT(
            item, 2, Rf_ScalarReal((double) properties.total_global_memory)
        );
        SET_VECTOR_ELT(
            item, 3,
            Rf_ScalarInteger(properties.compute_capability_major)
        );
        SET_VECTOR_ELT(
            item, 4,
            Rf_ScalarInteger(properties.compute_capability_minor)
        );
        SET_VECTOR_ELT(
            item, 5, Rf_ScalarInteger(properties.multiprocessor_count)
        );
        SET_VECTOR_ELT(
            item, 6, Rf_ScalarInteger(properties.driver_version)
        );
        SET_VECTOR_ELT(
            item, 7, Rf_ScalarInteger(properties.runtime_version)
        );
        kdml_set_names(item, device_names, 8);
        SET_VECTOR_ELT(devices, device, item);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(result, 6, devices);
    UNPROTECT(1);
    kdml_set_names(result, result_names, 7);
    UNPROTECT(1);
    return result;
}

typedef struct {
    const kdml_data *data;
    const int *kernel_code;
    const double *bandwidth;
    int requested_backend;
    kdml_cpu_backend_context *cpu_context;
    kdml_score_backend backend;
    double *score;
    int status;
} kdml_score_run_context;

static SEXP kdml_run_score_unwind_body(void *pointer)
{
    kdml_score_run_context *context =
        (kdml_score_run_context *) pointer;
    *context->score = R_NegInf;
    if (!kdml_backend_initialize(
            &context->backend, context->requested_backend,
            context->data, context->cpu_context)) {
        context->status = 1;
        return R_NilValue;
    }
    if (!context->backend.ops->build(
            &context->backend, context->kernel_code,
            context->bandwidth, context->score) &&
        context->backend.error[0] != '\0') {
        context->status = 1;
    }
    return R_NilValue;
}

static void kdml_destroy_score_backend(void *pointer, Rboolean jump)
{
    kdml_score_run_context *context =
        (kdml_score_run_context *) pointer;
    (void) jump;
    kdml_backend_destroy(&context->backend);
}

SEXP kdml_score_call(SEXP spec)
{
    static const char *result_names[] = {"score", "similarity"};
    kdml_data data;
    SEXP kernel;
    SEXP bandwidth;
    SEXP return_similarity;
    SEXP result;
    SEXP score;
    SEXP similarity = R_NilValue;
    double *similarity_pointer = NULL;
    int requested_backend;

    kdml_setup_data(spec, &data, 0, 2);
    kernel = kdml_list_get(spec, "kernel");
    bandwidth = kdml_list_get(spec, "bandwidth");
    return_similarity = kdml_list_get(spec, "return_similarity");
    requested_backend = kdml_scalar_integer(
        kdml_list_get(spec, "backend"), "backend", KDML_BACKEND_CPU
    );
    if (requested_backend > KDML_BACKEND_AUTO) {
        Rf_error("Internal scoring backend code is invalid.");
    }
    kdml_validate_state(&data, kernel, bandwidth);
    if (TYPEOF(return_similarity) != LGLSXP ||
        XLENGTH(return_similarity) != 1 ||
        LOGICAL(return_similarity)[0] == NA_LOGICAL) {
        Rf_error("Internal `return_similarity` must be logical.");
    }

    result = PROTECT(Rf_allocVector(VECSXP, 2));
    score = PROTECT(Rf_allocVector(REALSXP, 1));
    kdml_set_list_element(result, 0, score);
    UNPROTECT(1);
    if (LOGICAL(return_similarity)[0]) {
        if (requested_backend == KDML_BACKEND_CUDA) {
            UNPROTECT(1);
            Rf_error("CUDA score evaluation cannot return a similarity matrix.");
        }
        similarity = PROTECT(Rf_allocMatrix(REALSXP, data.n, data.n));
        similarity_pointer = REAL(similarity);
        kdml_set_list_element(result, 1, similarity);
        UNPROTECT(1);
    } else {
        kdml_set_list_element(result, 1, R_NilValue);
    }
    if (similarity_pointer != NULL) {
        REAL(VECTOR_ELT(result, 0))[0] = kdml_score_full(
            &data, INTEGER(kernel), REAL(bandwidth), similarity_pointer
        );
    } else {
        kdml_cpu_backend_context cpu_context;
        kdml_score_run_context run_context;
        SEXP unwind_continuation;

        memset(&cpu_context, 0, sizeof(cpu_context));
        memset(&run_context, 0, sizeof(run_context));
        run_context.data = &data;
        run_context.kernel_code = INTEGER(kernel);
        run_context.bandwidth = REAL(bandwidth);
        run_context.requested_backend = requested_backend;
        run_context.cpu_context = &cpu_context;
        run_context.score = REAL(VECTOR_ELT(result, 0));

        unwind_continuation = PROTECT(R_MakeUnwindCont());
        R_UnwindProtect(
            kdml_run_score_unwind_body, &run_context,
            kdml_destroy_score_backend, &run_context,
            unwind_continuation
        );
        UNPROTECT(1);
        if (run_context.status != 0) {
            const char *message = run_context.backend.error[0] == '\0'
                ? "Scoring backend failed."
                : run_context.backend.error;
            UNPROTECT(1);
            Rf_error("%s", message);
        }
    }
    kdml_set_names(result, result_names, 2);
    UNPROTECT(1);
    return result;
}

SEXP kdml_distance_call(SEXP spec)
{
    kdml_data data;
    SEXP kernel;
    SEXP bandwidth;
    SEXP result;
    double *self;
    int first;
    int second;

    kdml_setup_data(spec, &data, 0, 1);
    kernel = kdml_list_get(spec, "kernel");
    bandwidth = kdml_list_get(spec, "bandwidth");
    kdml_validate_state(&data, kernel, bandwidth);

    result = PROTECT(Rf_allocMatrix(REALSXP, data.n, data.n));
    memset(REAL(result), 0,
           (size_t) data.n * (size_t) data.n * sizeof(double));
    self = (double *) R_alloc((size_t) data.n, sizeof(double));
    for (first = 0; first < data.n; ++first) {
        if ((first & 127) == 0) {
            R_CheckUserInterrupt();
        }
        self[first] = kdml_pair_similarity(
            &data, first, first, INTEGER(kernel), REAL(bandwidth)
        );
        if (!R_FINITE(self[first])) {
            UNPROTECT(1);
            Rf_error("The selected state produced a non-finite self-similarity.");
        }
    }
    for (first = 0; first < data.n - 1; ++first) {
        if ((first & 127) == 0) {
            R_CheckUserInterrupt();
        }
        for (second = first + 1; second < data.n; ++second) {
            const double pair = kdml_pair_similarity(
                &data, first, second, INTEGER(kernel), REAL(bandwidth)
            );
            const double distance = self[first] + self[second] - 2.0 * pair;
            if (!R_FINITE(pair) || !R_FINITE(distance)) {
                UNPROTECT(1);
                Rf_error("The selected state produced a non-finite distance.");
            }
            REAL(result)[(size_t) first +
                         (size_t) data.n * (size_t) second] = distance;
            REAL(result)[(size_t) second +
                         (size_t) data.n * (size_t) first] = distance;
        }
    }
    UNPROTECT(1);
    return result;
}

static void kdml_validate_mcmc_spec(SEXP spec, const kdml_data *data,
                                    kdml_mcmc_config *config)
{
    SEXP candidate_codes = kdml_list_get(spec, "candidate_codes");
    SEXP candidate_offsets = kdml_list_get(spec, "candidate_offsets");
    SEXP initial_kernel = kdml_list_get(spec, "initial_kernel");
    SEXP initial_theta = kdml_list_get(spec, "initial_theta");
    SEXP proposal_sd = kdml_list_get(spec, "proposal_sd");
    SEXP prior_mean = kdml_list_get(spec, "prior_mean");
    SEXP prior_sd = kdml_list_get(spec, "prior_sd");
    SEXP progress = kdml_list_get_optional(spec, "progress");
    SEXP progress_in_place =
        kdml_list_get_optional(spec, "progress_in_place");
    SEXP chain_threads = kdml_list_get_optional(spec, "chain_threads");
    int feature;

    if (TYPEOF(candidate_codes) != INTSXP ||
        TYPEOF(candidate_offsets) != INTSXP ||
        XLENGTH(candidate_offsets) != (R_xlen_t) data->p + 1 ||
        TYPEOF(initial_kernel) != INTSXP ||
        XLENGTH(initial_kernel) != data->p ||
        TYPEOF(initial_theta) != REALSXP ||
        XLENGTH(initial_theta) != data->p ||
        TYPEOF(proposal_sd) != REALSXP ||
        XLENGTH(proposal_sd) != data->p ||
        TYPEOF(prior_mean) != REALSXP ||
        XLENGTH(prior_mean) != data->p ||
        TYPEOF(prior_sd) != REALSXP ||
        XLENGTH(prior_sd) != data->p) {
        Rf_error("Internal MCMC vectors have invalid dimensions.");
    }
    if (INTEGER(candidate_offsets)[0] != 0 ||
        INTEGER(candidate_offsets)[data->p] != XLENGTH(candidate_codes)) {
        Rf_error("Internal candidate offsets are invalid.");
    }

    for (feature = 0; feature < data->p; ++feature) {
        const int begin = INTEGER(candidate_offsets)[feature];
        const int end = INTEGER(candidate_offsets)[feature + 1];
        int position;
        int initial_found = 0;

        if (begin < 0 || end <= begin || end > XLENGTH(candidate_codes)) {
            Rf_error("Every feature needs at least one kernel candidate.");
        }
        for (position = begin; position < end; ++position) {
            int earlier;
            const int code = INTEGER(candidate_codes)[position];
            if (!kdml_kernel_is_valid(data->type[feature], code)) {
                Rf_error("A kernel candidate is incompatible with its feature.");
            }
            for (earlier = begin; earlier < position; ++earlier) {
                if (INTEGER(candidate_codes)[earlier] == code) {
                    Rf_error("Kernel candidate sets cannot contain duplicates.");
                }
            }
            if (code == INTEGER(initial_kernel)[feature]) {
                initial_found = 1;
            }
        }
        if (!initial_found) {
            Rf_error("An initial kernel is absent from its candidate set.");
        }
        if (!R_FINITE(REAL(initial_theta)[feature]) ||
            !R_FINITE(REAL(proposal_sd)[feature]) ||
            REAL(proposal_sd)[feature] <= 0.0 ||
            !R_FINITE(REAL(prior_mean)[feature]) ||
            !R_FINITE(REAL(prior_sd)[feature]) ||
            REAL(prior_sd)[feature] <= 0.0) {
            Rf_error("Internal MCMC parameters must be finite with positive scales.");
        }
    }

    config->candidate_codes = INTEGER(candidate_codes);
    config->candidate_offsets = INTEGER(candidate_offsets);
    config->initial_kernel = INTEGER(initial_kernel);
    config->initial_theta = REAL(initial_theta);
    config->proposal_sd = REAL(proposal_sd);
    config->prior_mean = REAL(prior_mean);
    config->prior_sd = REAL(prior_sd);
    config->chains = kdml_scalar_integer(
        kdml_list_get(spec, "chains"), "chains", 1
    );
    config->warmup = kdml_scalar_integer(
        kdml_list_get(spec, "warmup"), "warmup", 0
    );
    config->draws = kdml_scalar_integer(
        kdml_list_get(spec, "draws"), "draws", 1
    );
    config->thin = kdml_scalar_integer(
        kdml_list_get(spec, "thin"), "thin", 1
    );
    config->progress = progress == R_NilValue ? 0 :
        kdml_scalar_logical(progress, "progress");
    config->progress_in_place = progress_in_place == R_NilValue ? 0 :
        kdml_scalar_logical(progress_in_place, "progress_in_place");
    config->chain_threads = chain_threads == R_NilValue ? 1 :
        kdml_scalar_integer(chain_threads, "chain_threads", 1);
    if (config->chain_threads > config->chains) {
        config->chain_threads = config->chains;
    }
    if (config->draws > (INT_MAX - config->warmup) / config->thin) {
        Rf_error("Internal iteration count overflows the supported range.");
    }
    config->beta = kdml_scalar_double(kdml_list_get(spec, "beta"), "beta");
    if (config->beta <= 0.0) {
        Rf_error("Internal `beta` must be positive.");
    }
    config->adapt = kdml_scalar_logical(
        kdml_list_get(spec, "adapt"), "adapt"
    );
    config->target_accept = kdml_scalar_double(
        kdml_list_get(spec, "target_accept"), "target_accept"
    );
    if (config->target_accept <= 0.0 || config->target_accept >= 1.0) {
        Rf_error("Internal target acceptance probability is invalid.");
    }
}

static SEXP kdml_allocate_matrix_in_list(SEXP result, int index,
                                         SEXPTYPE type, int rows, int columns)
{
    SEXP value = PROTECT(Rf_allocMatrix(type, rows, columns));
    SET_VECTOR_ELT(result, index, value);
    UNPROTECT(1);
    return VECTOR_ELT(result, index);
}

static SEXP kdml_allocate_array_in_list(SEXP result, int index,
                                        SEXPTYPE type,
                                        int first, int second, int third)
{
    SEXP value = PROTECT(Rf_alloc3DArray(type, first, second, third));
    SET_VECTOR_ELT(result, index, value);
    UNPROTECT(1);
    return VECTOR_ELT(result, index);
}

typedef struct {
    const kdml_data *data;
    const kdml_mcmc_config *config;
    kdml_mcmc_output *output;
    int requested_backend;
    kdml_cpu_backend_context *cpu_context;
    kdml_score_backend backend;
    int *state_kernel;
    double *state_theta;
    double *state_bandwidth;
    double *state_proposal_sd;
    int *adapt_count;
    kdml_rng *rng;
    int *error_chain;
    kdml_progress progress;
    int status;
} kdml_run_context;

static SEXP kdml_run_mcmc_unwind_body(void *pointer)
{
    kdml_run_context *context = (kdml_run_context *) pointer;
    if (!kdml_backend_initialize(
            &context->backend, context->requested_backend,
            context->data, context->cpu_context)) {
        context->status = 4;
        return R_NilValue;
    }
    context->status = kdml_run_mcmc(
        context->data, context->config, context->output,
        &context->backend, &context->progress,
        context->state_kernel, context->state_theta,
        context->state_bandwidth, context->state_proposal_sd,
        context->adapt_count, context->rng, context->error_chain
    );
    if (context->status != 0) {
        kdml_progress_end(&context->progress, "failed");
    }
    return R_NilValue;
}

static void kdml_cleanup_mcmc(void *pointer, Rboolean jump)
{
    kdml_run_context *context = (kdml_run_context *) pointer;
    if (jump) {
        kdml_progress_end(&context->progress, "interrupted");
    }
    kdml_backend_destroy(&context->backend);
}

static uint64_t kdml_splitmix64(uint64_t value)
{
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static void kdml_seed_chains(kdml_rng *rng, int chains)
{
    int chain;

    GetRNGstate();
    for (chain = 0; chain < chains; ++chain) {
        const uint64_t high =
            (uint64_t) floor(unif_rand() * 4294967296.0);
        const uint64_t low =
            (uint64_t) floor(unif_rand() * 4294967296.0);
        uint64_t seed = kdml_splitmix64(
            (high << 32) ^ low ^
            (UINT64_C(0x9e3779b97f4a7c15) * (uint64_t) (chain + 1))
        );
        if (seed == 0) {
            seed = UINT64_C(0x2545f4914f6cdd1d);
        }
        rng[chain].state = seed;
        rng[chain].has_spare = 0;
        rng[chain].spare = 0.0;
    }
    PutRNGstate();
}

#ifdef _OPENMP
static void kdml_check_interrupt(void *unused)
{
    (void) unused;
    R_CheckUserInterrupt();
}

static void kdml_monitor_pause(void)
{
#ifdef _WIN32
    Sleep(25);
#else
    struct timeval delay;
    delay.tv_sec = 0;
    delay.tv_usec = 25000;
    (void) select(0, NULL, NULL, NULL, &delay);
#endif
}

static int kdml_atomic_read(volatile int *value)
{
    int answer;
# pragma omp atomic read
    answer = *value;
    return answer;
}

static void kdml_atomic_write(volatile int *value, int answer)
{
# pragma omp atomic write
    *value = answer;
}

static int kdml_run_parallel_cpu(const kdml_data *data,
                                 const kdml_mcmc_config *config,
                                 kdml_mcmc_output *output,
                                 kdml_rng *rng, kdml_progress *progress,
                                 int *error_chain, int *interrupted)
{
    const int workers =
        config->chain_threads < config->chains ?
            config->chain_threads : config->chains;
    kdml_cpu_backend_context *contexts =
        (kdml_cpu_backend_context *) R_alloc(
            (size_t) workers, sizeof(kdml_cpu_backend_context)
        );
    kdml_score_backend *backends =
        (kdml_score_backend *) R_alloc(
            (size_t) workers, sizeof(kdml_score_backend)
        );
    int *state_kernel = (int *) R_alloc(
        (size_t) data->p * (size_t) workers, sizeof(int)
    );
    double *state_theta = (double *) R_alloc(
        (size_t) data->p * (size_t) workers, sizeof(double)
    );
    double *state_bandwidth = (double *) R_alloc(
        (size_t) data->p * (size_t) workers, sizeof(double)
    );
    double *state_proposal_sd = (double *) R_alloc(
        (size_t) data->p * (size_t) workers, sizeof(double)
    );
    int *adapt_count = (int *) R_alloc(
        (size_t) data->p * (size_t) workers, sizeof(int)
    );
    int *statuses = (int *) R_alloc(
        (size_t) config->chains, sizeof(int)
    );
    volatile int *completed = (volatile int *) R_alloc(
        (size_t) config->chains, sizeof(int)
    );
    int *progress_snapshot = (int *) R_alloc(
        (size_t) config->chains, sizeof(int)
    );
    volatile int cancel = 0;
    volatile int next_chain = 0;
    volatile int workers_done = 0;
    int active_workers = workers;
    int dynamic_threads;
    int chain;
    int status = 0;

    memset(contexts, 0,
           (size_t) workers * sizeof(kdml_cpu_backend_context));
    memset(backends, 0,
           (size_t) workers * sizeof(kdml_score_backend));
    memset(statuses, 0, (size_t) config->chains * sizeof(int));
    for (chain = 0; chain < config->chains; ++chain) {
        completed[chain] = 0;
        progress_snapshot[chain] = 0;
    }
    for (chain = 0; chain < workers; ++chain) {
        kdml_backend_use_cpu(
            &backends[chain], data, &contexts[chain], 0
        );
    }
    kdml_initialize_mcmc_output(data, config, output);
    kdml_progress_group_begin(progress, progress_snapshot);
    *interrupted = 0;
    dynamic_threads = omp_get_dynamic();
    omp_set_dynamic(0);

# pragma omp parallel num_threads(workers + 1) \
    shared(next_chain, workers_done, cancel, completed, statuses, active_workers)
    {
        if (omp_get_thread_num() == 0) {
            active_workers = omp_get_num_threads() - 1;
            while (kdml_atomic_read(&workers_done) < active_workers) {
                int index;
                for (index = 0; index < config->chains; ++index) {
                    progress_snapshot[index] =
                        kdml_atomic_read(&completed[index]);
                }
                kdml_progress_group_tick(progress, progress_snapshot);
                if (!R_ToplevelExec(kdml_check_interrupt, NULL)) {
                    *interrupted = 1;
                    kdml_atomic_write(&cancel, 1);
                }
                kdml_monitor_pause();
            }
        } else {
            const int worker = omp_get_thread_num() - 1;
            for (;;) {
                int assigned;
# pragma omp atomic capture
                {
                    assigned = next_chain;
                    next_chain++;
                }
                if (assigned >= config->chains ||
                    kdml_atomic_read(&cancel)) {
                    break;
                }
                statuses[assigned] = kdml_run_chain(
                    data, config, output, &backends[worker], NULL,
                    assigned,
                    state_kernel + (size_t) worker * data->p,
                    state_theta + (size_t) worker * data->p,
                    state_bandwidth + (size_t) worker * data->p,
                    state_proposal_sd + (size_t) worker * data->p,
                    adapt_count + (size_t) worker * data->p,
                    &rng[assigned], &completed[assigned], &cancel, 0
                );
                if (statuses[assigned] != 0 &&
                    statuses[assigned] != 5) {
                    kdml_atomic_write(&cancel, 1);
                }
            }
# pragma omp atomic update
            workers_done++;
        }
    }
    omp_set_dynamic(dynamic_threads);

    if (progress != NULL) {
        for (chain = 0; chain < config->chains; ++chain) {
            progress_snapshot[chain] = completed[chain];
        }
        kdml_progress_group_end(
            progress, progress_snapshot,
            *interrupted ? "interrupted" :
                (cancel ? "failed" : "done")
        );
    }
    for (chain = 0; chain < workers; ++chain) {
        kdml_backend_destroy(&backends[chain]);
    }
    for (chain = 0; chain < config->chains; ++chain) {
        if (statuses[chain] != 0 && statuses[chain] != 5 &&
            status == 0) {
            status = statuses[chain];
            *error_chain = chain;
        }
    }
    if (*interrupted) {
        return 5;
    }
    if (active_workers == 0) {
        return 6;
    }
    return status;
}
#endif

SEXP kdml_mcmc_call(SEXP spec)
{
    static const char *result_names[] = {
        "theta", "bandwidth", "kernel_code", "score", "log_target",
        "proposal_sd", "kernel_attempts", "kernel_accepts",
        "bandwidth_attempts", "bandwidth_accepts",
        "sampling_kernel_attempts", "sampling_kernel_accepts",
        "sampling_bandwidth_attempts", "sampling_bandwidth_accepts",
        "invalid_proposals", "backend"
    };
    kdml_data data;
    kdml_mcmc_config config;
    kdml_mcmc_output output;
    kdml_cpu_backend_context cpu_context;
    SEXP result;
    SEXP value;
    int *state_kernel;
    double *state_theta;
    double *state_bandwidth;
    double *state_proposal_sd;
    int *adapt_count;
    kdml_rng *rng;
    int status;
    int error_chain = 0;
    int requested_backend;
    kdml_run_context run_context;
    SEXP unwind_continuation;

    kdml_setup_data(spec, &data, 1, 2);
    kdml_validate_mcmc_spec(spec, &data, &config);
    requested_backend = kdml_scalar_integer(
        kdml_list_get(spec, "backend"), "backend", KDML_BACKEND_CPU
    );
    if (requested_backend > KDML_BACKEND_AUTO) {
        Rf_error("Internal MCMC backend code is invalid.");
    }
    if (config.chain_threads > 1 && requested_backend != KDML_BACKEND_CPU) {
        Rf_error("Parallel chains require `backend = \"cpu\"`.");
    }
#ifndef _OPENMP
    if (config.chain_threads > 1) {
        Rf_error(
            "This KDML installation was built without OpenMP; "
            "use `chain_threads = 1` or reinstall with OpenMP support."
        );
    }
#endif
    if ((R_xlen_t) data.p > R_XLEN_T_MAX / config.draws ||
        (R_xlen_t) data.p * config.draws > R_XLEN_T_MAX / config.chains) {
        Rf_error("Requested retained-draw arrays are too large.");
    }
    result = PROTECT(Rf_allocVector(VECSXP, 16));
    value = kdml_allocate_array_in_list(
        result, 0, REALSXP, config.draws, data.p, config.chains
    );
    output.theta = REAL(value);
    value = kdml_allocate_array_in_list(
        result, 1, REALSXP, config.draws, data.p, config.chains
    );
    output.bandwidth = REAL(value);
    value = kdml_allocate_array_in_list(
        result, 2, INTSXP, config.draws, data.p, config.chains
    );
    output.kernel_code = INTEGER(value);
    value = kdml_allocate_matrix_in_list(
        result, 3, REALSXP, config.draws, config.chains
    );
    output.score = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 4, REALSXP, config.draws, config.chains
    );
    output.log_target = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 5, REALSXP, data.p, config.chains
    );
    output.final_proposal_sd = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 6, REALSXP, data.p, config.chains
    );
    output.kernel_attempts = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 7, REALSXP, data.p, config.chains
    );
    output.kernel_accepts = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 8, REALSXP, data.p, config.chains
    );
    output.bandwidth_attempts = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 9, REALSXP, data.p, config.chains
    );
    output.bandwidth_accepts = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 10, REALSXP, data.p, config.chains
    );
    output.sampling_kernel_attempts = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 11, REALSXP, data.p, config.chains
    );
    output.sampling_kernel_accepts = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 12, REALSXP, data.p, config.chains
    );
    output.sampling_bandwidth_attempts = REAL(value);
    value = kdml_allocate_matrix_in_list(
        result, 13, REALSXP, data.p, config.chains
    );
    output.sampling_bandwidth_accepts = REAL(value);
    value = PROTECT(Rf_allocVector(INTSXP, config.chains));
    SET_VECTOR_ELT(result, 14, value);
    UNPROTECT(1);
    output.invalid_proposals = INTEGER(VECTOR_ELT(result, 14));
    kdml_set_list_element(result, 15, R_NilValue);
    kdml_set_names(result, result_names, 16);

    state_kernel = (int *) R_alloc((size_t) data.p, sizeof(int));
    state_theta = (double *) R_alloc((size_t) data.p, sizeof(double));
    state_bandwidth = (double *) R_alloc((size_t) data.p, sizeof(double));
    state_proposal_sd = (double *) R_alloc((size_t) data.p, sizeof(double));
    adapt_count = (int *) R_alloc((size_t) data.p, sizeof(int));
    rng = (kdml_rng *) R_alloc((size_t) config.chains, sizeof(kdml_rng));
    kdml_seed_chains(rng, config.chains);

#ifdef _OPENMP
    if (config.chain_threads > 1) {
        int interrupted = 0;
        kdml_progress parallel_progress;
        kdml_progress_initialize(
            &parallel_progress, config.progress, config.progress_in_place,
            config.chains,
            (long long) (config.warmup + config.draws * config.thin)
        );
        status = kdml_run_parallel_cpu(
            &data, &config, &output, rng, &parallel_progress,
            &error_chain, &interrupted
        );
        if (interrupted) {
            UNPROTECT(1);
            Rf_error("MCMC sampling was interrupted.");
        }
        if (status != 0) {
            UNPROTECT(1);
            if (status == 6) {
                Rf_error(
                    "OpenMP could not start a CPU worker thread; "
                    "use `chain_threads = 1` or increase the OpenMP thread limit."
                );
            }
            if (status == 1) {
                Rf_error(
                    "Initial MCMC state for chain %d has a non-finite MSCV target; "
                    "adjust the initial bandwidths or kernels.",
                    error_chain + 1
                );
            }
            if (status == 2) {
                Rf_error(
                    "MCMC cache became numerically invalid in chain %d.",
                    error_chain + 1
                );
            }
            Rf_error(
                "MCMC retained an unexpected number of draws in chain %d.",
                error_chain + 1
            );
        }
        value = PROTECT(Rf_mkString("cpu"));
        SET_VECTOR_ELT(result, 15, value);
        UNPROTECT(1);
        UNPROTECT(1);
        return result;
    }
#endif

    memset(&cpu_context, 0, sizeof(cpu_context));

    memset(&run_context, 0, sizeof(run_context));
    run_context.data = &data;
    run_context.config = &config;
    run_context.output = &output;
    run_context.requested_backend = requested_backend;
    run_context.cpu_context = &cpu_context;
    run_context.state_kernel = state_kernel;
    run_context.state_theta = state_theta;
    run_context.state_bandwidth = state_bandwidth;
    run_context.state_proposal_sd = state_proposal_sd;
    run_context.adapt_count = adapt_count;
    run_context.rng = rng;
    run_context.error_chain = &error_chain;
    run_context.status = 0;
    kdml_progress_initialize(
        &run_context.progress, config.progress, config.progress_in_place,
        config.chains,
        config.warmup + config.draws * config.thin
    );

    unwind_continuation = PROTECT(R_MakeUnwindCont());
    R_UnwindProtect(
        kdml_run_mcmc_unwind_body, &run_context,
        kdml_cleanup_mcmc, &run_context, unwind_continuation
    );
    UNPROTECT(1);
    status = run_context.status;

    if (status != 0) {
        UNPROTECT(1);
        if (status == 4) {
            Rf_error(
                "%s",
                run_context.backend.error[0] == '\0'
                    ? "The scoring backend failed."
                    : run_context.backend.error
            );
        }
        if (status == 1) {
            Rf_error("Initial MCMC state for chain %d has a non-finite MSCV target; adjust the initial bandwidths or kernels.",
                     error_chain + 1);
        }
        if (status == 2) {
            Rf_error("MCMC cache became numerically invalid in chain %d.",
                     error_chain + 1);
        }
        Rf_error("MCMC retained an unexpected number of draws in chain %d.",
                 error_chain + 1);
    }

    value = PROTECT(Rf_mkString(run_context.backend.name));
    SET_VECTOR_ELT(result, 15, value);
    UNPROTECT(1);

    UNPROTECT(1);
    return result;
}
