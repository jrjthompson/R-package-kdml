#ifndef KDML_CUDA_API_H
#define KDML_CUDA_API_H

#include <stddef.h>

#if defined(_WIN32) && defined(__CUDACC__)
#define KDML_CUDA_API __declspec(dllexport)
#elif defined(KDML_CUDA_EXPORTS) && defined(__GNUC__)
#define KDML_CUDA_API __attribute__((visibility("default")))
#else
#define KDML_CUDA_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* These values deliberately match the internal KDML C representation. */
enum {
    KDML_CUDA_CONTINUOUS = 0,
    KDML_CUDA_NOMINAL = 1,
    KDML_CUDA_ORDINAL = 2
};

enum {
    KDML_CUDA_DKPS = 0,
    KDML_CUDA_DKSS = 1
};

typedef enum {
    KDML_CUDA_SUCCESS = 0,
    KDML_CUDA_UNAVAILABLE = 1,
    KDML_CUDA_INVALID_ARGUMENT = 2,
    KDML_CUDA_ALLOCATION_ERROR = 3,
    KDML_CUDA_RUNTIME_ERROR = 4,
    KDML_CUDA_INVALID_STATE = 5
} kdml_cuda_status;

typedef struct kdml_cuda_workspace kdml_cuda_workspace;

typedef struct {
    int device;
    char name[256];
    size_t total_global_memory;
    int compute_capability_major;
    int compute_capability_minor;
    int multiprocessor_count;
    int driver_version;
    int runtime_version;
} kdml_cuda_device_properties;

/* Runtime and build inspection.  The real CUDA translation unit returns 1
 * from kdml_cuda_compiled(); a non-CUDA stub may return 0 with the same ABI. */
KDML_CUDA_API int kdml_cuda_compiled(void);
KDML_CUDA_API int kdml_cuda_available(void);
KDML_CUDA_API int kdml_cuda_compile_version(void);
KDML_CUDA_API const char *kdml_cuda_compile_version_string(void);
KDML_CUDA_API int kdml_cuda_runtime_version(int *version);
KDML_CUDA_API int kdml_cuda_device_count(int *count);
KDML_CUDA_API int kdml_cuda_device_info(
    int device, kdml_cuda_device_properties *properties
);
KDML_CUDA_API const char *kdml_cuda_last_error(void);
KDML_CUDA_API const char *kdml_cuda_status_string(int status);

/* x is an n-by-p, column-major double matrix.  All input arrays are copied;
 * the caller need not retain them after this call.  A requested_device below
 * zero selects CUDA's current device. */
KDML_CUDA_API int kdml_cuda_create(
    int n, int p, int metric, const double *x, const int *type,
    const int *level_count, int continuous_count, int requested_device,
    kdml_cuda_workspace **workspace
);

/* Build establishes the current state and returns its mean leave-one-out
 * log-similarity.  A mathematically invalid state/proposal is represented by
 * KDML_CUDA_SUCCESS and score == -INFINITY; nonzero status means an API or
 * CUDA failure. */
KDML_CUDA_API int kdml_cuda_build(
    kdml_cuda_workspace *workspace, const int *kernel_code,
    const double *bandwidth, double *score
);
KDML_CUDA_API int kdml_cuda_propose(
    kdml_cuda_workspace *workspace, int feature, int proposed_kernel,
    double proposed_bandwidth, double *score
);
KDML_CUDA_API int kdml_cuda_accept(
    kdml_cuda_workspace *workspace, int feature, int proposed_kernel,
    double proposed_bandwidth
);
KDML_CUDA_API void kdml_cuda_destroy(kdml_cuda_workspace *workspace);

#ifdef __cplusplus
}
#endif

#endif
