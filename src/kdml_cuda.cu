#include "kdml_cuda_api.h"

#include <cuda_runtime.h>
#include <math_constants.h>

#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define KDML_CUDA_PI 3.141592653589793238462643383279502884
#define KDML_CUDA_SQRT_TWO 1.414213562373095048801688724209698079
#define KDML_CUDA_LOG_TWO 0.693147180559945309417232121458176568
#define KDML_CUDA_MAGIC UINT64_C(0x4b444d4c43554441)
#define KDML_CUDA_ERROR_LENGTH 1024
#define KDML_CUDA_BLOCK_X 16
#define KDML_CUDA_BLOCK_Y 16
#define KDML_CUDA_SCORE_THREADS 256

#define KDML_CUDA_STRINGIFY_INNER(value) #value
#define KDML_CUDA_STRINGIFY(value) KDML_CUDA_STRINGIFY_INNER(value)

#if defined(_MSC_VER)
static __declspec(thread) char kdml_cuda_error[KDML_CUDA_ERROR_LENGTH];
#elif defined(__GNUC__)
static __thread char kdml_cuda_error[KDML_CUDA_ERROR_LENGTH];
#else
static char kdml_cuda_error[KDML_CUDA_ERROR_LENGTH];
#endif

struct kdml_cuda_workspace {
    uint64_t magic;
    int n;
    int p;
    int metric;
    int continuous_count;
    int device;
    size_t pair_count;

    int built;
    int proposal_ready;
    int proposed_feature;
    int proposed_kernel;
    double proposed_bandwidth;

    int *host_type;
    int *host_level_count;
    int *host_kernel;
    double *host_bandwidth;

    double *device_x;
    int *device_type;
    int *device_level_count;
    int *device_kernel;
    double *device_bandwidth;
    double *device_current_rows;
    double *device_proposed_rows;
    double *device_current_product;
    double *device_proposed_product;
    int *device_current_zeros;
    int *device_proposed_zeros;
    int *device_invalid;
    double *device_score;
};

static void kdml_cuda_clear_error(void)
{
    kdml_cuda_error[0] = '\0';
}

static void kdml_cuda_set_error(const char *format, ...)
{
    va_list arguments;
    va_start(arguments, format);
#if defined(_MSC_VER)
    _vsnprintf_s(kdml_cuda_error, sizeof(kdml_cuda_error), _TRUNCATE,
                 format, arguments);
#else
    vsnprintf(kdml_cuda_error, sizeof(kdml_cuda_error), format, arguments);
#endif
    va_end(arguments);
    kdml_cuda_error[sizeof(kdml_cuda_error) - 1U] = '\0';
}

static int kdml_cuda_host_finite(double value)
{
#if defined(_MSC_VER)
    return _finite(value) != 0;
#else
    return isfinite(value) != 0;
#endif
}

static double kdml_cuda_host_infinity(void)
{
    union {
        uint64_t bits;
        double value;
    } infinity;
    infinity.bits = UINT64_C(0x7ff0000000000000);
    return infinity.value;
}

static int kdml_cuda_cuda_failure(cudaError_t error, const char *operation)
{
    int status = KDML_CUDA_RUNTIME_ERROR;
    if (error == cudaErrorMemoryAllocation) {
        status = KDML_CUDA_ALLOCATION_ERROR;
    } else if (error == cudaErrorNoDevice ||
               error == cudaErrorInsufficientDriver) {
        status = KDML_CUDA_UNAVAILABLE;
    }
    kdml_cuda_set_error("%s failed: %s", operation, cudaGetErrorString(error));
    return status;
}

static int kdml_cuda_valid_workspace(const kdml_cuda_workspace *workspace)
{
    return workspace != NULL && workspace->magic == KDML_CUDA_MAGIC;
}

static int kdml_cuda_select_workspace_device(kdml_cuda_workspace *workspace)
{
    cudaError_t error;
    if (!kdml_cuda_valid_workspace(workspace)) {
        kdml_cuda_set_error("CUDA workspace is null or invalid.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    error = cudaSetDevice(workspace->device);
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "cudaSetDevice");
    }
    return KDML_CUDA_SUCCESS;
}

static int kdml_cuda_kernel_valid(int type, int code)
{
    if (type == KDML_CUDA_CONTINUOUS) {
        return code >= 0 && code <= 10;
    }
    if (type == KDML_CUDA_NOMINAL) {
        return code >= 0 && code <= 1;
    }
    if (type == KDML_CUDA_ORDINAL) {
        return code >= 0 && code <= 4;
    }
    return 0;
}

static double kdml_cuda_bandwidth_upper(int type, int code, int level_count)
{
    if (type == KDML_CUDA_CONTINUOUS) {
        return kdml_cuda_host_infinity();
    }
    if (type == KDML_CUDA_NOMINAL && code == 1) {
        return ((double) level_count - 1.0) / (double) level_count;
    }
    return 1.0;
}

static int kdml_cuda_bandwidth_valid(int type, int code, int level_count,
                                     double bandwidth)
{
    const double upper = kdml_cuda_bandwidth_upper(type, code, level_count);
    return kdml_cuda_kernel_valid(type, code) &&
           kdml_cuda_host_finite(bandwidth) && bandwidth > 0.0 &&
           bandwidth < upper;
}

static double kdml_cuda_host_self_value(int type, int code, int level_count,
                                        int metric, double bandwidth)
{
    double value;
    if (!kdml_cuda_bandwidth_valid(type, code, level_count, bandwidth)) {
        return -kdml_cuda_host_infinity();
    }
    if (type == KDML_CUDA_CONTINUOUS) {
        switch (code) {
        case 0:
            value = 1.0 / sqrt(2.0 * KDML_CUDA_PI);
            break;
        case 1:
            value = 0.75;
            break;
        case 2:
            value = 0.5;
            break;
        case 3:
            value = 1.0;
            break;
        case 4:
            value = 15.0 / 16.0;
            break;
        case 5:
            value = 35.0 / 32.0;
            break;
        case 6:
            value = 70.0 / 81.0;
            break;
        case 7:
            value = KDML_CUDA_PI / 4.0;
            break;
        case 8:
            value = 0.25;
            break;
        case 9:
            value = 1.0 / KDML_CUDA_PI;
            break;
        case 10:
            value = 0.5 * sin(KDML_CUDA_PI / 4.0);
            break;
        default:
            return -kdml_cuda_host_infinity();
        }
        return metric == KDML_CUDA_DKPS ? value / bandwidth : value;
    }
    if (type == KDML_CUDA_NOMINAL) {
        return code == 0 ? 1.0 : 1.0 - bandwidth;
    }
    switch (code) {
    case 0:
        return 1.0 - bandwidth;
    case 1:
        return 1.0;
    case 2:
        return bandwidth;
    case 3:
        return pow(1.0 - bandwidth, (double) (level_count - 1));
    case 4:
        return 1.0;
    default:
        return -kdml_cuda_host_infinity();
    }
}

static int kdml_cuda_state_self_finite(const kdml_cuda_workspace *workspace,
                                       const int *kernel_code,
                                       const double *bandwidth,
                                       int replacement_feature,
                                       int replacement_kernel,
                                       double replacement_bandwidth)
{
    int feature;
    int has_continuous = 0;
    double continuous = workspace->metric == KDML_CUDA_DKPS ? 1.0 : 0.0;
    double categorical = 0.0;

    for (feature = 0; feature < workspace->p; ++feature) {
        const int code = feature == replacement_feature
                             ? replacement_kernel : kernel_code[feature];
        const double bw = feature == replacement_feature
                              ? replacement_bandwidth : bandwidth[feature];
        const int type = workspace->host_type[feature];
        const double value = kdml_cuda_host_self_value(
            type, code, workspace->host_level_count[feature],
            workspace->metric, bw
        );
        if (!kdml_cuda_host_finite(value)) {
            return 0;
        }
        if (type == KDML_CUDA_CONTINUOUS) {
            has_continuous = 1;
            if (workspace->metric == KDML_CUDA_DKPS) {
                continuous *= value;
                if (!kdml_cuda_host_finite(continuous)) {
                    return 0;
                }
            } else {
                continuous += value;
            }
        } else {
            categorical += value;
        }
    }
    if (workspace->metric == KDML_CUDA_DKPS && !has_continuous) {
        continuous = 0.0;
    }
    return kdml_cuda_host_finite(continuous + categorical);
}

__device__ static double kdml_cuda_continuous_kernel(int code, double z)
{
    const double absolute = fabs(z);
    double base;
    double tail;

    switch (code) {
    case 0:
        return exp(-0.5 * z * z) / sqrt(2.0 * KDML_CUDA_PI);
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
                   ? (KDML_CUDA_PI / 4.0) * cos(KDML_CUDA_PI * z / 2.0)
                   : 0.0;
    case 8:
        tail = exp(-absolute);
        return tail / ((1.0 + tail) * (1.0 + tail));
    case 9:
        tail = exp(-absolute);
        return (2.0 * tail) /
               (KDML_CUDA_PI * (1.0 + tail * tail));
    case 10:
        base = absolute / KDML_CUDA_SQRT_TWO;
        return 0.5 * exp(-base) * sin(base + KDML_CUDA_PI / 4.0);
    default:
        return nan("");
    }
}

__device__ static double kdml_cuda_feature_value(
    int metric, int type, int level_count, double first, double second,
    int kernel_code, double bandwidth)
{
    const double difference = first - second;
    const double distance = fabs(difference);
    double value;

    if (type == KDML_CUDA_CONTINUOUS) {
        value = kdml_cuda_continuous_kernel(
            kernel_code, difference / bandwidth
        );
        return metric == KDML_CUDA_DKPS ? value / bandwidth : value;
    }
    if (type == KDML_CUDA_NOMINAL) {
        if (kernel_code == 0) {
            return difference == 0.0 ? 1.0 : bandwidth;
        }
        if (difference == 0.0) {
            return 1.0 - bandwidth;
        }
        return bandwidth / ((double) level_count - 1.0);
    }

    {
        const int ordinal_distance = (int) floor(distance + 0.5);
        const int ordinal_range = level_count - 1;
        switch (kernel_code) {
        case 0:
            return ordinal_distance == 0
                       ? 1.0 - bandwidth
                       : 0.5 * (1.0 - bandwidth) *
                             pow(bandwidth, (double) ordinal_distance);
        case 1:
            return pow(bandwidth,
                       (double) ordinal_distance *
                       (double) ordinal_distance);
        case 2:
            return ordinal_distance == 0
                       ? bandwidth
                       : (1.0 - bandwidth) *
                             exp(-(double) ordinal_distance *
                                 KDML_CUDA_LOG_TWO);
        case 3:
            if (ordinal_distance > ordinal_range) {
                return 0.0;
            }
            value = lgamma((double) ordinal_range + 1.0) -
                    lgamma((double) ordinal_distance + 1.0) -
                    lgamma((double) (ordinal_range - ordinal_distance) +
                           1.0) +
                    (double) ordinal_distance * log(bandwidth) +
                    (double) (ordinal_range - ordinal_distance) *
                        log1p(-bandwidth);
            return exp(value);
        case 4:
            return pow(bandwidth, (double) ordinal_distance);
        default:
            return nan("");
        }
    }
}

__device__ static double kdml_cuda_atomic_add(double *address, double value)
{
#if __CUDA_ARCH__ >= 600
    return atomicAdd(address, value);
#else
    unsigned long long int *integer_address =
        (unsigned long long int *) address;
    unsigned long long int old = *integer_address;
    unsigned long long int assumed;
    do {
        assumed = old;
        old = atomicCAS(
            integer_address, assumed,
            __double_as_longlong(value + __longlong_as_double(assumed))
        );
    } while (assumed != old);
    return __longlong_as_double(old);
#endif
}

__device__ static size_t kdml_cuda_pair_index(int n, int first, int second)
{
    return ((size_t) first *
            (2U * (size_t) n - (size_t) first - 1U)) / 2U +
           (size_t) (second - first - 1);
}

__global__ static void kdml_cuda_build_kernel(
    int n, int p, int metric, int continuous_count,
    const double *x, const int *type, const int *level_count,
    const int *kernel_code, const double *bandwidth,
    double *row_sums, double *nonzero_product, int *zero_count,
    int *invalid)
{
    int first_start = (int) blockIdx.y * blockDim.y + threadIdx.y;
    int second_start = (int) blockIdx.x * blockDim.x + threadIdx.x;
    int first_stride = blockDim.y * gridDim.y;
    int second_stride = blockDim.x * gridDim.x;
    int first;

    for (first = first_start; first < n - 1; first += first_stride) {
        int second;
        for (second = second_start; second < n; second += second_stride) {
            int feature;
            int zeros = 0;
            double product = 1.0;
            double continuous_sum = 0.0;
            double categorical_sum = 0.0;
            double similarity;
            size_t pair;

            if (second <= first) {
                continue;
            }
            pair = kdml_cuda_pair_index(n, first, second);
            for (feature = 0; feature < p; ++feature) {
                const double value = kdml_cuda_feature_value(
                    metric, type[feature], level_count[feature],
                    x[(size_t) first + (size_t) n * (size_t) feature],
                    x[(size_t) second + (size_t) n * (size_t) feature],
                    kernel_code[feature], bandwidth[feature]
                );
                if (!isfinite(value)) {
                    atomicExch(invalid, 1);
                    return;
                }
                if (type[feature] == KDML_CUDA_CONTINUOUS) {
                    if (metric == KDML_CUDA_DKPS) {
                        if (value == 0.0) {
                            ++zeros;
                        } else {
                            product *= value;
                        }
                    } else {
                        continuous_sum += value;
                    }
                } else {
                    categorical_sum += value;
                }
            }
            if (metric == KDML_CUDA_DKPS) {
                if (continuous_count > 0) {
                    nonzero_product[pair] = product;
                    zero_count[pair] = zeros;
                }
                continuous_sum = continuous_count == 0 || zeros > 0
                                     ? 0.0 : product;
            }
            similarity = continuous_sum + categorical_sum;
            if (!isfinite(similarity)) {
                atomicExch(invalid, 1);
                return;
            }
            kdml_cuda_atomic_add(row_sums + first, similarity);
            kdml_cuda_atomic_add(row_sums + second, similarity);
        }
    }
}

__global__ static void kdml_cuda_probe_kernel(void)
{
}

__global__ static void kdml_cuda_proposal_kernel(
    int n, int p, int metric, int continuous_count, int feature,
    int proposed_kernel, double proposed_bandwidth,
    const double *x, const int *type, const int *level_count,
    const int *current_kernel, const double *current_bandwidth,
    const double *current_product, const int *current_zeros,
    double *proposed_product, int *proposed_zeros,
    double *proposed_rows, int *invalid)
{
    int first_start = (int) blockIdx.y * blockDim.y + threadIdx.y;
    int second_start = (int) blockIdx.x * blockDim.x + threadIdx.x;
    int first_stride = blockDim.y * gridDim.y;
    int second_stride = blockDim.x * gridDim.x;
    int first;

    for (first = first_start; first < n - 1; first += first_stride) {
        int second;
        for (second = second_start; second < n; second += second_stride) {
            double first_value;
            double second_value;
            double old_value;
            double new_value;
            double difference;
            size_t pair;

            if (second <= first) {
                continue;
            }
            first_value = x[
                (size_t) first + (size_t) n * (size_t) feature
            ];
            second_value = x[
                (size_t) second + (size_t) n * (size_t) feature
            ];
            old_value = kdml_cuda_feature_value(
                metric, type[feature], level_count[feature],
                first_value, second_value, current_kernel[feature],
                current_bandwidth[feature]
            );
            new_value = kdml_cuda_feature_value(
                metric, type[feature], level_count[feature],
                first_value, second_value, proposed_kernel,
                proposed_bandwidth
            );
            if (!isfinite(old_value) || !isfinite(new_value)) {
                atomicExch(invalid, 1);
                return;
            }
            pair = kdml_cuda_pair_index(n, first, second);
            if (metric == KDML_CUDA_DKPS &&
                type[feature] == KDML_CUDA_CONTINUOUS) {
                int other;
                int new_zero_count = 0;
                double new_product = 1.0;
                double old_continuous;
                double new_continuous;

                /* Use the same feature-order multiplication as a fresh
                 * build.  Repeated divide/multiply cache updates drift enough
                 * over long chains to change MSCV and MH decisions. */
                for (other = 0; other < p; ++other) {
                    double other_value;
                    if (type[other] != KDML_CUDA_CONTINUOUS) {
                        continue;
                    }
                    if (other == feature) {
                        other_value = new_value;
                    } else {
                        other_value = kdml_cuda_feature_value(
                            metric, type[other], level_count[other],
                            x[(size_t) first +
                              (size_t) n * (size_t) other],
                            x[(size_t) second +
                              (size_t) n * (size_t) other],
                            current_kernel[other],
                            current_bandwidth[other]
                        );
                    }
                    if (!isfinite(other_value)) {
                        atomicExch(invalid, 1);
                        return;
                    }
                    if (other_value == 0.0) {
                        ++new_zero_count;
                    } else {
                        new_product *= other_value;
                    }
                }
                old_continuous =
                    continuous_count == 0 || current_zeros[pair] > 0
                        ? 0.0 : current_product[pair];
                new_continuous =
                    continuous_count == 0 || new_zero_count > 0
                        ? 0.0 : new_product;
                proposed_product[pair] = new_product;
                proposed_zeros[pair] = new_zero_count;
                difference = new_continuous - old_continuous;
            } else {
                difference = new_value - old_value;
            }
            if (!isfinite(difference)) {
                atomicExch(invalid, 1);
                return;
            }
            kdml_cuda_atomic_add(proposed_rows + first, difference);
            kdml_cuda_atomic_add(proposed_rows + second, difference);
        }
    }
}

__global__ static void kdml_cuda_score_kernel(const double *row_sums, int n,
                                               const int *input_invalid,
                                               double *score)
{
    __shared__ double totals[KDML_CUDA_SCORE_THREADS];
    __shared__ int invalid[KDML_CUDA_SCORE_THREADS];
    const int thread = threadIdx.x;
    const double denominator = (double) n - 1.0;
    double local_total = 0.0;
    int local_invalid = *input_invalid != 0;
    int row;

    for (row = thread; row < n; row += blockDim.x) {
        const double leave_one_out = row_sums[row] / denominator;
        if (!isfinite(leave_one_out) || leave_one_out <= 0.0) {
            local_invalid = 1;
        } else {
            local_total += log(leave_one_out);
        }
    }
    totals[thread] = local_total;
    invalid[thread] = local_invalid;
    __syncthreads();

    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        if (thread < stride) {
            totals[thread] += totals[thread + stride];
            invalid[thread] |= invalid[thread + stride];
        }
        __syncthreads();
    }
    if (thread == 0) {
        const double answer = totals[0] / (double) n;
        *score = invalid[0] || !isfinite(answer)
                     ? -CUDART_INF : answer;
    }
}

static dim3 kdml_cuda_pair_grid(int n)
{
    unsigned int blocks =
        ((unsigned int) n + KDML_CUDA_BLOCK_X - 1U) /
        KDML_CUDA_BLOCK_X;
    if (blocks > 65535U) {
        blocks = 65535U;
    }
    return dim3(blocks, blocks, 1U);
}

static int kdml_cuda_launch_score(kdml_cuda_workspace *workspace,
                                  const double *device_rows,
                                  double *score)
{
    cudaError_t error;
    kdml_cuda_score_kernel<<<1, KDML_CUDA_SCORE_THREADS>>>(
        device_rows, workspace->n, workspace->device_invalid,
        workspace->device_score
    );
    error = cudaGetLastError();
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "CUDA score kernel launch");
    }
    error = cudaMemcpy(score, workspace->device_score, sizeof(double),
                       cudaMemcpyDeviceToHost);
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "CUDA score transfer");
    }
    return KDML_CUDA_SUCCESS;
}

static int kdml_cuda_allocate_device(void **pointer, size_t bytes,
                                     const char *description)
{
    cudaError_t error = cudaMalloc(pointer, bytes);
    if (error != cudaSuccess) {
        int status = kdml_cuda_cuda_failure(error, description);
        return status;
    }
    return KDML_CUDA_SUCCESS;
}

static int kdml_cuda_probe_selected_device(void)
{
    cudaError_t error;

    /* Discard a stale per-thread launch status before testing this module's
     * own executable image.  Synchronization also detects PTX/cubin and
     * device-runtime failures that are deferred past kernel launch. */
    (void) cudaGetLastError();
    kdml_cuda_probe_kernel<<<1, 1>>>();
    error = cudaGetLastError();
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(
            error, "CUDA execution-probe kernel launch"
        );
    }
    error = cudaDeviceSynchronize();
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(
            error, "CUDA execution-probe synchronization"
        );
    }
    return KDML_CUDA_SUCCESS;
}

extern "C" int kdml_cuda_compiled(void)
{
    return 1;
}

extern "C" int kdml_cuda_available(void)
{
    int count = 0;
    int device = 0;
    int status = kdml_cuda_device_count(&count);
    cudaError_t error;
    if (status != KDML_CUDA_SUCCESS || count < 1) {
        if (status == KDML_CUDA_SUCCESS) {
            kdml_cuda_set_error("No CUDA-capable device is available.");
        }
        return 0;
    }
    error = cudaGetDevice(&device);
    if (error != cudaSuccess) {
        (void) kdml_cuda_cuda_failure(error, "cudaGetDevice");
        return 0;
    }
    error = cudaSetDevice(device);
    if (error != cudaSuccess) {
        (void) kdml_cuda_cuda_failure(error, "cudaSetDevice");
        return 0;
    }
    status = kdml_cuda_probe_selected_device();
    if (status != KDML_CUDA_SUCCESS) {
        return 0;
    }
    kdml_cuda_clear_error();
    return 1;
}

extern "C" int kdml_cuda_compile_version(void)
{
    return CUDART_VERSION;
}

extern "C" const char *kdml_cuda_compile_version_string(void)
{
    return KDML_CUDA_STRINGIFY(CUDART_VERSION);
}

extern "C" int kdml_cuda_runtime_version(int *version)
{
    cudaError_t error;
    if (version == NULL) {
        kdml_cuda_set_error("CUDA runtime version output is null.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    *version = 0;
    error = cudaRuntimeGetVersion(version);
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "cudaRuntimeGetVersion");
    }
    kdml_cuda_clear_error();
    return KDML_CUDA_SUCCESS;
}

extern "C" int kdml_cuda_device_count(int *count)
{
    cudaError_t error;
    if (count == NULL) {
        kdml_cuda_set_error("CUDA device count output is null.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    *count = 0;
    error = cudaGetDeviceCount(count);
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "cudaGetDeviceCount");
    }
    kdml_cuda_clear_error();
    return KDML_CUDA_SUCCESS;
}

extern "C" int kdml_cuda_device_info(
    int device, kdml_cuda_device_properties *properties)
{
    cudaDeviceProp cuda_properties;
    cudaError_t error;
    int count;
    int status;

    if (properties == NULL) {
        kdml_cuda_set_error("CUDA device properties output is null.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    memset(properties, 0, sizeof(*properties));
    status = kdml_cuda_device_count(&count);
    if (status != KDML_CUDA_SUCCESS) {
        return status;
    }
    if (device < 0 || device >= count) {
        kdml_cuda_set_error("CUDA device index %d is outside [0, %d).",
                            device, count);
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    error = cudaGetDeviceProperties(&cuda_properties, device);
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "cudaGetDeviceProperties");
    }
    properties->device = device;
    strncpy(properties->name, cuda_properties.name,
            sizeof(properties->name) - 1U);
    properties->name[sizeof(properties->name) - 1U] = '\0';
    properties->total_global_memory = cuda_properties.totalGlobalMem;
    properties->compute_capability_major = cuda_properties.major;
    properties->compute_capability_minor = cuda_properties.minor;
    properties->multiprocessor_count = cuda_properties.multiProcessorCount;
    error = cudaDriverGetVersion(&properties->driver_version);
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "cudaDriverGetVersion");
    }
    error = cudaRuntimeGetVersion(&properties->runtime_version);
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "cudaRuntimeGetVersion");
    }
    kdml_cuda_clear_error();
    return KDML_CUDA_SUCCESS;
}

extern "C" const char *kdml_cuda_last_error(void)
{
    return kdml_cuda_error;
}

extern "C" const char *kdml_cuda_status_string(int status)
{
    switch (status) {
    case KDML_CUDA_SUCCESS:
        return "success";
    case KDML_CUDA_UNAVAILABLE:
        return "CUDA unavailable";
    case KDML_CUDA_INVALID_ARGUMENT:
        return "invalid argument";
    case KDML_CUDA_ALLOCATION_ERROR:
        return "CUDA allocation error";
    case KDML_CUDA_RUNTIME_ERROR:
        return "CUDA runtime error";
    case KDML_CUDA_INVALID_STATE:
        return "invalid CUDA workspace state";
    default:
        return "unknown CUDA status";
    }
}

extern "C" int kdml_cuda_create(
    int n, int p, int metric, const double *x, const int *type,
    const int *level_count, int continuous_count, int requested_device,
    kdml_cuda_workspace **output)
{
    kdml_cuda_workspace *workspace;
    size_t cells;
    size_t pair_count;
    int observed_continuous = 0;
    int feature;
    size_t cell;
    int count;
    int device;
    int status;
    cudaError_t error;

    kdml_cuda_clear_error();
    if (output == NULL) {
        kdml_cuda_set_error("CUDA workspace output is null.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    *output = NULL;
    if (n < 2 || p < 1 || x == NULL || type == NULL ||
        level_count == NULL ||
        (metric != KDML_CUDA_DKPS && metric != KDML_CUDA_DKSS)) {
        kdml_cuda_set_error("CUDA workspace dimensions, metric, or inputs are invalid.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    if ((size_t) n > SIZE_MAX / (size_t) p) {
        kdml_cuda_set_error("CUDA data matrix dimensions overflow size_t.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    cells = (size_t) n * (size_t) p;
    if ((size_t) n > SIZE_MAX / (size_t) (n - 1)) {
        kdml_cuda_set_error("CUDA pair count overflows size_t.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    pair_count = (size_t) n * (size_t) (n - 1) / 2U;
    if (cells > SIZE_MAX / sizeof(double) ||
        (size_t) p > SIZE_MAX / sizeof(double) ||
        (size_t) p > SIZE_MAX / sizeof(int) ||
        (size_t) n > SIZE_MAX / sizeof(double) ||
        pair_count > SIZE_MAX / sizeof(double) ||
        pair_count > SIZE_MAX / sizeof(int)) {
        kdml_cuda_set_error("CUDA allocation sizes overflow size_t.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    for (feature = 0; feature < p; ++feature) {
        if (type[feature] < KDML_CUDA_CONTINUOUS ||
            type[feature] > KDML_CUDA_ORDINAL) {
            kdml_cuda_set_error("Feature %d has invalid type code %d.",
                                feature, type[feature]);
            return KDML_CUDA_INVALID_ARGUMENT;
        }
        if (type[feature] == KDML_CUDA_CONTINUOUS) {
            ++observed_continuous;
        } else if (level_count[feature] < 2) {
            kdml_cuda_set_error(
                "Categorical feature %d has fewer than two levels.",
                feature
            );
            return KDML_CUDA_INVALID_ARGUMENT;
        }
    }
    if (continuous_count != observed_continuous) {
        kdml_cuda_set_error(
            "Continuous feature count %d does not match metadata count %d.",
            continuous_count, observed_continuous
        );
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    for (cell = 0; cell < cells; ++cell) {
        if (!kdml_cuda_host_finite(x[cell])) {
            kdml_cuda_set_error("CUDA data matrix contains a non-finite value.");
            return KDML_CUDA_INVALID_ARGUMENT;
        }
    }
    for (feature = 0; feature < p; ++feature) {
        int row;
        if (type[feature] == KDML_CUDA_CONTINUOUS) {
            continue;
        }
        for (row = 0; row < n; ++row) {
            const double value = x[
                (size_t) row + (size_t) n * (size_t) feature
            ];
            if (value < 0.0 ||
                value > (double) level_count[feature] - 1.0 ||
                value != floor(value)) {
                kdml_cuda_set_error(
                    "Categorical feature %d has invalid encoded data.",
                    feature
                );
                return KDML_CUDA_INVALID_ARGUMENT;
            }
        }
    }

    status = kdml_cuda_device_count(&count);
    if (status != KDML_CUDA_SUCCESS) {
        return status;
    }
    if (count < 1) {
        kdml_cuda_set_error("No CUDA-capable device is available.");
        return KDML_CUDA_UNAVAILABLE;
    }
    if (requested_device < 0) {
        error = cudaGetDevice(&device);
        if (error != cudaSuccess) {
            return kdml_cuda_cuda_failure(error, "cudaGetDevice");
        }
    } else {
        device = requested_device;
    }
    if (device < 0 || device >= count) {
        kdml_cuda_set_error("CUDA device index %d is outside [0, %d).",
                            device, count);
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    error = cudaSetDevice(device);
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "cudaSetDevice");
    }
    status = kdml_cuda_probe_selected_device();
    if (status != KDML_CUDA_SUCCESS) {
        return status;
    }

    workspace = (kdml_cuda_workspace *) calloc(1, sizeof(*workspace));
    if (workspace == NULL) {
        kdml_cuda_set_error("Could not allocate the CUDA host workspace.");
        return KDML_CUDA_ALLOCATION_ERROR;
    }
    workspace->magic = KDML_CUDA_MAGIC;
    workspace->n = n;
    workspace->p = p;
    workspace->metric = metric;
    workspace->continuous_count = continuous_count;
    workspace->device = device;
    workspace->pair_count = pair_count;

    workspace->host_type = (int *) malloc((size_t) p * sizeof(int));
    workspace->host_level_count =
        (int *) malloc((size_t) p * sizeof(int));
    workspace->host_kernel = (int *) malloc((size_t) p * sizeof(int));
    workspace->host_bandwidth =
        (double *) malloc((size_t) p * sizeof(double));
    if (workspace->host_type == NULL ||
        workspace->host_level_count == NULL ||
        workspace->host_kernel == NULL ||
        workspace->host_bandwidth == NULL) {
        kdml_cuda_set_error("Could not allocate CUDA host state arrays.");
        kdml_cuda_destroy(workspace);
        return KDML_CUDA_ALLOCATION_ERROR;
    }
    memcpy(workspace->host_type, type, (size_t) p * sizeof(int));
    memcpy(workspace->host_level_count, level_count,
           (size_t) p * sizeof(int));

#define KDML_CUDA_ALLOCATE(member, bytes, description)                 \
    do {                                                               \
        status = kdml_cuda_allocate_device(                            \
            (void **) &(workspace->member), (bytes), (description)     \
        );                                                             \
        if (status != KDML_CUDA_SUCCESS) {                             \
            kdml_cuda_destroy(workspace);                              \
            return status;                                             \
        }                                                              \
    } while (0)

    KDML_CUDA_ALLOCATE(device_x, cells * sizeof(double),
                       "CUDA data allocation");
    KDML_CUDA_ALLOCATE(device_type, (size_t) p * sizeof(int),
                       "CUDA type allocation");
    KDML_CUDA_ALLOCATE(device_level_count, (size_t) p * sizeof(int),
                       "CUDA level-count allocation");
    KDML_CUDA_ALLOCATE(device_kernel, (size_t) p * sizeof(int),
                       "CUDA kernel-state allocation");
    KDML_CUDA_ALLOCATE(device_bandwidth, (size_t) p * sizeof(double),
                       "CUDA bandwidth-state allocation");
    KDML_CUDA_ALLOCATE(device_current_rows, (size_t) n * sizeof(double),
                       "CUDA current row-cache allocation");
    KDML_CUDA_ALLOCATE(device_proposed_rows, (size_t) n * sizeof(double),
                       "CUDA proposal row-cache allocation");
    if (metric == KDML_CUDA_DKPS && continuous_count > 0) {
        KDML_CUDA_ALLOCATE(device_current_product,
                           pair_count * sizeof(double),
                           "CUDA current product-cache allocation");
        KDML_CUDA_ALLOCATE(device_proposed_product,
                           pair_count * sizeof(double),
                           "CUDA proposal product-cache allocation");
        KDML_CUDA_ALLOCATE(device_current_zeros,
                           pair_count * sizeof(int),
                           "CUDA current zero-cache allocation");
        KDML_CUDA_ALLOCATE(device_proposed_zeros,
                           pair_count * sizeof(int),
                           "CUDA proposal zero-cache allocation");
    }
    KDML_CUDA_ALLOCATE(device_invalid, sizeof(int),
                       "CUDA validity-flag allocation");
    KDML_CUDA_ALLOCATE(device_score, sizeof(double),
                       "CUDA score allocation");
#undef KDML_CUDA_ALLOCATE

    error = cudaMemcpy(workspace->device_x, x, cells * sizeof(double),
                       cudaMemcpyHostToDevice);
    if (error == cudaSuccess) {
        error = cudaMemcpy(workspace->device_type, type,
                           (size_t) p * sizeof(int),
                           cudaMemcpyHostToDevice);
    }
    if (error == cudaSuccess) {
        error = cudaMemcpy(workspace->device_level_count, level_count,
                           (size_t) p * sizeof(int),
                           cudaMemcpyHostToDevice);
    }
    if (error != cudaSuccess) {
        status = kdml_cuda_cuda_failure(error,
                                        "CUDA immutable-data transfer");
        kdml_cuda_destroy(workspace);
        return status;
    }
    *output = workspace;
    kdml_cuda_clear_error();
    return KDML_CUDA_SUCCESS;
}

extern "C" int kdml_cuda_build(
    kdml_cuda_workspace *workspace, const int *kernel_code,
    const double *bandwidth, double *score)
{
    cudaError_t error;
    int status;
    dim3 block(KDML_CUDA_BLOCK_X, KDML_CUDA_BLOCK_Y, 1U);
    dim3 grid;

    kdml_cuda_clear_error();
    if (!kdml_cuda_valid_workspace(workspace) || kernel_code == NULL ||
        bandwidth == NULL || score == NULL) {
        kdml_cuda_set_error("CUDA build arguments are null or invalid.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    *score = -kdml_cuda_host_infinity();
    workspace->built = 0;
    workspace->proposal_ready = 0;
    if (!kdml_cuda_state_self_finite(
            workspace, kernel_code, bandwidth, -1, 0, 0.0)) {
        return KDML_CUDA_SUCCESS;
    }
    status = kdml_cuda_select_workspace_device(workspace);
    if (status != KDML_CUDA_SUCCESS) {
        return status;
    }
    error = cudaMemcpy(workspace->device_kernel, kernel_code,
                       (size_t) workspace->p * sizeof(int),
                       cudaMemcpyHostToDevice);
    if (error == cudaSuccess) {
        error = cudaMemcpy(workspace->device_bandwidth, bandwidth,
                           (size_t) workspace->p * sizeof(double),
                           cudaMemcpyHostToDevice);
    }
    if (error == cudaSuccess) {
        error = cudaMemset(workspace->device_current_rows, 0,
                           (size_t) workspace->n * sizeof(double));
    }
    if (error == cudaSuccess) {
        error = cudaMemset(workspace->device_invalid, 0, sizeof(int));
    }
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error,
                                      "CUDA build-state preparation");
    }

    grid = kdml_cuda_pair_grid(workspace->n);
    kdml_cuda_build_kernel<<<grid, block>>>(
        workspace->n, workspace->p, workspace->metric,
        workspace->continuous_count, workspace->device_x,
        workspace->device_type, workspace->device_level_count,
        workspace->device_kernel, workspace->device_bandwidth,
        workspace->device_current_rows,
        workspace->device_current_product,
        workspace->device_current_zeros,
        workspace->device_invalid
    );
    error = cudaGetLastError();
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "CUDA build kernel launch");
    }
    status = kdml_cuda_launch_score(
        workspace, workspace->device_current_rows, score
    );
    if (status != KDML_CUDA_SUCCESS) {
        return status;
    }
    if (kdml_cuda_host_finite(*score)) {
        memcpy(workspace->host_kernel, kernel_code,
               (size_t) workspace->p * sizeof(int));
        memcpy(workspace->host_bandwidth, bandwidth,
               (size_t) workspace->p * sizeof(double));
        workspace->built = 1;
    }
    kdml_cuda_clear_error();
    return KDML_CUDA_SUCCESS;
}

extern "C" int kdml_cuda_propose(
    kdml_cuda_workspace *workspace, int feature, int proposed_kernel,
    double proposed_bandwidth, double *score)
{
    cudaError_t error;
    int status;
    dim3 block(KDML_CUDA_BLOCK_X, KDML_CUDA_BLOCK_Y, 1U);
    dim3 grid;

    kdml_cuda_clear_error();
    if (!kdml_cuda_valid_workspace(workspace) || score == NULL) {
        kdml_cuda_set_error("CUDA proposal arguments are null or invalid.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    *score = -kdml_cuda_host_infinity();
    workspace->proposal_ready = 0;
    if (!workspace->built) {
        kdml_cuda_set_error("CUDA proposal requires a valid built state.");
        return KDML_CUDA_INVALID_STATE;
    }
    if (feature < 0 || feature >= workspace->p) {
        kdml_cuda_set_error("CUDA proposal feature %d is out of range.",
                            feature);
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    if (!kdml_cuda_state_self_finite(
            workspace, workspace->host_kernel,
            workspace->host_bandwidth, feature, proposed_kernel,
            proposed_bandwidth)) {
        return KDML_CUDA_SUCCESS;
    }
    status = kdml_cuda_select_workspace_device(workspace);
    if (status != KDML_CUDA_SUCCESS) {
        return status;
    }
    error = cudaMemcpy(workspace->device_proposed_rows,
                       workspace->device_current_rows,
                       (size_t) workspace->n * sizeof(double),
                       cudaMemcpyDeviceToDevice);
    if (error == cudaSuccess) {
        error = cudaMemset(workspace->device_invalid, 0, sizeof(int));
    }
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error,
                                      "CUDA proposal-state preparation");
    }

    grid = kdml_cuda_pair_grid(workspace->n);
    kdml_cuda_proposal_kernel<<<grid, block>>>(
        workspace->n, workspace->p, workspace->metric,
        workspace->continuous_count, feature, proposed_kernel,
        proposed_bandwidth, workspace->device_x,
        workspace->device_type, workspace->device_level_count,
        workspace->device_kernel, workspace->device_bandwidth,
        workspace->device_current_product,
        workspace->device_current_zeros,
        workspace->device_proposed_product,
        workspace->device_proposed_zeros,
        workspace->device_proposed_rows,
        workspace->device_invalid
    );
    error = cudaGetLastError();
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error,
                                      "CUDA proposal kernel launch");
    }
    status = kdml_cuda_launch_score(
        workspace, workspace->device_proposed_rows, score
    );
    if (status != KDML_CUDA_SUCCESS) {
        return status;
    }
    if (kdml_cuda_host_finite(*score)) {
        workspace->proposal_ready = 1;
        workspace->proposed_feature = feature;
        workspace->proposed_kernel = proposed_kernel;
        workspace->proposed_bandwidth = proposed_bandwidth;
    }
    kdml_cuda_clear_error();
    return KDML_CUDA_SUCCESS;
}

extern "C" int kdml_cuda_accept(
    kdml_cuda_workspace *workspace, int feature, int proposed_kernel,
    double proposed_bandwidth)
{
    cudaError_t error;
    int status;
    double *row_swap;
    double *product_swap;
    int *zero_swap;

    kdml_cuda_clear_error();
    if (!kdml_cuda_valid_workspace(workspace)) {
        kdml_cuda_set_error("CUDA acceptance workspace is null or invalid.");
        return KDML_CUDA_INVALID_ARGUMENT;
    }
    if (!workspace->built || !workspace->proposal_ready) {
        kdml_cuda_set_error("CUDA acceptance has no valid pending proposal.");
        return KDML_CUDA_INVALID_STATE;
    }
    if (feature != workspace->proposed_feature ||
        proposed_kernel != workspace->proposed_kernel ||
        proposed_bandwidth != workspace->proposed_bandwidth) {
        kdml_cuda_set_error(
            "CUDA acceptance does not match the pending proposal."
        );
        return KDML_CUDA_INVALID_STATE;
    }
    status = kdml_cuda_select_workspace_device(workspace);
    if (status != KDML_CUDA_SUCCESS) {
        return status;
    }
    error = cudaMemcpy(
        workspace->device_kernel + feature, &proposed_kernel,
        sizeof(int), cudaMemcpyHostToDevice
    );
    if (error == cudaSuccess) {
        error = cudaMemcpy(
            workspace->device_bandwidth + feature, &proposed_bandwidth,
            sizeof(double), cudaMemcpyHostToDevice
        );
    }
    if (error != cudaSuccess) {
        return kdml_cuda_cuda_failure(error, "CUDA proposal commit");
    }
    row_swap = workspace->device_current_rows;
    workspace->device_current_rows = workspace->device_proposed_rows;
    workspace->device_proposed_rows = row_swap;
    if (workspace->metric == KDML_CUDA_DKPS &&
        workspace->host_type[feature] == KDML_CUDA_CONTINUOUS) {
        product_swap = workspace->device_current_product;
        workspace->device_current_product =
            workspace->device_proposed_product;
        workspace->device_proposed_product = product_swap;
        zero_swap = workspace->device_current_zeros;
        workspace->device_current_zeros = workspace->device_proposed_zeros;
        workspace->device_proposed_zeros = zero_swap;
    }
    workspace->host_kernel[feature] = proposed_kernel;
    workspace->host_bandwidth[feature] = proposed_bandwidth;
    workspace->proposal_ready = 0;
    kdml_cuda_clear_error();
    return KDML_CUDA_SUCCESS;
}

extern "C" void kdml_cuda_destroy(kdml_cuda_workspace *workspace)
{
    if (!kdml_cuda_valid_workspace(workspace)) {
        return;
    }
    (void) cudaSetDevice(workspace->device);
    (void) cudaFree(workspace->device_x);
    (void) cudaFree(workspace->device_type);
    (void) cudaFree(workspace->device_level_count);
    (void) cudaFree(workspace->device_kernel);
    (void) cudaFree(workspace->device_bandwidth);
    (void) cudaFree(workspace->device_current_rows);
    (void) cudaFree(workspace->device_proposed_rows);
    (void) cudaFree(workspace->device_current_product);
    (void) cudaFree(workspace->device_proposed_product);
    (void) cudaFree(workspace->device_current_zeros);
    (void) cudaFree(workspace->device_proposed_zeros);
    (void) cudaFree(workspace->device_invalid);
    (void) cudaFree(workspace->device_score);
    free(workspace->host_type);
    free(workspace->host_level_count);
    free(workspace->host_kernel);
    free(workspace->host_bandwidth);
    workspace->magic = 0;
    free(workspace);
}
