#if defined(KDML_CUDA_SIDECAR) && !defined(_WIN32) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include "kdml_cuda_api.h"

#include <stdio.h>
#include <string.h>

#if defined(KDML_CUDA_SIDECAR)

#if defined(_WIN32)
#include <stdint.h>
#include <wchar.h>
#include <windows.h>
#else
#include <dlfcn.h>
#include <stdlib.h>
#endif

typedef struct {
    int (*compiled)(void);
    int (*available)(void);
    int (*compile_version)(void);
    const char *(*compile_version_string)(void);
    int (*runtime_version)(int *);
    int (*device_count)(int *);
    int (*device_info)(int, kdml_cuda_device_properties *);
    const char *(*last_error)(void);
    const char *(*status_string)(int);
    int (*create)(int, int, int, const double *, const int *, const int *,
                  int, int, kdml_cuda_workspace **);
    int (*build)(kdml_cuda_workspace *, const int *, const double *, double *);
    int (*propose)(kdml_cuda_workspace *, int, int, double, double *);
    int (*accept)(kdml_cuda_workspace *, int, int, double);
    void (*destroy)(kdml_cuda_workspace *);
} kdml_cuda_function_table;

static kdml_cuda_function_table kdml_cuda_api;
static int kdml_cuda_load_attempted = 0;
static int kdml_cuda_loaded = 0;
static char kdml_cuda_loader_error[512] = "CUDA sidecar has not been loaded";

static void kdml_cuda_set_loader_error(const char *message)
{
    if (message == NULL || message[0] == '\0') {
        message = "unknown dynamic-loader error";
    }
    (void) snprintf(kdml_cuda_loader_error,
                    sizeof(kdml_cuda_loader_error), "%s", message);
}

#if defined(_WIN32)

typedef HMODULE kdml_cuda_library;

static kdml_cuda_library kdml_cuda_open_library(void)
{
    HMODULE package_module = NULL;
    wchar_t package_path[32768];
    wchar_t *separator;
    const wchar_t sidecar_name[] = L"kdmlcuda.dll";
    DWORD length;

    if (!GetModuleHandleExW(
            GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            (LPCWSTR) (uintptr_t) &kdml_cuda_open_library,
            &package_module)) {
        kdml_cuda_set_loader_error(
            "could not locate the loaded kdml package DLL");
        return NULL;
    }

    length = GetModuleFileNameW(package_module, package_path,
                                (DWORD) (sizeof(package_path) /
                                         sizeof(package_path[0])));
    if (length == 0 ||
        length >= (DWORD) (sizeof(package_path) / sizeof(package_path[0]))) {
        kdml_cuda_set_loader_error(
            "could not determine the kdml package DLL path");
        return NULL;
    }

    separator = wcsrchr(package_path, L'\\');
    if (separator == NULL) {
        separator = wcsrchr(package_path, L'/');
    }
    if (separator == NULL) {
        kdml_cuda_set_loader_error(
            "could not determine the kdml package DLL directory");
        return NULL;
    }
    ++separator;
    if ((size_t) (separator - package_path) +
            (sizeof(sidecar_name) / sizeof(sidecar_name[0])) >
        (sizeof(package_path) / sizeof(package_path[0]))) {
        kdml_cuda_set_loader_error("the CUDA sidecar path is too long");
        return NULL;
    }
    (void) wcscpy(separator, sidecar_name);

    package_module = LoadLibraryW(package_path);
    if (package_module == NULL) {
        char message[256];
        (void) snprintf(message, sizeof(message),
                        "could not load kdmlcuda.dll (Windows error %lu)",
                        (unsigned long) GetLastError());
        kdml_cuda_set_loader_error(message);
    }
    return package_module;
}

static void *kdml_cuda_find_symbol(kdml_cuda_library library,
                                   const char *name)
{
    FARPROC address = GetProcAddress(library, name);
    void *answer = NULL;
    if (address != NULL) {
        size_t count = sizeof(address) < sizeof(answer) ?
            sizeof(address) : sizeof(answer);
        memcpy(&answer, &address, count);
    }
    return answer;
}

#else

typedef void *kdml_cuda_library;

static kdml_cuda_library kdml_cuda_open_library(void)
{
    Dl_info package_info;
    char *sidecar_path;
    char *separator;
    void *library;
    const char sidecar_name[] = "kdmlcuda.so";
    size_t required;

    if (dladdr((void *) &kdml_cuda_open_library, &package_info) == 0 ||
        package_info.dli_fname == NULL) {
        kdml_cuda_set_loader_error(
            "could not locate the loaded kdml package shared library");
        return NULL;
    }

    required = strlen(package_info.dli_fname) + sizeof(sidecar_name) + 1;
    sidecar_path = (char *) malloc(required);
    if (sidecar_path == NULL) {
        kdml_cuda_set_loader_error(
            "could not allocate the CUDA sidecar path");
        return NULL;
    }
    (void) snprintf(sidecar_path, required, "%s", package_info.dli_fname);
    separator = strrchr(sidecar_path, '/');
    if (separator == NULL) {
        free(sidecar_path);
        kdml_cuda_set_loader_error(
            "could not determine the kdml shared-library directory");
        return NULL;
    }
    ++separator;
    (void) snprintf(separator,
                    required - (size_t) (separator - sidecar_path),
                    "%s", sidecar_name);

    library = dlopen(sidecar_path, RTLD_NOW | RTLD_LOCAL);
    free(sidecar_path);
    if (library == NULL) {
        kdml_cuda_set_loader_error(dlerror());
    }
    return library;
}

static void *kdml_cuda_find_symbol(kdml_cuda_library library,
                                   const char *name)
{
    (void) dlerror();
    return dlsym(library, name);
}

#endif

static int kdml_cuda_resolve(void *destination, size_t destination_size,
                             kdml_cuda_library library, const char *name)
{
    void *address = kdml_cuda_find_symbol(library, name);
    const size_t count = destination_size < sizeof(address) ?
        destination_size : sizeof(address);
    if (address == NULL) {
        char message[512];
        (void) snprintf(message, sizeof(message),
                        "CUDA sidecar is missing required symbol `%s`", name);
        kdml_cuda_set_loader_error(message);
        return 0;
    }
    memset(destination, 0, destination_size);
    memcpy(destination, &address, count);
    return 1;
}

#define KDML_CUDA_RESOLVE(field, name)                                      \
    do {                                                                     \
        if (!kdml_cuda_resolve(&kdml_cuda_api.field,                         \
                               sizeof(kdml_cuda_api.field), library, name)) { \
            memset(&kdml_cuda_api, 0, sizeof(kdml_cuda_api));                \
            return 0;                                                        \
        }                                                                    \
    } while (0)

static int kdml_cuda_load(void)
{
    kdml_cuda_library library;
    if (kdml_cuda_load_attempted) {
        return kdml_cuda_loaded;
    }
    kdml_cuda_load_attempted = 1;
    library = kdml_cuda_open_library();
    if (library == NULL) {
        return 0;
    }

    KDML_CUDA_RESOLVE(compiled, "kdmlcuda_compiled");
    KDML_CUDA_RESOLVE(available, "kdmlcuda_available");
    KDML_CUDA_RESOLVE(compile_version, "kdmlcuda_compile_version");
    KDML_CUDA_RESOLVE(compile_version_string,
                      "kdmlcuda_compile_version_string");
    KDML_CUDA_RESOLVE(runtime_version, "kdmlcuda_runtime_version");
    KDML_CUDA_RESOLVE(device_count, "kdmlcuda_device_count");
    KDML_CUDA_RESOLVE(device_info, "kdmlcuda_device_info");
    KDML_CUDA_RESOLVE(last_error, "kdmlcuda_last_error");
    KDML_CUDA_RESOLVE(status_string, "kdmlcuda_status_string");
    KDML_CUDA_RESOLVE(create, "kdmlcuda_create");
    KDML_CUDA_RESOLVE(build, "kdmlcuda_build");
    KDML_CUDA_RESOLVE(propose, "kdmlcuda_propose");
    KDML_CUDA_RESOLVE(accept, "kdmlcuda_accept");
    KDML_CUDA_RESOLVE(destroy, "kdmlcuda_destroy");

    kdml_cuda_loaded = 1;
    return 1;
}

#undef KDML_CUDA_RESOLVE

int kdml_cuda_compiled(void)
{
    return 1;
}

int kdml_cuda_available(void)
{
    return kdml_cuda_load() ? kdml_cuda_api.available() : 0;
}

int kdml_cuda_compile_version(void)
{
    return kdml_cuda_load() ? kdml_cuda_api.compile_version() : 0;
}

const char *kdml_cuda_compile_version_string(void)
{
    return kdml_cuda_load() ? kdml_cuda_api.compile_version_string() :
        "unavailable";
}

int kdml_cuda_runtime_version(int *version)
{
    if (version != NULL) {
        *version = 0;
    }
    return kdml_cuda_load() ? kdml_cuda_api.runtime_version(version) :
        KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_device_count(int *count)
{
    if (count != NULL) {
        *count = 0;
    }
    return kdml_cuda_load() ? kdml_cuda_api.device_count(count) :
        KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_device_info(int device,
                          kdml_cuda_device_properties *properties)
{
    if (properties != NULL) {
        memset(properties, 0, sizeof(*properties));
    }
    return kdml_cuda_load() ?
        kdml_cuda_api.device_info(device, properties) :
        KDML_CUDA_UNAVAILABLE;
}

const char *kdml_cuda_last_error(void)
{
    if (kdml_cuda_load() && kdml_cuda_api.last_error != NULL) {
        return kdml_cuda_api.last_error();
    }
    return kdml_cuda_loader_error;
}

const char *kdml_cuda_status_string(int status)
{
    if (kdml_cuda_load() && kdml_cuda_api.status_string != NULL) {
        return kdml_cuda_api.status_string(status);
    }
    switch (status) {
    case KDML_CUDA_SUCCESS:
        return "success";
    case KDML_CUDA_UNAVAILABLE:
        return "CUDA unavailable";
    case KDML_CUDA_INVALID_ARGUMENT:
        return "invalid CUDA argument";
    case KDML_CUDA_ALLOCATION_ERROR:
        return "CUDA allocation error";
    case KDML_CUDA_RUNTIME_ERROR:
        return "CUDA runtime error";
    case KDML_CUDA_INVALID_STATE:
        return "invalid CUDA state";
    default:
        return "unknown CUDA status";
    }
}

int kdml_cuda_create(int n, int p, int metric,
                     const double *x, const int *type,
                     const int *level_count, int continuous_count,
                     int requested_device,
                     kdml_cuda_workspace **workspace)
{
    if (workspace != NULL) {
        *workspace = NULL;
    }
    return kdml_cuda_load() ?
        kdml_cuda_api.create(n, p, metric, x, type, level_count,
                             continuous_count, requested_device, workspace) :
        KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_build(kdml_cuda_workspace *workspace,
                    const int *kernel_code, const double *bandwidth,
                    double *score)
{
    return kdml_cuda_load() ?
        kdml_cuda_api.build(workspace, kernel_code, bandwidth, score) :
        KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_propose(kdml_cuda_workspace *workspace,
                      int feature, int proposed_kernel,
                      double proposed_bandwidth, double *score)
{
    return kdml_cuda_load() ?
        kdml_cuda_api.propose(workspace, feature, proposed_kernel,
                              proposed_bandwidth, score) :
        KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_accept(kdml_cuda_workspace *workspace,
                     int feature, int proposed_kernel,
                     double proposed_bandwidth)
{
    return kdml_cuda_load() ?
        kdml_cuda_api.accept(workspace, feature, proposed_kernel,
                             proposed_bandwidth) :
        KDML_CUDA_UNAVAILABLE;
}

void kdml_cuda_destroy(kdml_cuda_workspace *workspace)
{
    if (workspace != NULL && kdml_cuda_load()) {
        kdml_cuda_api.destroy(workspace);
    }
}

#else

int kdml_cuda_compiled(void)
{
    return 0;
}

int kdml_cuda_available(void)
{
    return 0;
}

int kdml_cuda_compile_version(void)
{
    return 0;
}

const char *kdml_cuda_compile_version_string(void)
{
    return "not compiled";
}

int kdml_cuda_runtime_version(int *version)
{
    if (version != NULL) {
        *version = 0;
    }
    return KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_device_count(int *count)
{
    if (count != NULL) {
        *count = 0;
    }
    return KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_device_info(int device,
                          kdml_cuda_device_properties *properties)
{
    (void) device;
    if (properties != NULL) {
        memset(properties, 0, sizeof(*properties));
    }
    return KDML_CUDA_UNAVAILABLE;
}

const char *kdml_cuda_last_error(void)
{
    return "kdml was installed without CUDA support";
}

const char *kdml_cuda_status_string(int status)
{
    switch (status) {
    case KDML_CUDA_SUCCESS:
        return "success";
    case KDML_CUDA_UNAVAILABLE:
        return "CUDA unavailable";
    case KDML_CUDA_INVALID_ARGUMENT:
        return "invalid CUDA argument";
    case KDML_CUDA_ALLOCATION_ERROR:
        return "CUDA allocation error";
    case KDML_CUDA_RUNTIME_ERROR:
        return "CUDA runtime error";
    case KDML_CUDA_INVALID_STATE:
        return "invalid CUDA state";
    default:
        return "unknown CUDA status";
    }
}

int kdml_cuda_create(int n, int p, int metric,
                     const double *x, const int *type,
                     const int *level_count, int continuous_count,
                     int requested_device,
                     kdml_cuda_workspace **workspace)
{
    (void) n;
    (void) p;
    (void) metric;
    (void) x;
    (void) type;
    (void) level_count;
    (void) continuous_count;
    (void) requested_device;
    if (workspace != NULL) {
        *workspace = NULL;
    }
    return KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_build(kdml_cuda_workspace *workspace,
                    const int *kernel_code, const double *bandwidth,
                    double *score)
{
    (void) workspace;
    (void) kernel_code;
    (void) bandwidth;
    (void) score;
    return KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_propose(kdml_cuda_workspace *workspace,
                      int feature, int proposed_kernel,
                      double proposed_bandwidth, double *score)
{
    (void) workspace;
    (void) feature;
    (void) proposed_kernel;
    (void) proposed_bandwidth;
    (void) score;
    return KDML_CUDA_UNAVAILABLE;
}

int kdml_cuda_accept(kdml_cuda_workspace *workspace,
                     int feature, int proposed_kernel,
                     double proposed_bandwidth)
{
    (void) workspace;
    (void) feature;
    (void) proposed_kernel;
    (void) proposed_bandwidth;
    return KDML_CUDA_UNAVAILABLE;
}

void kdml_cuda_destroy(kdml_cuda_workspace *workspace)
{
    (void) workspace;
}

#endif
