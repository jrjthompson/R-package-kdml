main_library <- paste0("kdml", SHLIB_EXT)
cuda_marker <- "kdml_cuda_enabled"
cuda_library <- if (WINDOWS) "kdmlcuda.dll" else "kdmlcuda.so"
files <- main_library
cuda_enabled <- file.exists(cuda_marker)
if (cuda_enabled) {
    files <- c(files, cuda_library)
}

missing <- files[!file.exists(files)]
if (length(missing)) {
    stop("Required kdml shared libraries were not produced: ",
         paste(missing, collapse = ", "))
}

destination <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
if (!cuda_enabled) {
    stale_sidecar <- file.path(destination, cuda_library)
    if (file.exists(stale_sidecar) && !file.remove(stale_sidecar)) {
        stop("Failed to remove a stale CUDA sidecar: ", stale_sidecar)
    }
}
copied <- file.copy(files, destination, overwrite = TRUE)
if (!all(copied)) {
    stop("Failed to install one or more kdml shared libraries: ",
         paste(files[!copied], collapse = ", "))
}

if (file.exists("symbols.rds")) {
    if (!file.copy("symbols.rds", destination, overwrite = TRUE)) {
        stop("Failed to install symbols.rds.")
    }
}
