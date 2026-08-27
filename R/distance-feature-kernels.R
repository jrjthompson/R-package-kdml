.kdml_validate_flag <- function(value, label) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", label, "` must be TRUE or FALSE.", call. = FALSE)
  }
  as.logical(value)
}

.kdml_validate_scalar_kernel <- function(value, choices, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% choices) {
    stop("Invalid ", label, " specified. Choose one of: ",
         paste(choices, collapse = ", "), call. = FALSE)
  }
  value
}

.kdml_fixed_distance <- function(df, distance, bandwidths, kernels,
                                 standardize, drop_unused = TRUE) {
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    stop("`stan` must be TRUE or FALSE.", call. = FALSE)
  }
  prep <- .kdml_prepare_data(df, drop_unused = drop_unused)
  state <- .kdml_validate_fixed_state(kernels, bandwidths, prep)
  spec <- list(
    x = prep$x,
    type = unname(prep$type),
    level_count = unname(prep$level_count),
    kernel = unname(state$kernel_code),
    bandwidth = unname(state$bandwidth),
    metric = if (distance == "dkps") 0L else 1L
  )
  distances <- .Call(C_kdml_distance_call, spec)
  dimnames(distances) <- list(prep$row_names, prep$row_names)

  tolerance <- 100 * .Machine$double.eps * max(1, max(abs(distances)))
  if (any(distances < -tolerance)) {
    warning("The selected kernel state produced negative squared distances; ",
            "negative values were truncated to zero.", call. = FALSE)
  }
  distances[distances < 0] <- 0
  if (isTRUE(standardize)) {
    maximum <- max(distances)
    if (is.finite(maximum) && maximum > 0) {
      distances <- distances / maximum
    }
  }

  list(
    distances = distances,
    bandwidths = matrix(
      state$bandwidth, nrow = 1L,
      dimnames = list(NULL, prep$feature_names)
    ),
    kernels = state$kernel,
    standardized = isTRUE(standardize)
  )
}
