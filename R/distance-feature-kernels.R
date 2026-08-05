.kdml_expand_type_kernels <- function(value, features, label) {
  if (length(features) == 0L) {
    return(setNames(character(), character()))
  }
  if (!is.character(value) || anyNA(value) || any(!nzchar(value))) {
    stop("`", label, "` must be a character vector.", call. = FALSE)
  }
  if (length(value) == 1L && is.null(names(value))) {
    return(setNames(rep(value, length(features)), features))
  }
  if (!is.null(names(value))) {
    if (any(!nzchar(names(value))) || anyDuplicated(names(value)) ||
        !setequal(names(value), features)) {
      stop("Named `", label, "` must contain every ", label,
           " feature exactly once.", call. = FALSE)
    }
    return(value[features])
  }
  if (length(value) != length(features)) {
    stop("`", label, "` must be scalar or have one value per applicable feature.",
         call. = FALSE)
  }
  setNames(value, features)
}

.kdml_feature_kernels <- function(df, kernels, cFUN, uFUN, oFUN) {
  prep <- .kdml_prepare_data(df)
  if (!is.null(kernels)) {
    return(.kdml_validate_kernels(kernels, prep)$kernel)
  }

  continuous <- prep$feature_names[prep$type == 0L]
  nominal <- prep$feature_names[prep$type == 1L]
  ordinal <- prep$feature_names[prep$type == 2L]
  answer <- c(
    .kdml_expand_type_kernels(cFUN, continuous, "cFUN"),
    .kdml_expand_type_kernels(uFUN, nominal, "uFUN"),
    .kdml_expand_type_kernels(oFUN, ordinal, "oFUN")
  )
  answer[prep$feature_names]
}

.kdml_fixed_distance <- function(df, distance, bandwidths, kernels,
                                 standardize) {
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    stop("`stan` must be TRUE or FALSE.", call. = FALSE)
  }
  prep <- .kdml_prepare_data(df)
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
    squared = TRUE,
    standardized = isTRUE(standardize)
  )
}
