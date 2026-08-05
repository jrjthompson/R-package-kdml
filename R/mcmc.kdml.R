.kdml_kernel_catalogue <- function() {
  list(
    continuous = c(
      c_gaussian = 0L,
      c_epanechnikov = 1L,
      c_uniform = 2L,
      c_triangle = 3L,
      c_biweight = 4L,
      c_triweight = 5L,
      c_tricube = 6L,
      c_cosine = 7L,
      c_logistic = 8L,
      c_sigmoid = 9L,
      c_silverman = 10L
    ),
    nominal = c(
      u_aitken = 0L,
      u_aitchisonaitken = 1L
    ),
    ordinal = c(
      o_wangvanryzin = 0L,
      o_habbema = 1L,
      o_aitken = 2L,
      o_aitchisonaitken = 3L,
      o_liracine = 4L
    )
  )
}

.kdml_type_name <- function(type) {
  c("continuous", "nominal", "ordinal")[type + 1L]
}

.kdml_canonical_kernel_name <- function(x) {
  x[x == "u.aitchisonaitken"] <- "u_aitchisonaitken"
  x
}

.kdml_prepare_data <- function(df) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }
  if (nrow(df) < 2L) {
    stop("`df` must contain at least two observations.", call. = FALSE)
  }
  if (ncol(df) < 1L) {
    stop("`df` must contain at least one feature.", call. = FALSE)
  }

  feature_names <- names(df)
  if (is.null(feature_names) || anyNA(feature_names) ||
      any(!nzchar(feature_names)) || anyDuplicated(feature_names)) {
    stop("Every feature must have a unique, non-empty name.", call. = FALSE)
  }

  n <- nrow(df)
  p <- ncol(df)
  x <- matrix(0, nrow = n, ncol = p,
              dimnames = list(rownames(df), feature_names))
  type <- integer(p)
  level_count <- integer(p)
  scales <- rep(1, p)
  ranges <- numeric(p)
  categorical_levels <- vector("list", p)

  for (j in seq_len(p)) {
    column <- df[[j]]

    if (is.ordered(column)) {
      if (anyNA(column)) {
        stop("Missing values are not supported; found one in `",
             feature_names[j], "`.", call. = FALSE)
      }
      column <- droplevels(column)
      observed_levels <- levels(column)
      if (length(observed_levels) < 2L) {
        stop("Ordinal feature `", feature_names[j],
             "` must have at least two observed levels.", call. = FALSE)
      }
      x[, j] <- as.integer(column) - 1L
      type[j] <- 2L
      level_count[j] <- length(observed_levels)
      categorical_levels[[j]] <- observed_levels
      ranges[j] <- length(observed_levels) - 1L
    } else if (is.factor(column)) {
      if (anyNA(column)) {
        stop("Missing values are not supported; found one in `",
             feature_names[j], "`.", call. = FALSE)
      }
      column <- droplevels(column)
      observed_levels <- levels(column)
      if (length(observed_levels) < 2L) {
        stop("Nominal feature `", feature_names[j],
             "` must have at least two observed levels.", call. = FALSE)
      }
      x[, j] <- as.integer(column) - 1L
      type[j] <- 1L
      level_count[j] <- length(observed_levels)
      categorical_levels[[j]] <- observed_levels
      ranges[j] <- length(observed_levels) - 1L
    } else if (is.numeric(column)) {
      column <- as.numeric(column)
      if (any(!is.finite(column))) {
        stop("Continuous feature `", feature_names[j],
             "` contains a missing or non-finite value.", call. = FALSE)
      }
      scale_j <- stats::sd(column)
      range_j <- diff(range(column))
      if (!is.finite(scale_j) || scale_j <= 0 || !is.finite(range_j) || range_j <= 0) {
        stop("Continuous feature `", feature_names[j],
             "` must vary across observations.", call. = FALSE)
      }
      x[, j] <- column
      type[j] <- 0L
      scales[j] <- scale_j
      ranges[j] <- range_j
    } else {
      stop("Feature `", feature_names[j],
           "` must be numeric, a factor, or an ordered factor.", call. = FALSE)
    }
  }

  storage.mode(x) <- "double"
  names(categorical_levels) <- feature_names
  names(scales) <- names(ranges) <- names(type) <- names(level_count) <- feature_names

  list(
    x = x,
    type = type,
    type_name = .kdml_type_name(type),
    level_count = level_count,
    scale = scales,
    range = ranges,
    levels = categorical_levels,
    feature_names = feature_names,
    row_names = rownames(df)
  )
}

.kdml_prepare_new_data <- function(df, training) {
  if (!is.data.frame(df)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (nrow(df) < 1L) {
    stop("`data` must contain at least one observation.", call. = FALSE)
  }
  if (is.null(names(df)) || anyDuplicated(names(df)) ||
      !setequal(names(df), training$feature_names)) {
    stop("`data` must contain exactly the features used to fit the sampler.",
         call. = FALSE)
  }

  df <- df[, training$feature_names, drop = FALSE]
  x <- matrix(0, nrow = nrow(df), ncol = length(training$feature_names),
              dimnames = list(rownames(df), training$feature_names))

  for (j in seq_along(training$feature_names)) {
    feature <- training$feature_names[j]
    column <- df[[j]]
    type <- training$type[j]

    if (type == 0L) {
      if (!is.numeric(column)) {
        stop("Continuous feature `", feature, "` must remain numeric.",
             call. = FALSE)
      }
      column <- as.numeric(column)
      if (any(!is.finite(column))) {
        stop("Feature `", feature, "` contains a missing or non-finite value.",
             call. = FALSE)
      }
      x[, j] <- column
    } else {
      if ((type == 1L && (!is.factor(column) || is.ordered(column))) ||
          (type == 2L && !is.ordered(column))) {
        stop("Categorical class for feature `", feature,
             "` does not match the fitted data.", call. = FALSE)
      }
      encoded <- match(as.character(column), training$levels[[j]])
      if (anyNA(encoded)) {
        stop("Feature `", feature,
             "` contains a level not observed during fitting.", call. = FALSE)
      }
      x[, j] <- encoded - 1L
    }
  }

  storage.mode(x) <- "double"
  x
}

.kdml_normalize_candidates <- function(kernel_candidates, prep) {
  catalogue <- .kdml_kernel_catalogue()
  p <- length(prep$type)
  candidates <- lapply(prep$type_name, function(type) names(catalogue[[type]]))
  names(candidates) <- prep$feature_names

  if (!is.null(kernel_candidates)) {
    if (!is.list(kernel_candidates)) {
      stop("`kernel_candidates` must be NULL or a list of character vectors.",
           call. = FALSE)
    }
    if (is.null(names(kernel_candidates))) {
      if (length(kernel_candidates) != p) {
        stop("An unnamed `kernel_candidates` list must have one element per feature.",
             call. = FALSE)
      }
      names(kernel_candidates) <- prep$feature_names
    } else {
      if (any(!nzchar(names(kernel_candidates))) || anyDuplicated(names(kernel_candidates)) ||
          any(!names(kernel_candidates) %in% prep$feature_names)) {
        stop("`kernel_candidates` contains invalid or duplicated feature names.",
             call. = FALSE)
      }
    }
    for (feature in names(kernel_candidates)) {
      candidates[[feature]] <- kernel_candidates[[feature]]
    }
  }

  codes <- vector("list", p)
  for (j in seq_len(p)) {
    value <- candidates[[j]]
    if (!is.character(value) || length(value) < 1L || anyNA(value)) {
      stop("Kernel candidates for `", prep$feature_names[j],
           "` must be a non-empty character vector.", call. = FALSE)
    }
    value <- unique(.kdml_canonical_kernel_name(value))
    allowed <- catalogue[[prep$type_name[j]]]
    if (any(!value %in% names(allowed))) {
      invalid <- value[!value %in% names(allowed)]
      stop("Kernel `", invalid[1L], "` is not valid for ",
           prep$type_name[j], " feature `", prep$feature_names[j], "`.",
           call. = FALSE)
    }
    candidates[[j]] <- value
    codes[[j]] <- unname(allowed[value])
  }

  offsets <- integer(p + 1L)
  offsets[-1L] <- cumsum(lengths(codes))
  list(
    names = candidates,
    codes = as.integer(unlist(codes, use.names = FALSE)),
    offsets = offsets
  )
}

.kdml_normalize_initial_kernels <- function(initial_kernels, candidates, prep) {
  defaults <- c(
    continuous = "c_gaussian",
    nominal = "u_aitken",
    ordinal = "o_wangvanryzin"
  )
  answer <- unname(defaults[prep$type_name])
  names(answer) <- prep$feature_names

  for (j in seq_along(answer)) {
    if (!answer[j] %in% candidates$names[[j]]) {
      answer[j] <- candidates$names[[j]][1L]
    }
  }

  if (!is.null(initial_kernels)) {
    if (!is.character(initial_kernels) || anyNA(initial_kernels)) {
      stop("`initial_kernels` must be a character vector.", call. = FALSE)
    }
    if (is.null(names(initial_kernels))) {
      if (length(initial_kernels) != length(answer)) {
        stop("Unnamed `initial_kernels` must have one value per feature.",
             call. = FALSE)
      }
      names(initial_kernels) <- prep$feature_names
    } else if (any(!names(initial_kernels) %in% prep$feature_names) ||
               anyDuplicated(names(initial_kernels))) {
      stop("`initial_kernels` contains invalid or duplicated feature names.",
           call. = FALSE)
    }
    answer[names(initial_kernels)] <- initial_kernels
  }

  answer <- .kdml_canonical_kernel_name(answer)
  catalogue <- .kdml_kernel_catalogue()
  codes <- integer(length(answer))
  for (j in seq_along(answer)) {
    if (!answer[j] %in% candidates$names[[j]]) {
      stop("Initial kernel `", answer[j], "` is not a candidate for feature `",
           prep$feature_names[j], "`.", call. = FALSE)
    }
    codes[j] <- unname(catalogue[[prep$type_name[j]]][answer[j]])
  }
  list(names = answer, codes = codes)
}

.kdml_expand_numeric <- function(x, prep, label, positive = FALSE) {
  if (!is.numeric(x) || anyNA(x)) {
    stop("`", label, "` must be numeric.", call. = FALSE)
  }
  if (length(x) == 1L) {
    x <- rep(as.numeric(x), length(prep$type))
  } else {
    if (!is.null(names(x))) {
      if (anyDuplicated(names(x)) ||
          !setequal(names(x), prep$feature_names)) {
        stop("Named `", label, "` must contain every feature exactly once.",
             call. = FALSE)
      }
      x <- x[prep$feature_names]
    } else if (length(x) != length(prep$type)) {
      stop("`", label, "` must be scalar or have one value per feature.",
           call. = FALSE)
    }
    x <- as.numeric(x)
  }
  if (any(!is.finite(x)) || (positive && any(x <= 0))) {
    qualifier <- if (positive) "finite and strictly positive" else "finite"
    stop("Every `", label, "` value must be ", qualifier, ".", call. = FALSE)
  }
  names(x) <- prep$feature_names
  x
}

.kdml_bandwidth_upper <- function(type, kernel_code, level_count) {
  if (type == 0L) {
    return(Inf)
  }
  if (type == 1L && kernel_code == 1L) {
    return((level_count - 1) / level_count)
  }
  1
}

.kdml_default_bandwidths <- function(prep, initial_codes) {
  answer <- numeric(length(prep$type))
  for (j in seq_along(answer)) {
    if (prep$type[j] == 0L) {
      answer[j] <- prep$scale[j]
      if (initial_codes[j] %in% c(1L:7L, 10L)) {
        answer[j] <- max(answer[j], 1.05 * prep$range[j])
      }
    } else {
      answer[j] <- 0.5 * .kdml_bandwidth_upper(
        prep$type[j], initial_codes[j], prep$level_count[j]
      )
    }
  }
  names(answer) <- prep$feature_names
  answer
}

.kdml_bandwidth_to_theta <- function(bandwidth, prep, kernel_codes) {
  theta <- numeric(length(bandwidth))
  for (j in seq_along(theta)) {
    if (prep$type[j] == 0L) {
      if (bandwidth[j] <= 0) {
        stop("Continuous bandwidths must be strictly positive.", call. = FALSE)
      }
      theta[j] <- log(bandwidth[j] / prep$scale[j])
    } else {
      upper <- .kdml_bandwidth_upper(
        prep$type[j], kernel_codes[j], prep$level_count[j]
      )
      if (bandwidth[j] <= 0 || bandwidth[j] >= upper) {
        stop("Initial categorical bandwidth for `", prep$feature_names[j],
             "` must lie strictly between 0 and ", signif(upper, 8), ".",
             call. = FALSE)
      }
      theta[j] <- stats::qlogis(bandwidth[j] / upper)
    }
  }
  if (any(!is.finite(theta))) {
    stop("Initial bandwidths produced a non-finite transformed state.",
         call. = FALSE)
  }
  names(theta) <- prep$feature_names
  theta
}

.kdml_scalar_integer <- function(x, label, minimum) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop("`", label, "` must be an integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.kdml_kernel_names_from_codes <- function(codes, prep) {
  catalogue <- .kdml_kernel_catalogue()
  answer <- character(length(codes))
  for (j in seq_along(answer)) {
    choices <- catalogue[[prep$type_name[j]]]
    answer[j] <- names(choices)[match(codes[j], unname(choices))]
  }
  names(answer) <- prep$feature_names
  answer
}

#'
#' @export
mcmc.kdml <- function(df, distance = c("dkps", "dkss"), chains = 1L,
                      chain_threads = 1L,
                      warmup = 1000L, draws = 1000L, thin = 1L,
                      beta = 1, proposal_sd = 0.15,
                      prior_mean = 0, prior_sd = 4,
                      adapt = TRUE, target_accept = 0.44,
                      kernel_candidates = NULL, initial_kernels = NULL,
                      initial_bandwidths = NULL,
                      backend = c("cpu", "cuda", "auto"),
                      progress = interactive()) {
  call <- match.call()
  distance <- match.arg(distance)
  backend <- match.arg(backend)
  prep <- .kdml_prepare_data(df)

  chains <- .kdml_scalar_integer(chains, "chains", 1L)
  chain_threads <- .kdml_scalar_integer(
    chain_threads, "chain_threads", 1L
  )
  chain_threads <- min(chain_threads, chains)
  if (chain_threads > 1L && backend != "cpu") {
    stop("`chain_threads > 1` requires `backend = \"cpu\"`.",
         call. = FALSE)
  }
  warmup <- .kdml_scalar_integer(warmup, "warmup", 0L)
  draws <- .kdml_scalar_integer(draws, "draws", 1L)
  thin <- .kdml_scalar_integer(thin, "thin", 1L)
  if (draws > floor((.Machine$integer.max - warmup) / thin)) {
    stop("The requested warmup, draws, and thinning exceed the supported iteration count.",
         call. = FALSE)
  }
  if (length(beta) != 1L || !is.numeric(beta) || !is.finite(beta) || beta <= 0) {
    stop("`beta` must be finite and strictly positive.", call. = FALSE)
  }
  if (!is.logical(adapt) || length(adapt) != 1L || is.na(adapt)) {
    stop("`adapt` must be TRUE or FALSE.", call. = FALSE)
  }
  adapt <- as.logical(adapt)
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  progress <- as.logical(progress)
  if (length(target_accept) != 1L || !is.numeric(target_accept) ||
      !is.finite(target_accept) || target_accept <= 0 || target_accept >= 1) {
    stop("`target_accept` must lie strictly between 0 and 1.", call. = FALSE)
  }

  proposal_sd <- .kdml_expand_numeric(
    proposal_sd, prep, "proposal_sd", positive = TRUE
  )
  prior_mean <- .kdml_expand_numeric(prior_mean, prep, "prior_mean")
  prior_sd <- .kdml_expand_numeric(prior_sd, prep, "prior_sd", positive = TRUE)

  candidates <- .kdml_normalize_candidates(kernel_candidates, prep)
  initial <- .kdml_normalize_initial_kernels(initial_kernels, candidates, prep)

  if (is.null(initial_bandwidths)) {
    initial_bandwidths <- .kdml_default_bandwidths(prep, initial$codes)
  } else {
    initial_bandwidths <- .kdml_expand_numeric(
      initial_bandwidths, prep, "initial_bandwidths", positive = TRUE
    )
  }
  initial_theta <- .kdml_bandwidth_to_theta(
    initial_bandwidths, prep, initial$codes
  )

  metric_code <- if (distance == "dkps") 0L else 1L
  spec <- list(
    x = prep$x,
    type = unname(prep$type),
    level_count = unname(prep$level_count),
    scale = unname(prep$scale),
    candidate_codes = candidates$codes,
    candidate_offsets = candidates$offsets,
    initial_kernel = unname(initial$codes),
    initial_theta = unname(initial_theta),
    proposal_sd = unname(proposal_sd),
    prior_mean = unname(prior_mean),
    prior_sd = unname(prior_sd),
    metric = metric_code,
    chains = chains,
    chain_threads = chain_threads,
    warmup = warmup,
    draws = draws,
    thin = thin,
    beta = as.numeric(beta),
    adapt = adapt,
    target_accept = as.numeric(target_accept),
    backend = match(backend, c("cpu", "cuda", "auto")) - 1L,
    progress = progress,
    progress_in_place = progress && interactive()
  )

  raw <- .Call(C_kdml_mcmc_call, spec)
  draw_names <- as.character(seq_len(draws))
  chain_names <- paste0("chain", seq_len(chains))
  draw_dimnames <- list(
    draw = draw_names,
    feature = prep$feature_names,
    chain = chain_names
  )
  dimnames(raw$theta) <- draw_dimnames
  dimnames(raw$bandwidth) <- draw_dimnames
  dimnames(raw$kernel_code) <- draw_dimnames
  dimnames(raw$score) <- list(draw = draw_names, chain = chain_names)
  dimnames(raw$log_target) <- list(draw = draw_names, chain = chain_names)

  diagnostic_names <- list(feature = prep$feature_names, chain = chain_names)
  diagnostic_fields <- c(
    "proposal_sd", "kernel_attempts", "kernel_accepts",
    "bandwidth_attempts", "bandwidth_accepts",
    "sampling_kernel_attempts", "sampling_kernel_accepts",
    "sampling_bandwidth_attempts", "sampling_bandwidth_accepts"
  )
  for (field in diagnostic_fields) {
    dimnames(raw[[field]]) <- diagnostic_names
  }
  names(raw$invalid_proposals) <- chain_names

  kernel_names <- array(NA_character_, dim = dim(raw$kernel_code),
                        dimnames = draw_dimnames)
  catalogue <- .kdml_kernel_catalogue()
  for (j in seq_along(prep$type)) {
    choices <- catalogue[[prep$type_name[j]]]
    kernel_names[, j, ] <- names(choices)[
      match(raw$kernel_code[, j, ], unname(choices))
    ]
  }

  map_index <- which.max(raw$log_target)
  map_location <- arrayInd(map_index, dim(raw$log_target))
  map_draw <- map_location[1L]
  map_chain <- map_location[2L]
  map_codes <- raw$kernel_code[map_draw, , map_chain]
  names(map_codes) <- prep$feature_names

  fit <- structure(
    list(
      call = call,
      distance = distance,
      backend = raw$backend,
      theta = raw$theta,
      bandwidth = raw$bandwidth,
      kernel = kernel_names,
      kernel_code = raw$kernel_code,
      score = raw$score,
      log_target = raw$log_target,
      map = list(
        draw = map_draw,
        chain = map_chain,
        theta = setNames(raw$theta[map_draw, , map_chain], prep$feature_names),
        bandwidth = setNames(raw$bandwidth[map_draw, , map_chain], prep$feature_names),
        kernel = .kdml_kernel_names_from_codes(map_codes, prep),
        kernel_code = map_codes,
        score = raw$score[map_draw, map_chain],
        log_target = raw$log_target[map_draw, map_chain]
      ),
      diagnostics = raw[diagnostic_fields],
      invalid_proposals = raw$invalid_proposals,
      candidates = candidates$names,
      initial = list(
        kernel = initial$names,
        kernel_code = setNames(initial$codes, prep$feature_names),
        bandwidth = initial_bandwidths,
        theta = initial_theta
      ),
      control = list(
        chains = chains,
        chain_threads = chain_threads,
        warmup = warmup,
        draws = draws,
        thin = thin,
        beta = beta,
        prior_mean = prior_mean,
        prior_sd = prior_sd,
        adapt = adapt,
        target_accept = target_accept,
        backend = backend,
        progress = progress
      ),
      preprocessing = prep,
      target_version = "kdml-mscv-v1"
    ),
    class = "kdml_mcmc"
  )
  fit$draws <- as_draws_array.kdml_mcmc(
    fit,
    variables = c(
      "theta", "bandwidth", "kernel_code", "score", "log_target"
    )
  )
  fit
}

#' @rdname mcmc.kdml
#' @export
mcmc.dkps <- function(df, ...) {
  mcmc.kdml(df = df, distance = "dkps", ...)
}

#' @rdname mcmc.kdml
#' @export
mcmc.dkss <- function(df, ...) {
  mcmc.kdml(df = df, distance = "dkss", ...)
}

#'
as_draws_array.kdml_mcmc <- function(
    x,
    variables = c("theta", "bandwidth", "score", "log_target"),
    ...) {
  if (!inherits(x, "kdml_mcmc")) {
    stop("`x` must be a `kdml_mcmc` object.", call. = FALSE)
  }
  choices <- c("theta", "bandwidth", "kernel_code", "score", "log_target")
  variables <- unique(match.arg(variables, choices, several.ok = TRUE))

  if (!is.null(x$draws) &&
      posterior::is_draws_array(x$draws)) {
    retained <- posterior::variables(x$draws)
    selected <- unlist(lapply(variables, function(group) {
      retained[
        retained == group |
          startsWith(retained, paste0(group, "["))
      ]
    }), use.names = FALSE)
    if (!length(selected)) {
      stop("The stored draws do not contain the requested variables.",
           call. = FALSE)
    }
    return(posterior::subset_draws(x$draws, variable = selected))
  }

  # Compatibility path for kdml_mcmc objects saved before the canonical
  # posterior draws component was introduced.
  draw_count <- dim(x$theta)[1L]
  chain_count <- dim(x$theta)[3L]
  iteration_names <- dimnames(x$theta)[[1L]]
  chain_names <- dimnames(x$theta)[[3L]]
  feature_names <- dimnames(x$theta)[[2L]]
  pieces <- vector("list", length(variables))

  for (index in seq_along(variables)) {
    group <- variables[index]
    if (group %in% c("theta", "bandwidth", "kernel_code")) {
      value <- aperm(x[[group]], c(1L, 3L, 2L))
      storage.mode(value) <- "double"
      dimnames(value) <- list(
        iteration = iteration_names,
        chain = chain_names,
        variable = paste0(group, "[", feature_names, "]")
      )
    } else {
      value <- array(
        as.numeric(x[[group]]),
        dim = c(draw_count, chain_count, 1L),
        dimnames = list(
          iteration = iteration_names,
          chain = chain_names,
          variable = group
        )
      )
    }
    pieces[[index]] <- value
  }

  variable_count <- sum(vapply(pieces, function(value) dim(value)[3L],
                               integer(1)))
  combined <- array(
    NA_real_,
    dim = c(draw_count, chain_count, variable_count),
    dimnames = list(
      iteration = iteration_names,
      chain = chain_names,
      variable = unlist(
        lapply(pieces, function(value) dimnames(value)[[3L]]),
        use.names = FALSE
      )
    )
  )
  offset <- 0L
  for (value in pieces) {
    width <- dim(value)[3L]
    combined[, , offset + seq_len(width)] <- value
    offset <- offset + width
  }

  posterior::as_draws_array(combined)
}

#'
#' @export
kdml.diagnostics <- function(
    object,
    variables = "theta") {
  if (!inherits(object, "kdml_mcmc")) {
    stop("`object` must be a `kdml_mcmc` object.", call. = FALSE)
  }
  draws <- posterior::as_draws_array(object, variables = variables)
  answer <- posterior::summarise_draws(
    draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"
  )
  answer <- as.data.frame(answer, stringsAsFactors = FALSE)
  rownames(answer) <- NULL
  answer
}

#'
#' @export
cuda_available <- function() {
  .Call(C_kdml_cuda_available_call)
}

#' @rdname cuda_available
#' @export
cuda_info <- function() {
  .Call(C_kdml_cuda_info_call)
}

#' @export
print.kdml_mcmc <- function(x, ...) {
  chain_threads <- if (is.null(x$control$chain_threads)) {
    1L
  } else {
    x$control$chain_threads
  }
  cat("KDML MCMC fit\n")
  cat("  distance: ", toupper(x$distance), "\n", sep = "")
  cat("  backend: ", x$backend, "\n", sep = "")
  cat("  features: ", length(x$preprocessing$feature_names), "\n", sep = "")
  cat("  chains: ", x$control$chains,
      ", chain threads: ", chain_threads,
      ", retained draws per chain: ", x$control$draws, "\n", sep = "")
  cat("  retained MAP: draw ", x$map$draw,
      ", chain ", x$map$chain, "\n", sep = "")
  cat("  MAP log target: ",
      format(x$map$log_target, digits = 6), "\n", sep = "")
  cat("  MAP MSCV score: ",
      format(x$map$score, digits = 6), "\n", sep = "")
  cat("  MAP kernels:\n")
  print(x$map$kernel)
  invisible(x)
}

.kdml_validate_kernels <- function(kernels, prep) {
  if (!is.character(kernels) || anyNA(kernels)) {
    stop("`kernels` must be a character vector.", call. = FALSE)
  }
  if (!is.null(names(kernels))) {
    if (anyDuplicated(names(kernels)) ||
        !setequal(names(kernels), prep$feature_names)) {
      stop("Named `kernels` must contain every feature exactly once.",
           call. = FALSE)
    }
    kernels <- kernels[prep$feature_names]
  } else if (length(kernels) != length(prep$type)) {
    stop("`kernels` must have one value per feature.", call. = FALSE)
  }
  kernels <- .kdml_canonical_kernel_name(kernels)
  catalogue <- .kdml_kernel_catalogue()
  codes <- integer(length(kernels))
  for (j in seq_along(kernels)) {
    allowed <- catalogue[[prep$type_name[j]]]
    if (!kernels[j] %in% names(allowed)) {
      stop("Kernel `", kernels[j], "` is invalid for feature `",
           prep$feature_names[j], "`.", call. = FALSE)
    }
    codes[j] <- unname(allowed[kernels[j]])
  }
  list(kernel = kernels, kernel_code = codes)
}

.kdml_validate_fixed_state <- function(kernels, bandwidths, prep) {
  kernel_state <- .kdml_validate_kernels(kernels, prep)
  kernels <- kernel_state$kernel
  codes <- kernel_state$kernel_code
  bandwidths <- .kdml_expand_numeric(
    bandwidths, prep, "bandwidths", positive = TRUE
  )
  for (j in seq_along(bandwidths)) {
    upper <- .kdml_bandwidth_upper(prep$type[j], codes[j], prep$level_count[j])
    if (bandwidths[j] >= upper) {
      stop("Bandwidth for `", prep$feature_names[j],
           "` must be below ", signif(upper, 8), ".", call. = FALSE)
    }
  }
  list(kernel = kernels, kernel_code = codes, bandwidth = bandwidths)
}

.kdml_mcmc_evaluate <- function(df, distance = c("dkps", "dkss"),
                                kernels, bandwidths,
                                return_similarity = FALSE,
                                backend = c("cpu", "cuda", "auto")) {
  distance <- match.arg(distance)
  backend <- match.arg(backend)
  prep <- .kdml_prepare_data(df)
  state <- .kdml_validate_fixed_state(kernels, bandwidths, prep)
  spec <- list(
    x = prep$x,
    type = unname(prep$type),
    level_count = unname(prep$level_count),
    kernel = unname(state$kernel_code),
    bandwidth = unname(state$bandwidth),
    metric = if (distance == "dkps") 0L else 1L,
    return_similarity = isTRUE(return_similarity),
    backend = match(backend, c("cpu", "cuda", "auto")) - 1L
  )
  answer <- .Call(C_kdml_score_call, spec)
  if (!is.null(answer$similarity)) {
    dimnames(answer$similarity) <- list(prep$row_names, prep$row_names)
  }
  answer
}

#'
#' @export
kdml.distance <- function(object, data = NULL, state = c("map", "last"),
                          chain = 1L, draw = NULL, standardize = TRUE) {
  if (!inherits(object, "kdml_mcmc")) {
    stop("`object` must inherit from `kdml_mcmc`.", call. = FALSE)
  }
  state <- match.arg(state)
  chain <- .kdml_scalar_integer(chain, "chain", 1L)
  if (chain > object$control$chains) {
    stop("`chain` exceeds the number of fitted chains.", call. = FALSE)
  }

  if (is.null(draw)) {
    if (state == "map") {
      draw <- object$map$draw
      chain <- object$map$chain
    } else {
      draw <- object$control$draws
    }
  } else {
    draw <- .kdml_scalar_integer(draw, "draw", 1L)
    if (draw > object$control$draws) {
      stop("`draw` exceeds the number of retained draws.", call. = FALSE)
    }
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    stop("`standardize` must be TRUE or FALSE.", call. = FALSE)
  }

  prep <- object$preprocessing
  if (is.null(data)) {
    x <- prep$x
    row_names <- prep$row_names
  } else {
    x <- .kdml_prepare_new_data(data, prep)
    row_names <- rownames(data)
  }

  codes <- as.integer(object$kernel_code[draw, , chain])
  bandwidths <- as.numeric(object$bandwidth[draw, , chain])
  spec <- list(
    x = x,
    type = unname(prep$type),
    level_count = unname(prep$level_count),
    kernel = codes,
    bandwidth = bandwidths,
    metric = if (object$distance == "dkps") 0L else 1L
  )
  distances <- .Call(C_kdml_distance_call, spec)
  dimnames(distances) <- list(row_names, row_names)

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

  kernel_names <- .kdml_kernel_names_from_codes(codes, prep)
  bandwidth_matrix <- matrix(
    bandwidths, nrow = 1L,
    dimnames = list(NULL, prep$feature_names)
  )
  list(
    distances = distances,
    bandwidths = bandwidth_matrix,
    kernels = kernel_names,
    state = list(draw = draw, chain = chain),
    squared = TRUE,
    standardized = isTRUE(standardize)
  )
}
