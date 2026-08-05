.kdml_mcmc_dimensions <- function(object) {
  if (!inherits(object, "kdml_mcmc")) {
    stop("`object` must be a `kdml_mcmc` object.", call. = FALSE)
  }
  dimensions <- dim(object$theta)
  if (length(dimensions) != 3L || any(!is.finite(dimensions)) ||
      any(dimensions < 1L)) {
    stop("The `kdml_mcmc` object has malformed `theta` draws.",
         call. = FALSE)
  }

  draw_names <- dimnames(object$theta)[[1L]]
  feature_names <- dimnames(object$theta)[[2L]]
  chain_names <- dimnames(object$theta)[[3L]]
  if (is.null(draw_names)) {
    draw_names <- as.character(seq_len(dimensions[1L]))
  }
  if (is.null(feature_names)) {
    feature_names <- paste0("feature", seq_len(dimensions[2L]))
  }
  if (is.null(chain_names)) {
    chain_names <- paste0("chain", seq_len(dimensions[3L]))
  }

  list(
    draws = dimensions[1L],
    features = dimensions[2L],
    chains = dimensions[3L],
    draw_names = draw_names,
    feature_names = feature_names,
    chain_names = chain_names
  )
}

.kdml_score_summary <- function(object) {
  fields <- c("score", "log_target")
  labels <- c("MSCV score", "log target")
  draws <- posterior::as_draws_array(object, variables = fields)
  summaries <- suppressWarnings(posterior::summarise_draws(
    draws, "mean", "sd", "min", "max"
  ))
  summaries <- as.data.frame(summaries, stringsAsFactors = FALSE)
  answer <- lapply(seq_along(fields), function(index) {
    field <- fields[index]
    values <- summaries[summaries$variable == field, , drop = FALSE]
    map <- if (!is.null(object$map[[field]])) {
      as.numeric(object$map[[field]])[1L]
    } else {
      NA_real_
    }
    data.frame(
      quantity = labels[index],
      mean = values$mean,
      sd = values$sd,
      min = values$min,
      max = values$max,
      map_state = map,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, answer)
}

.kdml_kernel_diagnostics <- function(object, dimensions) {
  expected <- c(dimensions$draws, dimensions$features, dimensions$chains)
  if (is.null(object$kernel) ||
      !identical(as.integer(dim(object$kernel)), as.integer(expected))) {
    stop("The `kdml_mcmc` object has malformed kernel draws.",
         call. = FALSE)
  }

  occupancy <- list()
  transitions <- list()
  occupancy_index <- 0L
  transition_index <- 0L
  for (feature_index in seq_len(dimensions$features)) {
    feature <- dimensions$feature_names[feature_index]
    candidates <- NULL
    if (!is.null(object$candidates) && !is.null(object$candidates[[feature]])) {
      candidates <- as.character(object$candidates[[feature]])
    }
    observed <- as.character(object$kernel[, feature_index, ])
    kernel_levels <- unique(c(candidates, observed[!is.na(observed)]))

    for (chain_index in seq_len(dimensions$chains)) {
      chain <- dimensions$chain_names[chain_index]
      values <- as.character(object$kernel[, feature_index, chain_index])
      valid <- !is.na(values)
      counts <- table(factor(values[valid], levels = kernel_levels))
      denominator <- sum(counts)
      proportions <- if (denominator > 0L) {
        as.numeric(counts) / denominator
      } else {
        rep(NA_real_, length(counts))
      }

      if (length(kernel_levels)) {
        occupancy_index <- occupancy_index + 1L
        occupancy[[occupancy_index]] <- data.frame(
          feature = feature,
          chain = chain,
          kernel = kernel_levels,
          count = as.integer(counts),
          proportion = proportions,
          stringsAsFactors = FALSE
        )
      }

      if (length(values) > 1L) {
        complete_pair <- !is.na(values[-1L]) &
          !is.na(values[-length(values)])
        transition_count <- sum(
          values[-1L][complete_pair] !=
            values[-length(values)][complete_pair]
        )
        opportunities <- sum(complete_pair)
      } else {
        transition_count <- 0L
        opportunities <- 0L
      }
      transition_index <- transition_index + 1L
      transitions[[transition_index]] <- data.frame(
        feature = feature,
        chain = chain,
        transitions = as.integer(transition_count),
        opportunities = as.integer(opportunities),
        transition_rate = if (opportunities > 0L) {
          transition_count / opportunities
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  }

  occupancy_frame <- if (length(occupancy)) {
    do.call(rbind, occupancy)
  } else {
    data.frame(
      feature = character(), chain = character(), kernel = character(),
      count = integer(), proportion = numeric(),
      stringsAsFactors = FALSE
    )
  }
  transition_frame <- do.call(rbind, transitions)
  rownames(occupancy_frame) <- NULL
  rownames(transition_frame) <- NULL
  list(occupancy = occupancy_frame, transitions = transition_frame)
}

.kdml_diagnostic_matrix <- function(field, object, dimensions) {
  value <- object$diagnostics[[field]]
  expected <- c(dimensions$features, dimensions$chains)
  if (is.null(value) ||
      !identical(as.integer(dim(value)), as.integer(expected))) {
    stop("The `kdml_mcmc` object has malformed `", field,
         "` diagnostics.", call. = FALSE)
  }
  as.numeric(value)
}

.kdml_acceptance_diagnostics <- function(object, dimensions) {
  if (is.null(object$diagnostics) || !is.list(object$diagnostics)) {
    stop("The `kdml_mcmc` object has malformed move diagnostics.",
         call. = FALSE)
  }
  fields <- c(
    "proposal_sd",
    "kernel_attempts", "kernel_accepts",
    "bandwidth_attempts", "bandwidth_accepts",
    "sampling_kernel_attempts", "sampling_kernel_accepts",
    "sampling_bandwidth_attempts", "sampling_bandwidth_accepts"
  )
  values <- lapply(
    fields, .kdml_diagnostic_matrix,
    object = object, dimensions = dimensions
  )
  names(values) <- fields

  answer <- data.frame(
    feature = rep(dimensions$feature_names, dimensions$chains),
    chain = rep(dimensions$chain_names, each = dimensions$features),
    stringsAsFactors = FALSE
  )
  for (field in fields) {
    answer[[field]] <- values[[field]]
  }
  rate <- function(accepts, attempts) {
    ifelse(attempts > 0, accepts / attempts, NA_real_)
  }
  answer$kernel_acceptance <- rate(
    answer$kernel_accepts, answer$kernel_attempts
  )
  answer$bandwidth_acceptance <- rate(
    answer$bandwidth_accepts, answer$bandwidth_attempts
  )
  answer$sampling_kernel_acceptance <- rate(
    answer$sampling_kernel_accepts, answer$sampling_kernel_attempts
  )
  answer$sampling_bandwidth_acceptance <- rate(
    answer$sampling_bandwidth_accepts,
    answer$sampling_bandwidth_attempts
  )

  invalid <- object$invalid_proposals
  if (is.null(invalid) || length(invalid) != dimensions$chains) {
    invalid <- rep(NA_integer_, dimensions$chains)
  }
  answer$invalid_proposals <- rep(
    as.integer(invalid), each = dimensions$features
  )
  answer
}

.kdml_convergence_warnings <- function(diagnostics, dimensions,
                                       rhat_threshold,
                                       min_ess_per_chain) {
  warnings <- character()
  if (dimensions$chains < 2L) {
    warnings <- c(
      warnings,
      "Only one chain was run; cross-chain convergence cannot be assessed reliably."
    )
  }

  unavailable_rhat <- diagnostics$feature[!is.finite(diagnostics$rhat)]
  if (length(unavailable_rhat)) {
    warnings <- c(
      warnings,
      paste0(
        "R-hat is unavailable for: ",
        paste(unavailable_rhat, collapse = ", "), "."
      )
    )
  }
  high_rhat <- diagnostics$feature[
    is.finite(diagnostics$rhat) &
      diagnostics$rhat > rhat_threshold
  ]
  if (length(high_rhat)) {
    warnings <- c(
      warnings,
      paste0(
        "R-hat exceeds ", format(rhat_threshold), " for: ",
        paste(high_rhat, collapse = ", "), "."
      )
    )
  }

  ess_target <- min_ess_per_chain * dimensions$chains
  unavailable_ess <- diagnostics$feature[
    !is.finite(diagnostics$ess_bulk) |
      !is.finite(diagnostics$ess_tail)
  ]
  if (length(unavailable_ess)) {
    warnings <- c(
      warnings,
      paste0(
        "Bulk or tail ESS is unavailable for: ",
        paste(unavailable_ess, collapse = ", "), "."
      )
    )
  }
  low_ess <- diagnostics$feature[
    (is.finite(diagnostics$ess_bulk) &
       diagnostics$ess_bulk < ess_target) |
      (is.finite(diagnostics$ess_tail) &
       diagnostics$ess_tail < ess_target)
  ]
  if (length(low_ess)) {
    warnings <- c(
      warnings,
      paste0(
        "Bulk or tail ESS is below ", format(ess_target),
        " (", format(min_ess_per_chain), " per chain) for: ",
        paste(low_ess, collapse = ", "), "."
      )
    )
  }
  unique(warnings)
}

#'
#' @export
as_draws.kdml_mcmc <- function(
    x,
    variables = c("theta", "bandwidth", "score", "log_target"),
    ...) {
  posterior::as_draws_array(x, variables = variables, ...)
}

#'
#' @export
summary.kdml_mcmc <- function(object, rhat_threshold = 1.01,
                              min_ess_per_chain = 100, ...) {
  dimensions <- .kdml_mcmc_dimensions(object)
  if (length(rhat_threshold) != 1L || !is.numeric(rhat_threshold) ||
      !is.finite(rhat_threshold) || rhat_threshold < 1) {
    stop("`rhat_threshold` must be one finite number greater than or equal to 1.",
         call. = FALSE)
  }
  if (length(min_ess_per_chain) != 1L ||
      !is.numeric(min_ess_per_chain) ||
      !is.finite(min_ess_per_chain) || min_ess_per_chain <= 0) {
    stop("`min_ess_per_chain` must be one finite, strictly positive number.",
         call. = FALSE)
  }

  theta_diagnostics <- suppressWarnings(kdml.diagnostics(
    object, variables = "theta"
  ))
  theta_diagnostics$feature <- dimensions$feature_names
  theta_diagnostics <- theta_diagnostics[
    c("feature", setdiff(names(theta_diagnostics), "feature"))
  ]
  kernel <- .kdml_kernel_diagnostics(object, dimensions)
  acceptance <- .kdml_acceptance_diagnostics(object, dimensions)
  convergence_warnings <- .kdml_convergence_warnings(
    theta_diagnostics, dimensions, rhat_threshold, min_ess_per_chain
  )

  control <- object$control
  metadata <- list(
    distance = if (!is.null(object$distance)) object$distance else NA_character_,
    backend = if (!is.null(object$backend)) object$backend else NA_character_,
    features = dimensions$features,
    chains = dimensions$chains,
    chain_threads = if (!is.null(control$chain_threads)) {
      control$chain_threads
    } else {
      1L
    },
    warmup = if (!is.null(control$warmup)) control$warmup else NA_integer_,
    retained_draws_per_chain = dimensions$draws,
    thin = if (!is.null(control$thin)) control$thin else NA_integer_,
    beta = if (!is.null(control$beta)) control$beta else NA_real_,
    adapt = if (!is.null(control$adapt)) control$adapt else NA
  )

  structure(
    list(
      call = object$call,
      metadata = metadata,
      score = .kdml_score_summary(object),
      map = object$map,
      theta = theta_diagnostics,
      convergence = list(
        ok = !length(convergence_warnings),
        warnings = convergence_warnings,
        thresholds = list(
          rhat = rhat_threshold,
          min_ess_per_chain = min_ess_per_chain
        )
      ),
      kernel_occupancy = kernel$occupancy,
      kernel_transitions = kernel$transitions,
      acceptance = acceptance
    ),
    class = "summary.kdml_mcmc"
  )
}

.kdml_pooled_kernel_occupancy <- function(occupancy) {
  if (!nrow(occupancy)) {
    return(occupancy)
  }
  key <- interaction(
    occupancy$feature, occupancy$kernel,
    drop = TRUE, lex.order = TRUE
  )
  count <- rowsum(occupancy$count, key, reorder = FALSE)
  first <- !duplicated(key)
  answer <- occupancy[first, c("feature", "kernel"), drop = FALSE]
  answer$count <- as.integer(count[, 1L])
  totals <- ave(answer$count, answer$feature, FUN = sum)
  answer$proportion <- ifelse(totals > 0, answer$count / totals, NA_real_)
  rownames(answer) <- NULL
  answer
}

.kdml_pooled_acceptance <- function(acceptance) {
  features <- unique(acceptance$feature)
  answer <- lapply(features, function(feature) {
    rows <- acceptance$feature == feature
    kernel_attempts <- sum(
      acceptance$sampling_kernel_attempts[rows], na.rm = TRUE
    )
    kernel_accepts <- sum(
      acceptance$sampling_kernel_accepts[rows], na.rm = TRUE
    )
    bandwidth_attempts <- sum(
      acceptance$sampling_bandwidth_attempts[rows], na.rm = TRUE
    )
    bandwidth_accepts <- sum(
      acceptance$sampling_bandwidth_accepts[rows], na.rm = TRUE
    )
    data.frame(
      feature = feature,
      kernel = if (kernel_attempts > 0) {
        kernel_accepts / kernel_attempts
      } else {
        NA_real_
      },
      bandwidth = if (bandwidth_attempts > 0) {
        bandwidth_accepts / bandwidth_attempts
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, answer)
}

#' @rdname summary.kdml_mcmc
#' @export
print.summary.kdml_mcmc <- function(x, digits = max(3L, getOption("digits") - 3L),
                                    ...) {
  cat("KDML MCMC summary\n")
  cat("  distance: ", toupper(x$metadata$distance), "\n", sep = "")
  cat("  backend: ", x$metadata$backend, "\n", sep = "")
  cat("  features: ", x$metadata$features, "\n", sep = "")
  cat("  chains: ", x$metadata$chains,
      ", chain threads: ", x$metadata$chain_threads,
      ", retained draws per chain: ",
      x$metadata$retained_draws_per_chain, "\n", sep = "")
  cat("  warmup: ", x$metadata$warmup,
      ", thin: ", x$metadata$thin,
      ", beta: ", format(x$metadata$beta, digits = digits),
      ", adaptation: ", x$metadata$adapt, "\n", sep = "")

  cat("\nRetained score summary:\n")
  print(x$score, row.names = FALSE, digits = digits)
  if (!is.null(x$map)) {
    cat("\nRetained MAP estimate: draw ", x$map$draw,
        ", chain ", x$map$chain,
        ", log target ", format(x$map$log_target, digits = digits),
        ", MSCV score ", format(x$map$score, digits = digits), "\n",
        sep = "")
    if (!is.null(x$map$kernel) && !is.null(x$map$bandwidth)) {
      map_features <- names(x$map$kernel)
      if (is.null(map_features)) {
        map_features <- if (length(x$map$kernel) == nrow(x$theta)) {
          x$theta$feature
        } else {
          paste0("feature", seq_along(x$map$kernel))
        }
      }
      map <- data.frame(
        feature = map_features,
        kernel = as.character(x$map$kernel),
        bandwidth = as.numeric(x$map$bandwidth),
        stringsAsFactors = FALSE
      )
      print(map, row.names = FALSE, digits = digits)
    }
  }

  cat("\nTheta diagnostics:\n")
  print(
    x$theta[
      c("feature", "mean", "sd", "rhat", "ess_bulk", "ess_tail",
        "mcse_mean")
    ],
    row.names = FALSE, digits = digits
  )
  if (x$convergence$ok) {
    cat("\nNo convergence warnings at the configured thresholds.\n")
  } else {
    cat("\nConvergence warnings:\n")
    for (warning in x$convergence$warnings) {
      cat("  - ", warning, "\n", sep = "")
    }
  }

  occupancy <- .kdml_pooled_kernel_occupancy(x$kernel_occupancy)
  cat("\nPooled kernel occupancy:\n")
  print(occupancy, row.names = FALSE, digits = digits)

  transitions <- stats::aggregate(
    x$kernel_transitions[c("transitions", "opportunities")],
    by = list(feature = x$kernel_transitions$feature),
    FUN = sum
  )
  transitions$transition_rate <- ifelse(
    transitions$opportunities > 0,
    transitions$transitions / transitions$opportunities,
    NA_real_
  )
  cat("\nKernel transitions:\n")
  print(transitions, row.names = FALSE, digits = digits)

  cat("\nPost-warmup acceptance:\n")
  print(
    .kdml_pooled_acceptance(x$acceptance),
    row.names = FALSE, digits = digits
  )
  invisible(x)
}

.kdml_select_dimension <- function(value, available, argument) {
  if (is.null(value)) {
    return(seq_along(available))
  }
  if (!length(value)) {
    stop("`", argument, "` must select at least one value.",
         call. = FALSE)
  }
  if (is.character(value)) {
    missing <- setdiff(value, available)
    if (length(missing)) {
      stop("Unknown ", argument, ": ", paste(missing, collapse = ", "),
           ".", call. = FALSE)
    }
    return(match(unique(value), available))
  }
  if (is.numeric(value) && length(value) &&
      all(is.finite(value)) && all(value == floor(value)) &&
      all(value >= 1L) && all(value <= length(available))) {
    return(unique(as.integer(value)))
  }
  stop("`", argument, "` must select valid names or indices.",
       call. = FALSE)
}

.kdml_plot_layout <- function(panel_count) {
  columns <- ceiling(sqrt(panel_count))
  rows <- ceiling(panel_count / columns)
  graphics::par(mfrow = c(rows, columns))
}

.kdml_do_plot <- function(fun, arguments, dots) {
  if (length(dots)) {
    if (is.null(names(dots)) || any(!nzchar(names(dots)))) {
      stop("Additional graphical parameters in `...` must be named.",
           call. = FALSE)
    }
    arguments[names(dots)] <- dots
  }
  do.call(fun, arguments)
}

.kdml_bandwidth_density_panel <- function(value, feature_name,
                                          log_bandwidth, dots) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  title <- paste(feature_name, "bandwidth posterior")

  if (!length(value)) {
    graphics::plot.new()
    graphics::title(main = title)
    graphics::text(0.5, 0.5, "No finite draws")
    return(invisible(NULL))
  }

  if (isTRUE(log_bandwidth)) {
    value <- value[value > 0]
    if (!length(value)) {
      graphics::plot.new()
      graphics::title(main = title)
      graphics::text(0.5, 0.5, "No positive finite draws")
      return(invisible(NULL))
    }
    value <- log(value)
    x_label <- "Log bandwidth"
  } else {
    value <- value[value >= 0 & value <= 1]
    if (!length(value)) {
      graphics::plot.new()
      graphics::title(main = title)
      graphics::text(0.5, 0.5, "No draws in [0, 1]")
      return(invisible(NULL))
    }
    x_label <- "Bandwidth"
  }

  spread <- diff(range(value))
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(value)))
  if (length(value) >= 2L && is.finite(spread) && spread > tolerance) {
    estimate <- if (isTRUE(log_bandwidth)) {
      stats::density(value)
    } else {
      stats::density(value, from = 0, to = 1, cut = 0)
    }
    arguments <- list(
      x = estimate, type = "l", lwd = 2, col = "#3366AA",
      main = title, xlab = x_label, ylab = "Posterior density"
    )
    if (!isTRUE(log_bandwidth)) {
      arguments$xlim <- c(0, 1)
    }
    .kdml_do_plot(
      graphics::plot,
      arguments,
      dots
    )
  } else {
    centre <- value[1L]
    limits <- if (isTRUE(log_bandwidth)) {
      padding <- max(0.5, abs(centre) * 0.05)
      c(centre - padding, centre + padding)
    } else {
      c(0, 1)
    }
    .kdml_do_plot(
      graphics::plot,
      list(
        x = limits, y = c(0, 1),
        type = "n", main = title, xlab = x_label,
        ylab = "Posterior density"
      ),
      dots
    )
    graphics::abline(v = centre, lwd = 2, col = "#3366AA")
    graphics::text(centre, 0.5, "Constant draws", pos = 4, cex = 0.8)
  }
  invisible(NULL)
}

.kdml_trace_panel <- function(values, title, ylab, chain_names, dots) {
  values <- as.matrix(values)
  draw_count <- nrow(values)
  colors <- seq_len(ncol(values))
  arguments <- list(
    x = seq_len(draw_count),
    y = values,
    type = if (draw_count > 1L) "l" else "p",
    lty = 1,
    col = colors,
    xlab = "Retained draw",
    ylab = ylab,
    main = title
  )
  .kdml_do_plot(graphics::matplot, arguments, dots)
  if (length(chain_names) > 1L) {
    graphics::legend(
      "topright", legend = chain_names, col = colors, lty = 1,
      bty = "n", cex = 0.8
    )
  }
}

.kdml_safe_acf <- function(value, lag_max) {
  value <- as.numeric(value)
  answer <- rep(NA_real_, lag_max + 1L)
  answer[1L] <- 1
  if (length(value) > 1L && is.finite(stats::sd(value)) &&
      stats::sd(value) > 0) {
    computed <- suppressWarnings(stats::acf(
      value, plot = FALSE, lag.max = lag_max, na.action = stats::na.pass
    ))
    result <- as.numeric(computed$acf)
    answer[seq_along(result)] <- result
  }
  answer
}

.kdml_selected_theta_diagnostics <- function(object, dimensions,
                                             feature_index, chain_index) {
  variables <- paste0(
    "theta[", dimensions$feature_names[feature_index], "]"
  )
  draws <- posterior::as_draws_array(object, variables = "theta")
  draws <- posterior::subset_draws(
    draws,
    variable = variables,
    chain = chain_index
  )
  answer <- suppressWarnings(posterior::summarise_draws(
    draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"
  ))
  as.data.frame(answer, stringsAsFactors = FALSE)
}

.kdml_diagnostic_plot <- function(object, dimensions, feature_index,
                                  chain_index, dots) {
  diagnostic <- .kdml_selected_theta_diagnostics(
    object, dimensions, feature_index, chain_index
  )
  labels <- dimensions$feature_names[feature_index]
  graphics::par(mfrow = c(1L, 2L))

  rhat <- diagnostic$rhat
  rhat_height <- ifelse(is.finite(rhat), rhat, 0)
  rhat_limit <- max(c(1.05, rhat[is.finite(rhat)]), na.rm = TRUE)
  positions <- .kdml_do_plot(
    graphics::barplot,
    list(
      height = rhat_height, names.arg = labels, las = 2,
      ylim = c(0, rhat_limit * 1.05), ylab = "R-hat",
      main = "Rank-normalized R-hat", col = "grey70"
    ),
    dots
  )
  graphics::abline(h = 1.01, col = "red", lty = 2)
  if (any(!is.finite(rhat))) {
    graphics::text(
      positions[!is.finite(rhat)], rep(0, sum(!is.finite(rhat))),
      labels = "NA", pos = 3, cex = 0.8
    )
  }

  ess <- rbind(bulk = diagnostic$ess_bulk, tail = diagnostic$ess_tail)
  ess_height <- ess
  ess_height[!is.finite(ess_height)] <- 0
  ess_target <- 100 * length(chain_index)
  ess_limit <- max(c(ess_target, ess[is.finite(ess)]), na.rm = TRUE)
  .kdml_do_plot(
    graphics::barplot,
    list(
      height = ess_height, beside = TRUE, names.arg = labels, las = 2,
      ylim = c(0, max(1, ess_limit * 1.05)), ylab = "ESS",
      main = "Bulk and tail ESS", col = c("#4477AA", "#CC6677"),
      legend.text = rownames(ess)
    ),
    dots
  )
  graphics::abline(h = ess_target, col = "red", lty = 2)
}

#'
#' @export
plot.kdml_mcmc <- function(
    x,
    type = c("trace", "acf", "diagnostics", "score", "bandwidth",
             "kernel", "acceptance"),
    features = NULL,
    chains = NULL,
    lag.max = NULL,
    ...) {
  type <- match.arg(type)
  dimensions <- .kdml_mcmc_dimensions(x)
  feature_index <- .kdml_select_dimension(
    features, dimensions$feature_names, "features"
  )
  chain_index <- .kdml_select_dimension(
    chains, dimensions$chain_names, "chains"
  )
  dots <- list(...)
  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_parameters), add = TRUE)

  if (type == "trace") {
    .kdml_plot_layout(length(feature_index))
    for (index in feature_index) {
      values <- matrix(
        x$theta[, index, chain_index, drop = FALSE],
        nrow = dimensions$draws,
        ncol = length(chain_index)
      )
      .kdml_trace_panel(
        values, paste(dimensions$feature_names[index], "theta trace"),
        expression(theta), dimensions$chain_names[chain_index], dots
      )
    }
  } else if (type == "acf") {
    if (is.null(lag.max)) {
      lag.max <- min(20L, dimensions$draws - 1L)
    }
    if (length(lag.max) != 1L || !is.numeric(lag.max) ||
        !is.finite(lag.max) || lag.max < 0 || lag.max != floor(lag.max)) {
      stop("`lag.max` must be one non-negative integer.", call. = FALSE)
    }
    lag.max <- min(as.integer(lag.max), dimensions$draws - 1L)
    .kdml_plot_layout(length(feature_index))
    for (index in feature_index) {
      acf_values <- vapply(
        chain_index,
        function(chain) {
          .kdml_safe_acf(x$theta[, index, chain], lag.max)
        },
        numeric(lag.max + 1L)
      )
      acf_values <- matrix(
        acf_values,
        nrow = lag.max + 1L,
        ncol = length(chain_index)
      )
      colors <- seq_along(chain_index)
      offsets <- if (length(chain_index) > 1L) {
        seq(-0.2, 0.2, length.out = length(chain_index))
      } else {
        0
      }
      lag_values <- 0:lag.max
      acf_lags <- outer(lag_values, offsets, `+`)
      arguments <- list(
        x = acf_lags,
        y = acf_values,
        type = "h",
        lty = 1,
        lwd = 2,
        col = colors,
        xlim = c(-0.5, max(0.5, lag.max + 0.5)),
        ylim = c(-1, 1),
        xlab = "Lag",
        ylab = "Autocorrelation",
        main = paste(dimensions$feature_names[index], "theta ACF")
      )
      .kdml_do_plot(graphics::matplot, arguments, dots)
      graphics::abline(h = 0, col = "grey70")
      if (length(chain_index) > 1L) {
        graphics::legend(
          "topright", legend = dimensions$chain_names[chain_index],
          col = colors, lty = 1, lwd = 2, bty = "n", cex = 0.8
        )
      }
    }
  } else if (type == "diagnostics") {
    .kdml_diagnostic_plot(
      x, dimensions, feature_index, chain_index, dots
    )
  } else if (type == "score") {
    fields <- c("score", "log_target")
    labels <- c("MSCV score", "Log target")
    graphics::par(mfrow = c(2L, 1L))
    for (index in seq_along(fields)) {
      values <- x[[fields[index]]][, chain_index, drop = FALSE]
      .kdml_trace_panel(
        values, paste(labels[index], "trace"), labels[index],
        dimensions$chain_names[chain_index], dots
      )
    }
  } else if (type == "bandwidth") {
    .kdml_plot_layout(length(feature_index))
    for (index in feature_index) {
      value <- as.numeric(x$bandwidth[, index, chain_index, drop = FALSE])
      feature_type <- x$preprocessing$type
      log_bandwidth <- !is.null(feature_type) &&
        length(feature_type) >= index &&
        isTRUE(as.integer(feature_type[index]) == 0L)
      .kdml_bandwidth_density_panel(
        value = value,
        feature_name = dimensions$feature_names[index],
        log_bandwidth = log_bandwidth,
        dots = dots
      )
    }
  } else if (type == "kernel") {
    .kdml_plot_layout(length(feature_index))
    for (index in feature_index) {
      values <- x$kernel[, index, chain_index, drop = FALSE]
      levels <- unique(as.character(values))
      levels <- levels[!is.na(levels)]
      counts <- vapply(
        seq_along(chain_index),
        function(position) {
          table(factor(values[, 1L, position], levels = levels))
        },
        numeric(length(levels))
      )
      if (!is.matrix(counts)) {
        counts <- matrix(counts, ncol = length(chain_index))
      }
      proportions <- apply(counts, 2L, function(value) {
        if (sum(value) > 0) value / sum(value) else value
      })
      if (!is.matrix(proportions)) {
        proportions <- matrix(proportions, ncol = length(chain_index))
      }
      rownames(proportions) <- levels
      .kdml_do_plot(
        graphics::barplot,
        list(
          height = proportions, beside = TRUE,
          names.arg = dimensions$chain_names[chain_index],
          ylim = c(0, 1), ylab = "Posterior probability",
          main = paste(
            dimensions$feature_names[index], "kernel posterior"
          ),
          legend.text = levels
        ),
        dots
      )
    }
  } else {
    acceptance <- .kdml_acceptance_diagnostics(x, dimensions)
    acceptance <- acceptance[
      acceptance$feature %in% dimensions$feature_names[feature_index] &
        acceptance$chain %in% dimensions$chain_names[chain_index],
      , drop = FALSE
    ]
    pooled <- .kdml_pooled_acceptance(acceptance)
    rates <- rbind(
      kernel = pooled$kernel,
      bandwidth = pooled$bandwidth
    )
    plotted <- rates
    plotted[!is.finite(plotted)] <- 0
    graphics::par(mfrow = c(1L, 1L))
    .kdml_do_plot(
      graphics::barplot,
      list(
        height = plotted, beside = TRUE, names.arg = pooled$feature,
        ylim = c(0, 1), las = 2, ylab = "Acceptance rate",
        main = "Post-warmup acceptance",
        col = c("#4477AA", "#CC6677"), legend.text = rownames(rates)
      ),
      dots
    )
  }
  invisible(x)
}
