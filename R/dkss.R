dkss <- function(df, bw = "mscv", cFUN = "c_gaussian", uFUN = "u_aitken",
                 oFUN = "o_wangvanryzin", stan = TRUE, verbose = FALSE) {
  stan <- .kdml_validate_flag(stan, "stan")
  verbose <- .kdml_validate_flag(verbose, "verbose")

  if (inherits(df, "kdml_mcmc")) {
    if (!identical(df$distance, "dkss")) {
      stop("`dkss()` requires a fit produced by `mcmc.dkss()`.",
           call. = FALSE)
    }
    if (!missing(bw) || !missing(cFUN) || !missing(uFUN) || !missing(oFUN)) {
      stop("Do not supply `bw`, `cFUN`, `uFUN`, or `oFUN` ",
           "when `df` is an MCMC fit.", call. = FALSE)
    }
    if (verbose) {
      message("Computing DKSS distances from the fitted MCMC state.")
    }
    return(.kdml_map_distance(df, standardize = stan))
  }

  prep <- .kdml_prepare_data(df, drop_unused = FALSE)
  v_ck <- c(
    "c_gaussian", "c_epanechnikov", "c_uniform", "c_triangle",
    "c_biweight", "c_triweight", "c_tricube", "c_cosine",
    "c_logistic", "c_sigmoid", "c_silverman"
  )
  v_uk <- c("u_aitken", "u_aitchisonaitken")
  v_ok <- c(
    "o_wangvanryzin", "o_habbema", "o_aitken",
    "o_aitchisonaitken", "o_liracine"
  )
  cFUN <- .kdml_validate_scalar_kernel(cFUN, v_ck, "cFUN")
  uFUN <- .kdml_validate_scalar_kernel(uFUN, v_uk, "uFUN")
  oFUN <- .kdml_validate_scalar_kernel(oFUN, v_ok, "oFUN")

  type_order <- c(
    prep$feature_names[prep$type == 0L],
    prep$feature_names[prep$type == 1L],
    prep$feature_names[prep$type == 2L]
  )
  df_ordered <- df[, type_order, drop = FALSE]

  if (is.null(bw)) {
    stop("No bandwidth selection chosen. Input a numeric vector of bandwidths, ",
         "or choose from 'np' or 'mscv' bandwidth selection methods.",
         call. = FALSE)
  }
  if (is.numeric(bw)) {
    bws <- .kdml_expand_numeric(bw, prep, "bw", positive = TRUE)
  } else {
    if (!is.character(bw) || length(bw) != 1L || is.na(bw)) {
      stop("`bw` must be a numeric vector, 'np', or 'mscv'.", call. = FALSE)
    }
    if (bw == "mscv") {
      selected <- mscv.dkss(
        df, nstart = min(3L, ncol(df)), ckernel = cFUN,
        ukernel = uFUN, okernel = oFUN, verbose = verbose
      )
      bws <- setNames(as.numeric(selected$bw[, 1L]), rownames(selected$bw))
      bws <- .kdml_expand_numeric(bws, prep, "bw", positive = TRUE)
      if (verbose) {
        message("Bandwidth calculation complete. Computing distances.")
      }
    } else if (bw == "np") {
      if (cFUN == "c_gaussian") cker <- "gaussian"
      if (cFUN == "c_epanechnikov") cker <- "epanechnikov"
      if (cFUN == "c_uniform") cker <- "uniform"
      if (cFUN %in% setdiff(v_ck, c(
        "c_gaussian", "c_epanechnikov", "c_uniform"
      ))) {
        stop("Choose one of c_gaussian, c_epanechnikov, or c_uniform for ",
             "continuous kernels while using np.", call. = FALSE)
      }
      uker <- if (uFUN == "u_aitken") "liracine" else "aitchisonaitken"
      if (oFUN == "o_wangvanryzin") oker <- "wangvanryzin"
      if (oFUN == "o_liracine") oker <- "liracine"
      if (oFUN %in% c("o_habbema", "o_aitken", "o_aitchisonaitken")) {
        stop("Choose one of o_wangvanryzin or o_liracine for ordinal kernels ",
             "while using np.", call. = FALSE)
      }
      if (verbose) {
        selected <- npudensbw(
          df_ordered, ckertype = cker, ukertype = uker, okertype = oker,
          bwmethod = "cv.ml"
        )
      } else {
        invisible(capture.output(
          selected <- npudensbw(
            df_ordered, ckertype = cker, ukertype = uker, okertype = oker,
            bwmethod = "cv.ml"
          )
        ))
      }
      bws <- setNames(as.numeric(selected$bw), type_order)
      bws <- .kdml_expand_numeric(bws, prep, "bw", positive = TRUE)
      if (verbose) {
        message("Bandwidth calculation complete. Computing distances.")
      }
    } else {
      stop("Invalid bandwidth selection. Input a numeric vector, 'np', or ",
           "'mscv'.", call. = FALSE)
    }
  }

  selected_kernels <- setNames(
    c(cFUN, uFUN, oFUN)[prep$type + 1L], prep$feature_names
  )
  answer <- .kdml_fixed_distance(
    df, "dkss", bws, selected_kernels, standardize = stan,
    drop_unused = FALSE
  )
  dimnames(answer$distances) <- NULL
  answer$bandwidths <- answer$bandwidths[, type_order, drop = FALSE]
  if (verbose) {
    message("Completed distance calculation.")
  }
  list(distances = answer$distances, bandwidths = answer$bandwidths)
}
