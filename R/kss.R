kss <- function(df, bw = "np", npmethod = NULL, cFUN = "c_gaussian",
                uFUN = "u_aitken", oFUN = "o_wangvanryzin",
                nstart = NULL, stan = TRUE, verbose = FALSE) {
  stan <- .kdml_validate_flag(stan, "stan")
  verbose <- .kdml_validate_flag(verbose, "verbose")
  prep <- .kdml_prepare_data(df, drop_unused = FALSE)

  if (is.null(nstart)) {
    nstart <- min(3L, ncol(df))
  }
  nstart <- .kdml_scalar_integer(nstart, "nstart", 1L)

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
    stop("No bandwidth selection chosen. Input a numeric bandwidth vector or ",
         "choose 'np'.", call. = FALSE)
  }
  if (is.numeric(bw)) {
    bws <- .kdml_expand_numeric(bw, prep, "bw", positive = TRUE)
  } else {
    if (!is.character(bw) || length(bw) != 1L || is.na(bw) || bw != "np") {
      stop("Invalid bandwidth selection. Input a numeric vector or choose 'np'.",
           call. = FALSE)
    }
    if (is.null(npmethod)) {
      if (verbose) {
        warning(
          "No npmethod selected, defaulting to maximum likelihood ",
          "cross-validation (cv.ml). You may also use least-squares ",
          "cross-validation (cv.ls) or normal reference (normal-reference).",
          call. = FALSE
        )
      }
      npmethod <- "cv.ml"
    } else if (!is.character(npmethod) || length(npmethod) != 1L ||
               is.na(npmethod) ||
               !npmethod %in% c("cv.ml", "cv.ls", "normal-reference")) {
      stop(
        "`npmethod` must be 'cv.ml', 'cv.ls', or 'normal-reference'.",
        call. = FALSE
      )
    }

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

    selected <- npudensbw(
      df_ordered, ckertype = cker, ukertype = uker, okertype = oker,
      bwmethod = npmethod, nmulti = nstart
    )
    bws <- setNames(as.numeric(selected$bw), type_order)
    bws <- .kdml_expand_numeric(bws, prep, "bw", positive = TRUE)
    if (verbose) {
      message("Bandwidth calculation complete. Computing similarities.")
    }
  }

  selected_kernels <- setNames(
    c(cFUN, uFUN, oFUN)[prep$type + 1L], prep$feature_names
  )
  evaluated <- .kdml_mcmc_evaluate(
    df, distance = "dkss", kernels = selected_kernels,
    bandwidths = bws, return_similarity = TRUE, backend = "cpu",
    drop_unused = FALSE
  )
  similarities <- evaluated$similarity
  dimnames(similarities) <- NULL
  bandwidths <- matrix(
    bws[type_order], nrow = 1L,
    dimnames = list(NULL, type_order)
  )

  if (!stan) {
    if (verbose) {
      message("Completed similarity calculation.")
    }
    return(list(similarities = similarities, bandwidths = bandwidths))
  }

  minimum <- min(similarities)
  maximum <- max(similarities)
  standardized <- (similarities - minimum) / (maximum - minimum)
  if (verbose) {
    message("Completed similarity calculation.")
  }
  list(similarities = standardized, bandwidths = bandwidths)
}
