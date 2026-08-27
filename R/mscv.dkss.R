mscv.dkss <- function(df, nstart = NULL, ckernel = "c_gaussian",
                      ukernel = "u_aitken", okernel = "o_wangvanryzin", verbose = FALSE) {
  prep <- .kdml_prepare_data(df)
  verbose <- .kdml_validate_flag(verbose, "verbose")
  if (is.null(nstart)) {
    nstart <- min(3L, ncol(df))
    message("No nstart value given, defaulting to ", nstart)
  }
  nstart <- .kdml_scalar_integer(nstart, "nstart", 1L)
  verb <- verbose
  if (identical(ukernel, "u.aitchisonaitken")) {
    ukernel <- "u_aitchisonaitken"
  }
  v_ck <- c("c_gaussian", "c_epanechnikov", "c_uniform", "c_triangle",
            "c_biweight", "c_triweight", "c_tricube", "c_cosine", 
            "c_logistic", "c_sigmoid", "c_silverman")
  v_uk <- c("u_aitken", "u_aitchisonaitken")
  v_ok <- c("o_wangvanryzin", "o_habbema", "o_aitken", "o_aitchisonaitken", "o_liracine")
  
  ckernel <- .kdml_validate_scalar_kernel(ckernel, v_ck, "ckernel")
  ukernel <- .kdml_validate_scalar_kernel(ukernel, v_uk, "ukernel")
  okernel <- .kdml_validate_scalar_kernel(okernel, v_ok, "okernel")
  
  # Get column names by type
  con_cols <- prep$feature_names[prep$type == 0L]
  fac_cols <- prep$feature_names[prep$type == 1L]
  ord_cols <- prep$feature_names[prep$type == 2L]
  
  # Reorder the dataframe
  df_ordered <- prep$x[, c(con_cols, fac_cols, ord_cols), drop = FALSE]
  
  # store indices
  con_ind <- length(con_cols) #index of continuous variables
  fac_ind <- con_ind + length(fac_cols) #index of factors
  ord_ind <- fac_ind + length(ord_cols) #index of ordinal variables
  
  n <- nrow(df_ordered)
  N <- ncol(df_ordered)
  nominal_categories <- if (fac_ind > con_ind) {
    vapply(
      (con_ind + 1):fac_ind,
      function(i) length(unique(df_ordered[, i])),
      integer(1)
    )
  } else {
    integer()
  }
  if (ukernel == "u_aitchisonaitken" &&
      any(nominal_categories < 2L)) {
    stop("Aitchison-Aitken kernels require at least two observed categories.")
  }
  nominal_upper <- (nominal_categories - 1) / nominal_categories
  parameter_upper <- rep(1, N)
  if (ukernel == "u_aitchisonaitken" && fac_ind > con_ind) {
    parameter_upper[(con_ind + 1):fac_ind] <- nominal_upper
  }
  
  # MSCV optimization function
  mscv_opt <- function(lambda) {
    penalty <- sum(c(
      if (con_ind > 0) lambda[1:con_ind] <= 0 else numeric(),
      if (fac_ind > con_ind) sapply((con_ind + 1):fac_ind, function(i) {
        max_val <- if (ukernel == "u_aitchisonaitken") {
          nominal_upper[i - con_ind]
        } else 1
        lambda[i] < 0 || lambda[i] > max_val
      }) else numeric(),
      if (ord_ind > fac_ind) lambda[(fac_ind + 1):ord_ind] < 0 | lambda[(fac_ind + 1):ord_ind] > 1 else numeric()
    ))
    if (penalty > 0) return(-Inf)
    
    outK <- outL <- outell <- K <-  L <-  ell <- list()
    
    if(con_ind != 0){
      for (i in 1:con_ind) { outK[[i]] <- outer(df_ordered[, i], df_ordered[, i], "-") }
    }
    if(fac_ind > con_ind){
      for (i in (con_ind + 1):fac_ind) { outL[[i - con_ind]] <- outer(df_ordered[, i], df_ordered[, i], "-") }
      outL <- Filter(Negate(is.null), outL)
    }
    if(ord_ind > fac_ind){
      for (i in (fac_ind + 1):ord_ind) { outell[[i - fac_ind]] <- outer(df_ordered[, i], df_ordered[, i], "-") }
      outell <- Filter(Negate(is.null), outell)
    }
    
    # Kernel calculations for continuous variables
    if (con_ind != 0) {
      K <- lapply(1:con_ind, function(i) {
        out <- outK[[i]]
        lambda_val <- lambda[i]
        switch(ckernel,
               "c_gaussian" = (1/sqrt(2*pi)) * exp(-0.5 * (out / lambda_val)^2),
               "c_epanechnikov" = (3/4) * (1 - (out / lambda_val)^2) * (abs(out / lambda_val) <= 1),
               "c_uniform" = 0.5 * (abs(out / lambda_val) <= 1),
               "c_triangle" = (1 - abs(out / lambda_val)) * (abs(out / lambda_val) <= 1),
               "c_biweight" = (15/16) * (1 - (out / lambda_val)^2)^2 * (abs(out / lambda_val) <= 1),
               "c_triweight" = (35/32) * (1 - (out / lambda_val)^2)^3 * (abs(out / lambda_val) <= 1),
               "c_tricube" = (70/81) * (1 - abs(out / lambda_val)^3)^3 * (abs(out / lambda_val) <= 1),
               "c_cosine" = (pi/4) * cos((pi/2) * out / lambda_val) * (abs(out / lambda_val) <= 1),
               "c_logistic" = 1 / (exp(out / lambda_val) + 2 + exp(-out / lambda_val)),
               "c_sigmoid" = 2 / (pi * (exp(out / lambda_val) + exp(-out / lambda_val))),
               "c_silverman" = 0.5 * exp(-abs(out / lambda_val) / sqrt(2)) * sin(abs(out / lambda_val) / sqrt(2) + pi/4)
        )
      })
      D1 <- Reduce(`+`, K)
    } else {
      D1 <- 0
    }
    
    # Kernel calculations for unordered factors
    if (fac_ind > con_ind) {
      L <- lapply(1:length(outL), function(i) {
        out <- outL[[i]]
        lambda_val <- lambda[i + con_ind]
        switch(ukernel,
               "u_aitken" = ifelse(out == 0, 1, lambda_val),
               "u_aitchisonaitken" = ifelse(
                 out == 0,
                 1 - lambda_val,
                 lambda_val / (nominal_categories[i] - 1)
               )
        )
      })
      D2 <- Reduce(`+`, L)
    } else {
      D2 <- 0
    }
    
    # Kernel calculations for ordered factors
    if (ord_ind > fac_ind) {
      ell <- lapply(1:length(outell), function(i) {
        out <- outell[[i]]
        lambda_val <- lambda[i + fac_ind]
        switch(okernel,
               "o_habbema" = lambda_val^(abs(out)^2),
               "o_wangvanryzin" = ifelse(out == 0, 1 - lambda_val, 0.5 * (1 - lambda_val) * lambda_val^abs(out)),
               "o_aitken" = ifelse(out == 0, lambda_val, (1 - lambda_val) / (2^abs(out))),
               "o_aitchisonaitken" = choose(max(unique(out)), abs(out)) * lambda_val^abs(out) * (1 - lambda_val)^(max(unique(out)) - abs(out)),
               "o_liracine" = ifelse(out == 0, 1, lambda_val^abs(out))
        )
      })
      D3 <- Reduce(`+`, ell)
    } else {
      D3 <- 0
    }
    
    D <- D1 + D2 + D3
    diag(D) <- 0
    Fx <- mean(log(colSums(D)) - log(n-1))
    return(Fx)
  }
  max_val <- -Inf
  best_obj <- NULL
  max_fail <- 7
  for (i in seq_len(nstart)) {
    if (verb) message("start ", i, " of ", nstart)
    result <- NULL
    for (attempt in seq_len(max_fail)) {
      params <- runif(N, 1e-16, parameter_upper)
      candidate <- tryCatch(
        optim(
          par = params, mscv_opt,
          control = list(fnscale = -1, maxit = 10000)
        ),
        error = function(e) NULL
      )
      if (!is.null(candidate) && length(candidate$par) == N &&
          all(is.finite(candidate$par))) {
        verified_value <- mscv_opt(candidate$par)
        if (is.finite(verified_value)) {
          candidate$value <- verified_value
          result <- candidate
          break
        }
      }
      if (verb && attempt < max_fail) {
        message("Optimization failed, trying a new starting point...")
      }
    }
    if (is.null(result)) {
      stop(
        "Too many optimization failures; try different kernel functions or ",
        "transform variables that may be on different scales.",
        call. = FALSE
      )
    }
    
    if (result$value > max_val) {
      max_val <- result$value
      best_obj <- result
    }
  }
  if (verb) {
    message(if (best_obj$convergence == 0) {
      "Objective converged"
    } else {
      "Objective did not converge."
    })
  }
  bw <- data.frame(x = best_obj$par)
  rownames(bw) <- colnames(
    df[, c(con_cols, fac_cols, ord_cols), drop = FALSE]
  )
  list(bw = bw, fn_value = best_obj$value)
}
