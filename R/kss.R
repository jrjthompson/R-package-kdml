kss <- function(df, bw = "np", npmethod = NULL, cFUN = "c_gaussian", uFUN = "u_aitken", oFUN = "o_wangvanryzin", nstart = NULL, stan = TRUE, verbose = FALSE) {
  if (is.null(nstart)) nstart = min(3, ncol(df))
  
  v_ck <- c("c_gaussian", "c_epanechnikov", "c_uniform", "c_triangle",
            "c_biweight", "c_triweight", "c_tricube", "c_cosine", 
            "c_logistic", "c_sigmoid", "c_silverman")
  v_uk <- c("u_aitken", "u_aitchisonaitken")
  v_ok <- c("o_wangvanryzin", "o_habbema", "o_aitken", "o_aitchisonaitken", "o_liracine")
  
  if (!(cFUN %in% v_ck)) stop("Invalid cFUN specified. Choose one of: ", paste(v_ck, collapse = ", "))
  if (!(uFUN %in% v_uk)) stop("Invalid uFUN specified. Choose one of: ", paste(v_uk, collapse = ", "))
  if (!(oFUN %in% v_ok)) stop("Invalid oFUN specified. Choose one of: ", paste(v_ok, collapse = ", "))
  
  # Get column names by type
  con_cols <- names(df)[sapply(df, is.numeric)]
  fac_cols <- names(df)[sapply(df, function(x) is.factor(x) & !is.ordered(x))]
  ord_cols <- names(df)[sapply(df, is.ordered)]
  
  # Reorder the dataframe
  df_ordered <- df[, c(con_cols, fac_cols, ord_cols)]
  
  # Bandwidth selection options (lines 52 - 77)
  if (is.null(bw)) stop("No bandwidth selection chosen. Input a numeric vector of bandwidths, or choose from 'np' or 'mscv' bandwidth selection methods")
  if (is.numeric(bw)) {
    bws <- bw[match(c(con_cols, fac_cols, ord_cols), names(df))]
  } else if (bw == "np") {
    if (is.null(npmethod)) {
      if (verbose == TRUE) warning("No npmethod selected, defaulting to maximum likelihood cross validation (cv.ml). You may also least squares cross-validation (cv.ls), normal reference (normal-reference)")
      npmethod = "cv.ml"
    } else if (!is.null(npmethod) && (npmethod != "cv.ml" && npmethod != "cv.ls" && npmethod != "normal-reference")) {
      stop("npmethod chosen does not exist. Select one of maximum likelihood cross validation (cv.ml), least squares cross-validation (cv.ls) or normal reference (normal-reference)")
    }
    if (cFUN == "c_gaussian") cker <- "gaussian"
    if (cFUN == "c_epanechnikov") cker <- "epanechnikov"
    if (cFUN == "c_uniform") cker <- "uniform"
    if (cFUN %in% c("c_triangle", "c_biweight", "c_triweight", "c_tricube", "c_cosine", "c_logistic", "c_sigmoid", "c_silverman")) {
      stop("Choose one of c_gaussian, c_epanechnikov, or c_uniform for continuous kernels while using np") 
    }
    if (uFUN == "u_aitken") uker <- "liracine"
    if (uFUN == "u_aitchisonaitken") uker <- "aitchisonaitken"
    if (oFUN == "o_wangvanryzin") oker <- "wangvanryzin"
    if (oFUN == "o_liracine") oker <- "liracine"
    if (oFUN %in% c("o_habbema", "o_aitken", "o_aitchisonaitken")) {
      stop("Choose one of o_wangvanryzin or o_liracine for ordinal kernels while using np")
    }
    bws <- npudensbw(df_ordered, ckertype = cker, ukertype = uker, okertype = oker, bwmethod = npmethod, nmulti = nstart)$bw
    if (verbose == TRUE) print("Bandwidth calculation complete. Computing similarities.")
  } else {
    stop("Invalid bandwidth selection. Input either a numeric vector of bandwidths, or choose 'np' or 'mscv'.")
  }
  
  # store indices
  con_ind <- length(con_cols) #index of continuous variables
  fac_ind <- ifelse(length(fac_cols) > 0, length(con_cols) + length(fac_cols), 0) #index of factors
  ord_ind <- ifelse(length(ord_cols) > 0, length(con_cols) + length(fac_cols) + length(ord_cols), 0) #index of ordinal variables
  
  df <- df_ordered
  df <- data.matrix(df)
  
  if(length(bws) != ncol(df)) stop("Invalid length of bandwidth vector (bw).")
  similarities <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  c1 <- combn(nrow(df), 2)
  c2 <- rbind(1:nrow(df),1:nrow(df))
  combinations <- cbind(c1,c2)
  # All available kernel functions
  kernel_functions <- list(
    c_gaussian = function(A, B, bws, C) {sum(exp(-((A - B) / bws)^2 / 2))}, #Continuous Gaussian kernel
    c_uniform = function(A, B, bws, C) {sum(sapply(1:length(A), function(i) ifelse(abs((A[i] - B[i]) / bws[i]) <= 1, 1/2, 0)))}, #Continuous Uniform kernel
    c_epanechnikov = function(A, B, bws, C) {sum(ifelse(abs((A - B) / bws) <= 1, (3/4) * (1 - ((A - B) / bws)^2), 0))}, #Continuous Epanechnikov kernel
    c_triangle = function(A, B, bws, C) {sum(ifelse(abs((A - B) / bws) <= 1, (1 - abs((A - B) / bws)), 0))}, #Continuous Triangle kernel
    c_biweight = function(A, B, bws, C) {sum(ifelse(abs((A - B) / bws) <= 1, (15/32) * (3 - ((A - B) / bws)^2)^2, 0))}, #Continuous biweight kernel
    c_triweight = function(A, B, bws, C) {sum(ifelse(abs((A - B) / bws) <= 1, (35/32) * (1 - ((A - B) / bws)^2)^3, 0))}, #Continuous triweight kernel
    c_tricube = function(A, B, bws, C) {sum(ifelse(abs((A - B) / bws) <= 1, (70/81) * (1 - abs((A - B) / bws)^3)^3, 0))}, #Continuous tricube kernel
    c_cosine = function(A, B, bws, C) {sum(ifelse(abs((A - B) / bws) <= 1, (pi/4) * cos((pi/2) * ((A - B) / bws)), 0))}, #Continuous cosine kernel
    c_logistic = function(A, B, bws, C) {sum(1 / (exp((A - B) / bws) + 2 + exp(-((A - B) / bws))))}, #Continuous logistic kernel
    c_sigmoid = function(A, B, bws, C) {sum(2 / (pi * (exp((A - B) / bws) + exp(-((A - B) / bws)))))}, #Continuous Sigmoid kernel
    c_silverman = function(A, B, bws, C) {sum(0.5 * exp(-abs((A-B)/bws) / sqrt(2)) * sin(abs((A-B)/bws) / sqrt(2) + pi/4))}, #Continuous Silverman Kernel
    u_aitchisonaitken = function(A, B, bws, C) {sum(ifelse(A == B, 1 - bws, bws / (length(unique(C[,i])) - 1)))}, #Nominal Aitchison and Aitken kernel
    u_aitken = function(A, B, bws, C) {sum(ifelse(A == B, 1, bws))}, #Nominal Aitken kernel
    o_wangvanryzin = function(A, B, bws, C) {sum(ifelse(all(A == B), 1 - bws, (1/2) * (1 - bws) * (bws^abs(A - B))))}, #Ordinal Wang & van Ryzin kernel
    o_aitchisonaitken = function(A, B, bws, C) {sum(ifelse(A == B, 1, bws^(abs(A - B))))}, #Ordinal Aitchison and Aitken kernel
    o_aitken = function(A, B, bws, C) {sum(ifelse(A == B, bws, (1 - bws) / (2^abs(A - B))))}, #Ordinal Aitken kernel
    o_habbema = function(A, B, bws, C) {sum(bws^((abs(A - B))^2))}, #Ordinal Habbema kernel
    o_liracine = function(A, B, bws, C) {sum(ifelse(A == B, 1, bws^abs(A-B)))} #Ordinal Li & Racine Kernel
  )
  
  kernel.distance.sum <- function(con_ind, fac_ind, ord_ind, cFUN, uFUN, oFUN, A, B, bws, df) {
    get_function <- function(FUN_name, type) {
      if (!FUN_name %in% names(kernel_functions)) {
        stop(paste("Invalid", type, "function name:", FUN_name))
      }
      return(kernel_functions[[FUN_name]])
    }
    
    compute_distance <- function(FUN, ind) {
      if (length(ind) == 0) return(0)
      FUN(A[ind], B[ind], bws[ind], df[, ind])
    }
    d_c <- if (con_ind > 0) compute_distance(get_function(cFUN, "cFUN"), 1:con_ind) else 0
    d_u <- if (fac_ind > 0) compute_distance(get_function(uFUN, "uFUN"), (con_ind + 1):fac_ind) else 0
    d_o <- if (ord_ind > 0) compute_distance(get_function(oFUN, "oFUN"), (con_ind + fac_ind + 1):ord_ind) else 0
    return(d_c + d_u + d_o)
  }
  for (i in 1:ncol(combinations)) {
    row1 <- combinations[1, i]
    row2 <- combinations[2, i]
    distance <- kernel.distance.sum(con_ind, fac_ind, ord_ind, cFUN, uFUN, oFUN, df[row1, ], df[row2, ], bws, df)
    similarities[row1, row2] <- as.numeric(distance)
    similarities[row2, row1] <- as.numeric(distance)
  }
  bandwidths <- matrix(bws, 1, length(bws))
  colnames(bandwidths) <- colnames(df_ordered)
  if (stan == FALSE) {
    if (verbose == TRUE) print("Completed distance calculation.")
    return(list(similarities = similarities, bandwidths = bandwidths))
  }
  if (stan == TRUE) {
    min <- min(similarities)
    max <- max(similarities)
    standardized <- (similarities - min) / (max - min)
    if (verbose == TRUE) print("Completed distance calculation.")
    return(list(similarities = standardized, bandwidths = bandwidths))
  }
}

