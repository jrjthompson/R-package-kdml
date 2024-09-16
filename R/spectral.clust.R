spectral.clust <- function(S, k, nstart = 10, iter.max = 1000, is.sim = NULL, neighbours = 10){ 
  if (!is.matrix(S) || nrow(S) != ncol(S)) stop("S must be a square matrix")
  if (k <= 0 || k > nrow(S)) stop("Invalid value for k")
  if (neighbours < 1 || neighbours >= nrow(S)) stop("Invalid number of neighbours")
  
  D <- matrix(0, nrow = nrow(S), ncol = nrow(S))  
  if (is.null(is.sim)) stop("Must specify is.sim as TRUE for similarity matrix and FALSE for distance matrix")
  
  #adjacency matrix
  if (is.sim == TRUE) { 
    for (i in 1:nrow(S)) { 
      index <- order(S[i,], decreasing = TRUE)[2:(min(neighbours + 1, nrow(S)))] 
      D[i,][index] <- 1  
    } 
  } else if (is.sim == FALSE) { 
    for (i in 1:nrow(S)) { 
      index <- order(S[i,], decreasing = FALSE)[2:(min(neighbours + 1, nrow(S)))] 
      D[i,][index] <- 1  
    } 
  } 
  
  D <- D + t(D)  
  D[D == 2] <- 1 
  degrees <- colSums(D)  
  n <- nrow(D) 
  laplacian <- diag(n) - diag(degrees^(-1/2)) %*% D %*% diag(degrees^(-1/2)) 
  
  # Eigen decomposition
  eigen_result <- eigen(laplacian, symmetric = TRUE) 
  eigenvectors <- eigen_result$vectors[, (n-k+1):n]  # Use sorted eigenvalues if needed
  
  #cluster
  sc <- kmeans(eigenvectors, k, iter.max = iter.max, nstart = nstart) 
  return(list(clusters = sc$cluster, S = S)) 
}
