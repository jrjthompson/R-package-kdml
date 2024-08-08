confactord <- function(n = 200, 
                        popProb = c(0.5,0.5), 
                        numMixVar = c(1,1,1), 
                        numMixVarOl = c(1,1,1),  
                        olVarType = c(0.1,0.1,0.1), 
                        catLevels = c(2,4)) {
  # storing each entry for use
  pops = popProb
  dimCon = numMixVar[1]
  dimCat = numMixVar[2]
  dimOrd = numMixVar[3]
  dimConOl = numMixVarOl[1]
  dimCatOl = numMixVarOl[2]
  dimOrdOl = numMixVarOl[3]
  conOl = olVarType[1]
  catOl = olVarType[2]
  ordOl = olVarType[3]
  catLev = catLevels[1]
  ordLev = catLevels[2]
  
  # argument checking
  # only 2 populations
  if (length(popProb) != 2) stop('Error: Must be 2 elements in popProb vec')
  if (sum(popProb) != 1) stop('Error: Elements of popProb vec must sum to 1')
  if (dimConOl > dimCon) stop('Error: # continuous variables must be >= # continuous variables with overlap')
  if (dimCatOl > dimCat) stop('Error: # nominal variables must be >= # nominal variables with overlap')
  if (dimOrdOl > dimOrd) stop('Error: # ordinal variables must be >= # ordinal variables with overlap')
  if (conOl > 0.99 | conOl < 0.01 | catOl > 0.99 | catOl < 0.01 | ordOl > 0.99 | ordOl < 0.01) stop('Error: 
                              Overlap must be between 0.01 and 0.99 (inclusive).')
  if (catLev < 2 | ordLev < 2)  warning('There should be at least 2 levels for 
                                       nominal and ordinal variables in catLevels argument')
  
  dimConOl = rep(c(1,0), c(dimConOl,dimCon - dimConOl))
  dimOrdOl = rep(c(1,0), c(dimOrdOl,dimOrd - dimOrdOl))
  ol = 1e-2 # overlap defaults
  memb <- numeric(n)
  exco <- round(n * pops)
  diffco <- n - sum(exco)
  if (diffco != 0) {
    if (diffco > 0) {
      exco[which.max(exco)] <- exco[which.max(exco)] + diffco
    } else {
      exco[which.min(exco)] <- exco[which.min(exco)] + diffco
    }
  }
  for (i in 1:length(pops)) memb[sample(which(memb == 0), exco[i])] <- i

  centres = c(0,-2*qnorm(ol/2))
  cV = centres[memb]
  # Continuous variables
  if(dimCon > 0){
    # generate continuous variables
    conV = matrix(rnorm(n*dimCon,mean=cV,sd=1), nrow=n, ncol=dimCon)
    # consider user inputted cluster overlap for continuous variables and calculate variance for overlaps
    for (i in 1:dimCon) if (dimConOl[i] == 1) conV[,i] = conV[,i] +  rnorm(n, 0, sqrt((abs(diff(centres)) / (2 * qnorm(conOl/2)))^2 - 1))
    conV <- as.data.frame(conV)
  } else conV <- NULL
  
  # Nominal variable generation
  if (dimCat > 0) {
    # generate Dirichlet distribution probabilities
    dirichlet <- function(alpha) {
      samp <- rgamma(length(alpha), alpha, 1)
      samp / sum(samp)
    }
    
    # generate no cluster overlap nominal variables if specified
    if (dimCatOl < dimCat) {
      # population probabilities using Dirichlet distribution with specified alpha parameters
      popProb1 = dirichlet(rep((1 - ol) / catLev + ol / (catLev * 2), catLev))
      popProb2 = dirichlet(rep(ol / catLev + (1 - ol) / (catLev * 2), catLev))
      cNoErr = rep(NA, n * (dimCat - dimCatOl))
      memNoErr = rep(memb, dimCat - dimCatOl)
      cNoErr[memNoErr == 1] = sample(1:catLev, sum(memb == 1) * (dimCat - dimCatOl), TRUE, popProb1)
      cNoErr[memNoErr == 2] = sample(1:catLev, (n - sum(memb == 1)) * (dimCat - dimCatOl), TRUE, popProb2)
      cNoErr = matrix(cNoErr, ncol = dimCat - dimCatOl)
    } else cNoErr = c()
    
    # generate cluster overlap on nominal variables if specified
    if (dimCatOl > 0) {
      # population probabilities using Dirichlet distribution with specified alpha parameters AND user defined overlap
      p1Err = dirichlet(rep((1 - catOl) / catLev + catOl / (catLev * 2), catLev))
      p2Err = dirichlet(rep(catOl / catLev + (1 - catOl) / (catLev * 2), catLev))
      cwErr = rep(NA, n * dimCatOl)
      memwErr = rep(memb, dimCatOl)
      cwErr[memwErr == 1] = sample(1:catLev, sum(memb == 1) * dimCatOl, TRUE, p1Err)
      cwErr[memwErr == 2] = sample(1:catLev, (n - sum(memb == 1)) * dimCatOl, TRUE, p2Err)
      cwErr = matrix(cwErr, ncol = dimCatOl)
    } else cwErr = c()

    catV = as.data.frame(cbind(cNoErr, cwErr))
    for (i in 1:ncol(catV)) catV[, i] <- as.factor(catV[, i])
  } else catV <- NULL
  
  # Ordinal variables
  if(dimOrd > 0){
    # generate continuous variables
    ordV = matrix(rnorm(n*dimOrd,mean=cV,sd=1), nrow=n, ncol=dimOrd)

    # add cluster overlap if specified and convert into ordinal variables based on their quantile distributions
    for (i in 1:dimOrd) {
      cVal <- ordV[,i]
      if (dimOrdOl[i] == 1) {
        oErr <- 
        cVal <- cVal + rnorm(n, 0, sqrt((abs(diff(centres)) / (2 * qnorm(ordOl/2)))^2 - 1))
      }
      ordLevels <- quantile(c(cVal-1e-5, cVal+1e-5), probs = seq(0, 1, length.out = ordLev + 1))
      ordV[,i] = cut(cVal, breaks = ordLevels, labels = FALSE)
    }
    ordV <- as.data.frame(ordV)
    for(i in 1:ncol(ordV)) ordV[,i] <- as.ordered(ordV[,i])
  } else ordV <- NULL
  
  # configure
  mat <- list(conV, catV, ordV)
  nnmat <- mat[!sapply(mat, is.null)]
  comb <- do.call(cbind, nnmat)
  colnames(comb) <- paste0("V", 1:ncol(comb))
  return(list(data = comb, class = memb))
}

