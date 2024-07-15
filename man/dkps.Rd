\name{dkps}
\alias{dkps}
\title{
Distance using Kernel Product Similarity (DKPS) for Mixed-type Data
}
\description{
This function calculates the pairwise distances between mixed-type observations consisting of numeric (continuous), factor (nominal), and ordered factor (ordinal) variables using the method described in Ghashti, J. S. and Thompson, J. R. J (2023). This kernel metric learning methodology learns the bandwidths associated with each kernel function for each variable type and returns a distance matrix that can be utilized in any distance-based clustering algorithm.
}
\usage{
dkps(df, bw = "mscv", cFUN = "c_gaussian", uFUN = "u_aitken", 
      oFUN = "o_wangvanryzin", stan = TRUE)
}
\arguments{
  \item{df}{
    a \eqn{p}-variate data frame for which the pairwise distances between observations will be calculated. The data types may be continuous, nominal (unordered factors), ordinal (ordered factors), or any combination thereof. Columns of \code{df} should be of appropriate variable type prior to running the function.
}
  \item{bw}{
   a bandwidth specification method. This can be set as a vector of \eqn{p}-many bandwidths, with each element \eqn{i} corresponding to the bandwidth for column \eqn{i} in \code{df}. Alternatively, one of two character strings may be inputted for bandwidth selection methods. \code{mscv} specifies maximum-similarity cross-validation, and \code{np} specifies likelihood-cross validation which is calculated via \code{\link{npudensbw}} in package \code{\link{np}}. Defaults to \code{mscv}. See details
}
  \item{cFUN}{
   character string specifying the continuous kernel function. Options include \code{c_gaussian}, \code{c_epanechnikov}, \code{c_uniform}, \code{c_triangle}, \code{c_biweight}, \code{c_triweight}, \code{c_tricube}, \code{c_cosine}, \code{c_logistic}, \code{c_sigmoid}, and \code{c_silverman}. Note that if using \code{np} for \code{bw} selection above, continuous kernel types are restricted to either \code{c_gaussian}, \code{c_epanechnikov}, or \code{c_uniform}. Defaults to \code{c_gaussian}. See details
}
 \item{uFUN}{
   character string specifying the nominal kernel function for unordered factors. Options include \code{u_aitken} and \code{u_aitchisonaitken}. Defaults to \code{u_aitken}. See details
}
 \item{oFUN}{
   character string specifying the ordinal kernel function for ordered factors. Options include \code{o_aitken}, \code{o_aitchisonaitken}, \code{o_habbema}, \code{o_wangvanryzin}, and \code{o_liracine}. Note that if using \code{np} for \code{bw} selection above, ordinal kernel types are restricted to either \code{o_wangvanryzin} or \code{o_liracine}. Defaults to \code{o_wangvanryzin}. See details
}
 \item{stan}{
   a logical value which specifies whether to scale the resulting distance matrix between 0 and 1 using min-max normalization. If set to \code{FALSE}, there is no normalization. Defaults to \code{TRUE}.
}
}
\details{
  \code{dkps} implements the kernel summation distance (DKPS) as described by Ghashti and Thompson (2023). This approach uses product kernels for continuous variables, and summation kernels for nominal and ordinal data, which are then summed over all variable types to return the pairwise distance between mixed-type data. 
  
  Each kernel requires a bandwidth specification, which can either be a user defined numeric vector of length \eqn{p} from alternative methodologies for bandwidth selection, or through two bandwidth specification methods. The \code{mscv} bandwidth selection routine is based on the maximum-similarity cross-validation routine by Ghashti and Thompson (2023), invoked by the function \code{\link{mscv.dkps}}. The \code{np} bandwidth selection routine follows maximum-likelihood cross-validation techniques described by Li and Racine (2007) and Li and Racine (2003) for kernel density estimation of mixed-type data. Bandwidths will differ for each variable.

  Data contained in the data frame \code{df} may constitute any combinations of continuous, nominal, or ordinal data, which is to be specified in the data frame \code{df} using \code{\link{factor}} for nominal data, and \code{\link{ordered}} for ordinal data. Data can be entered in an arbitrary order and
  data types will be detected automatically. User-inputted vectors of bandwidths \code{bw} must be defined in the same order as the variables in the data frame \code{df}, as to ensure they sorted accordingly by the routine.

  The are many kernels which can be specified by the user. The majority of the continuous kernel functions may be found in Cameron and Trivedi (2005), Härdle et al. (2004) or Silverman (1986). Nominal kernels use a variation on Aitchison and Aitken's (1976) kernel, while ordinal kernels use a variation of
  the Wang and van Ryzin (1981) kernel. Both nominal and ordinal kernel functions can be found in Li and Racine (2007), Li and Racine (2003), Ouyan et al. (2006), and Titterington and Bowman (1985).
}

\value{
  \code{dkps} returns a \code{list} object, with the
  following components:
 \item{distances}{
 an \eqn{n \times n} numeric matrix containing pairwise distances between observations
 }
 \item{bandwidths}{
 a \eqn{p}-variate vector of bandwidth values returned based on the \code{bw} bandwidth specification method, sorted by variable type
 }
}
\references{

  Aitchison, J. and  C.G.G. Aitken (1976), \dQuote{Multivariate binary discrimination by the kernel method,} Biometrika, 63, 413-420.
  
  Cameron, A. and P. Trivedi (2005), \dQuote{Microeconometrics: Methods and Applications}, Cambridge University Press.

  Ghashti, J.S. and J.R.J Thompson (2023), \dQuote{Mixed-type Distance Shrinkage and Selection for Clustering via Kernel Metric Learning.} arXiv preprint arXiv:2306.01890.
  
  Härdle, W., and M. Müller and S. Sperlich and A. Werwatz (2004), \emph{Nonparametric and Semiparametric Models}, (Vol. 1). Berlin: Springer.
  
  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2003), \dQuote{Nonparametric estimation of
   distributions with categorical and continuous data,} Journal
  of Multivariate Analysis, 86, 266-292.

  Ouyang, D. and Q. Li and J.S. Racine (2006), \dQuote{Cross-validation
    and the estimation of probability distributions with categorical
    data,} Journal of Nonparametric Statistics, 18, 69-100.

  Silverman, B.W. (1986), \emph{Density Estimation,} London: Chapman and
  Hall.
  
  Titterington, D.M. and A.W. Bowman (1985), \dQuote{A comparative study of smoothing procedures for ordered categorical data}, Journal of Statistical Computation and Simulation, 21(3-4), 291-312.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
   for discrete distributions,}  Biometrika, 68, 301-309.

}
\author{
    John R. J. Thompson \email{john.thompson@ubc.ca}, Jesse S. Ghashti
    \email{jesse.ghashti@ubc.ca}
}

\seealso{
 \code{\link{mscv.dkps}}, \code{\link{dkss}}, \code{\link{mscv.dkss}}
}
\examples{
\dontrun{
# example data frame with mixed numeric, nominal, and ordinal data.
levels = c("Low", "Medium", "High")
df <- data.frame(
  x1 = runif(100, 0, 100),
  x2 = factor(sample(c("A", "B", "C"), 100, TRUE)),
  x3 = factor(sample(c("A", "B", "C"), 100, TRUE)),
  x4 = rnorm(100, 10, 3),
  x5 = ordered(sample(c("Low", "Medium", "High"), 100, TRUE), levels = levels),
  x6 = ordered(sample(c("Low", "Medium", "High"), 100, TRUE), levels = levels))

# minimal implementation requires just the data frame, and will automatically be
# defaulted to the mscv bandwidth specification technique and default kernel 
# function
d1 <- dkps(df = df)
# d$bandwidths to see the mscv obtained bandwidths
# d$distances to see the distance matrix


# try using the np package, which has few continuous and ordinal kernels to 
# choose from. Recommended using default kernel functions
d2 <- dkps(df = df, bw = "np")


# precomputed bandwidth example
# note that continuous variables requires bandwidths > 0
# ordinal variables requires bandwidths in [0,1]
# for nominal variables, u_aitken requires bandwidths in [0,1] 
# and u_aitchisonaitken in [0,(c-1)/c]
# where c is the number of unique values in the i-th column of df.
# any bandwidths outside this range will result in a warning message
bw_vec <- c(1.0, 0.5, 0.5, 5.0, 0.3, 0.3) 
d3 <- dkps(df = df, bw = bw_vec)


# user-specific kernel functions example
d5 <- dkps(df = df, bw = "mscv", cFUN = "c_epanechnikov", uFUN = "u_aitken", 
      oFUN = "o_habbema")
}
}

\keyword{ ~multivariate }
\keyword{ ~nonparametric }
\keyword{ ~optimize }

\concept{ ~distances }
\concept{ ~metrics }
\concept{ ~metriclearning }


