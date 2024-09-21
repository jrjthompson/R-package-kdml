\name{kss}
\alias{kss}
\title{
Kernel Summation Similarity Function (KSS) for Mixed-type Data
}
\description{
This function calculates the pairwise similarities between mixed-type 
observations  consisting of continuous (\code{numeric}), nominal 
(\code{factor}), and ordinal (\code{ordered}) variables using the method 
described in Ghashti (2024). This kernel similarity learning methodology 
calculates a kernel sum similarity function, with a variety of options for 
kernel functions associated with each variable type and returns a distance 
matrix that can be used in any distance-based algorithm.
}
\usage{
kss(df, bw = "np", npmethod = NULL, cFUN = "c_gaussian", 
    uFUN = "u_aitken", oFUN = "o_wangvanryzin", nstart = NULL, 
    stan = TRUE, verbose = FALSE)
}
\arguments{
  \item{df}{
    a \eqn{p}-variate data frame for which the pairwise similarities between 
    observations will be calculated. The data types may be continuous 
    (\code{numeric}), nominal (\code{factor}), and ordinal (\code{ordered}), 
    or any combination thereof. Columns of \code{df} should be of appropriate 
    variable type prior to running the function.
}
  \item{bw}{
   numeric bandwidth vector of length \eqn{p}, with each element \eqn{i} 
   corresponding to the bandwidth for column \eqn{i} in \code{df}. 
   Alternatively, a character strings may be inputted for bandwidth 
   selection methods. \code{np} specifies this techniques which calculate 
   bandwidths using \code{\link[np]{npudensbw}} in package \code{\link[np]{np}}. 
   Defaults to \code{np} with \code{npmethod} set to \code{cv.ml} for maximum 
   likelihood cross-validation. See details.
}

  \item{npmethod}{
   character value specifying the \code{np} bandwidth selection to be used for 
   calculating bandwidths. Options include \code{cv.ml} for maximum likelihood 
   cross-validation, \code{cv.ls} for least squares cross-validation, and 
   \code{normal-reference} for normal reference. If left as \code{NULL}, defaults 
   to \code{cv.ml}.
}

  \item{cFUN}{
   character value specifying the continuous kernel function. Options include 
   \code{c_gaussian}, \code{c_epanechnikov}, \code{c_uniform}, 
   \code{c_triangle}, \code{c_biweight}, \code{c_triweight}, \code{c_tricube}, 
   \code{c_cosine}, \code{c_logistic}, \code{c_sigmoid}, and \code{c_silverman}. 
   Note that if using \code{np} for \code{bw} selection above, continuous kernel 
   types are restricted to either \code{c_gaussian}, \code{c_epanechnikov}, or
   \code{c_uniform}. Defaults to \code{c_gaussian}. See details.
}
 \item{uFUN}{
   character value specifying the nominal kernel function for unordered 
   factors. Options include \code{u_aitken} and \code{u_aitchisonaitken}. 
   Defaults to \code{u_aitken}. See details.
}
 \item{oFUN}{
   character value specifying the ordinal kernel function for ordered factors.
   Options include \code{o_aitken}, \code{o_aitchisonaitken}, \code{o_habbema}, 
   \code{o_wangvanryzin}, and \code{o_liracine}. Note that if using \code{np} 
   for \code{bw} selection above, ordinal kernel types are restricted to either
   \code{o_wangvanryzin} or \code{o_liracine}. Defaults to \code{o_wangvanryzin}. 
   See details.
}
  \item{nstart}{
    integer value specifying the number of random starts for the \code{kmeans} 
    algorithm. Defaults to \code{10}.
}
\item{stan}{
   a logical value which specifies whether to scale the resulting distance 
   matrix between 0 and 1 using min-max normalization. If set to \code{FALSE}, 
   distances are unscaled. Defaults to \code{TRUE}.
}
 \item{verbose}{
   a logical value which specifies whether to print procedural steps to the 
   console. If set to \code{FALSE}, no output is printed to the console. 
   Defaults to \code{FALSE}.
}
}
\details{
  \code{kss} implements the kernel summation similarity function (KSS)
  as described by Ghashti (2024). This approach uses summation kernels for 
  continuous, nominal and ordinal data, which are then summed over all variable 
  types to return the pairwise similarities between mixed-type data. 
  
  There are several kernels to select from. The continuous kernel functions may 
  be found in Cameron and Trivedi (2005), Härdle et al. (2004) or Silverman 
  (1986). Nominal kernels use a variation on Aitchison and Aitken's (1976) 
  kernel, while ordinal kernels use a variation of the Wang and van Ryzin (1981)
  kernel. Both nominal and ordinal kernel functions can be found in Li and 
  Racine (2007), Li and Racine (2003), Ouyan et al. (2006), and Titterington and
  Bowman (1985).
  
  Each kernel requires a bandwidth specification, which can either be a user 
  defined numeric vector of length \eqn{p} from alternative methodologies for 
  bandwidth selection, or through one bandwidth selection method can be 
  specified. The \code{np} bandwidth selection methods follow three techniques 
  (\code{cv.ml}, \code{cv.ls} and \code{normal-reference}) described by Li and 
  Racine (2007) and Li and Racine (2003) for kernel density estimation of 
  mixed-type data. 

  Data contained in the data frame \code{df} may constitute any combinations of 
  continuous, nominal, or ordinal data, which is to be specified in the data 
  frame \code{df} using \code{\link{factor}} for nominal data, and 
  \code{\link{ordered}} for ordinal data. Data types can be in any order and 
  will be detected automatically. User-inputted vectors of
  bandwidths \code{bw} must be specified in the same order as the variables in 
  the data frame \code{df}, as to ensure they sorted accordingly by the routine.
}

\value{
  \code{kss} returns a \code{list} object, with the
  following components:
 \item{similarities}{
 an \eqn{n \times n} numeric matrix containing pairwise similarities between 
 observations
 }
 \item{bandwidths}{
 a \eqn{p}-variate vector of bandwidth values returned based on the \code{bw} 
 bandwidth specification method, sorted by variable type
 }
}
\references{

  Aitchison, J. and  C.G.G. Aitken (1976), \dQuote{Multivariate binary 
  discrimination by the kernel method}, Biometrika, 63, 413-420.
  
  Cameron, A. and P. Trivedi (2005), \dQuote{Microeconometrics: Methods and 
  Applications}, Cambridge University Press.
  
  Ghashti, J.S. (2024), \dQuote{Similarity Maximization and Shrinkage Approach 
  in Kernel Metric Learning for Clustering Mixed-type Data (T)}, University of 
  British Columbia.
  
  Härdle, W., and M. Müller and S. Sperlich and A. Werwatz (2004), 
  \dQuote{Nonparametric and Semiparametric Models}, (Vol. 1). Berlin: Springer.
  
  Li, Q. and J.S. Racine (2007), \dQuote{Nonparametric Econometrics: Theory
  and Practice}, Princeton University Press.

  Li, Q. and J.S. Racine (2003), \dQuote{Nonparametric estimation of
  distributions with categorical and continuous data}, Journal
  of Multivariate Analysis, 86, 266-292.

  Ouyang, D. and Q. Li and J.S. Racine (2006), \dQuote{Cross-validation
  and the estimation of probability distributions with categorical
  data}, Journal of Nonparametric Statistics, 18, 69-100.

  Silverman, B.W. (1986), \dQuote{Density Estimation}, London: Chapman and
  Hall.
  
  Titterington, D.M. and A.W. Bowman (1985), \dQuote{A comparative study of 
  smoothing procedures for ordered categorical data}, Journal of Statistical 
  Computation and Simulation, 21(3-4), 291-312.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
  for discrete distributions},  Biometrika, 68, 301-309.

}
\author{
    John R. J. Thompson \email{john.thompson@ubc.ca}, Jesse S. Ghashti
    \email{jesse.ghashti@ubc.ca}
}


\seealso{
 \code{\link{mscv.dkps}}, \code{\link{dkps}}, \code{\link{mscv.dkss}}, \code{\link{dkss}}, \code{link{spectral.clust}}
}
\examples{
\donttest{
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
s1 <- kss(df = df)
# s$bandwidths to see the mscv obtained bandwidths
# s$similarities to see the similarity matrix


# try using the np package, which has few continuous and ordinal kernels 
# to choose from. Recommended using default kernel functions
s2 <- kss(df = df, bw = "np") #defaults to npmethod "cv.ml"


# precomputed bandwidth example
# note that continuous variables requires bandwidths > 0
# ordinal variables requires bandwidths in [0,1]
# for nominal variables, u_aitken requires bandwidths in [0,1] 
# and u_aitchisonaitken in [0,(c-1)/c]
# where c is the number of unique values in the i-th column of df.
# any bandwidths outside this range will result in a warning message
bw_vec <- c(1.0, 0.5, 0.5, 5.0, 0.3, 0.3) 
s3 <- kss(df = df, bw = bw_vec)


# user-specific kernel functions example with "cv.ls" from np.
s4 <- kss(df = df, bw = "np", npmethod = "cv.ls", cFUN = "c_epanechnikov", 
    uFUN = "u_aitken", oFUN = "o_wangvanryzin")
}
}

\keyword{ ~multivariate }
\keyword{ ~nonparametric }
\keyword{ ~optimize }

\concept{ ~distances }
\concept{ ~metrics }
\concept{ ~metriclearning }