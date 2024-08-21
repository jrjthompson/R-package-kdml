\name{mscv.dkps}
\alias{mscv.dkps}
\title{
Maximum-similarity Cross-validated (MSCV) bandwidth selection method for the 
Distance using Kernel Product Similarities (DKPS)
}
\description{
This function calculates maximum-similarity cross-validated bandwidths for the 
distance using kernel summation similarity. This implementation uses the method 
described in Ghashti and Thompson (2024) for mixed-type data that includes any 
of numeric (continuous), factor (nominal), and ordered factor (ordinal) 
variables. \code{mscv.dkps} calculates the bandwidths associated with each 
kernel function for variable types and returns a numeric vector of bandwidths 
that can be used with the \code{dkps} pairwise distance calculation.
}

\usage{
mscv.dkps(df, nstart = NULL, ckernel = "c_gaussian", ukernel = "u_aitken",
          okernel = "o_wangvanryzin", verbose = FALSE) 
}
\arguments{
  \item{df}{
    a \eqn{p}-variate data frame. The data types may be continuous (\code{\link{numeric}}), 
    nominal (\code{\link{factor}}), ordinal (\code{\link{ordered}}), or any 
    combination thereof. Columns of \code{df} should be set to the appropriate 
    data type class.
}
  \item{nstart}{
   integer number of restarts for the process of finding extrema of the mscv 
   function from random initial bandwidth parameters (starting points). If the 
   default of \code{NULL} is used, then the number of restarts will be 
   \eqn{min(3,\text{ncol(df)})}.
}
  \item{ckernel}{
   character string specifying the continuous kernel function. Options include 
   \code{c_gaussian}, \code{c_epanechnikov}, \code{c_uniform}, \code{c_triangle}, 
   \code{c_biweight}, \code{c_triweight}, \code{c_tricube}, \code{c_cosine}, 
   \code{c_logistic}, \code{c_sigmoid}, and \code{c_silverman}. Note that if 
   using \code{np} for \code{bw} selection above, continuous kernel types are 
   restricted to either \code{c_gaussian}, \code{c_epanechnikov}, or 
   \code{c_uniform}. Defaults to \code{c_gaussian}. See details.
}
 \item{ukernel}{
   character string specifying the nominal kernel function for unordered factors. 
   Options include \code{u_aitken} and \code{u_aitchisonaitken}. Defaults to 
   \code{u_aitken}. See details.
}
 \item{okernel}{
   character string specifying the ordinal kernel function for ordered factors. 
   Options include \code{o_aitken}, \code{o_aitchisonaitken}, \code{o_habbema}, 
   \code{o_wangvanryzin}, and \code{o_liracine}. Note that if using \code{np} 
   for \code{bw} selection above, ordinal kernel types are restricted to either 
   \code{o_wangvanryzin} or \code{o_liracine}. Defaults to \code{o_wangvanryzin}.
   See details.
}
 \item{verbose}{
   a logical value which specifies whether to output the \eqn{i}-th iteration of
   the total number of \code{nstarts}, and output if the optimization procedure 
   converges. Defaults to \code{FALSE}.
}
}
\details{
  \code{mscv.dkps} implements the maximum-similarity cross-validation (MSCV) 
  technique for bandwidth selection pertaining to the \code{\link{dkps}} function,
  as described by Ghashti and Thompson (2024). This approach uses product kernels 
  for continuous variables, and summation kernels for nominal and ordinal data,
  which are then summed over all variable types to return the pairwise distance 
  between mixed-type data. 
  
  The maximization procedure for bandwidth selection is based on the objective
  \eqn{\text{arg}\max_{\boldsymbol{\lambda}}\left\{\frac{1}{n}\sum_{i=1}^n\log\left(\frac{1}{(n-1)}\sum_{\substack{j=1 \\ j \ne i}}^n\psi_{\boldsymbol{\lambda}}(\textbf{x}_i,\textbf{x}_j)\right)\right\},}
  where
  
  \eqn{\psi(\textbf{x}_i, \textbf{x}_j \ | \boldsymbol{\lambda}) = \prod_{k=1}^{p_c}\frac{1}{\lambda_k^c}K(x_{i,k}^c, x_{j,k}^c, \lambda_k^c) + \sum_{k=1}^{p_u}L(x_{i,k}^u,x_{j,k}^u,\lambda_k^u) + \sum_{k=1}^{p_o}\ell(x_{i,k}^o,x_{j,k}^o,\lambda_k^o).}
  
  \eqn{K(\cdot)}, \eqn{L(\cdot)}, and \eqn{\ell(\cdot)} are the continuous, 
  nominal, and ordinal kernel functions, repectively, with \eqn{\lambda_k}'s 
  representing kernel specifical bandwiths for the \eqn{k}-th variable, and 
  \eqn{p_c}, \eqn{p_u}, \eqn{p_o} the number of continuous, nominal, and ordinal
  variables in the data frame \code{df}. The resulting \code{bw} vector returned
  is the bandwidths that yield the highest objective function value. 

  Data contained in the data frame \code{df} may constitute any combinations of 
  continuous, nominal, or ordinal data, which is to be specified in the data 
  frame \code{df} using \code{\link{numeric}} for continuous data, \code{\link{factor}}
  for nominal data, and \code{\link{ordered}} for ordinal data. Data can be 
  entered in an arbitrary order and data types will be detected automatically. 
  User-inputted vectors of bandwidths \code{bw} must be defined in the same 
  order as the variables in the data frame \code{df}, as to ensure they sorted 
  accordingly by the routine.

  The are many kernels which can be specified by the user. Continuous kernel 
  functions may be found in Cameron and Trivedi (2005), Härdle et al. (2004) or 
  Silverman (1986). Nominal kernels use a variation on Aitchison and Aitken's 
  (1976) kernel. Ordinal kernels use a variation of the Wang and van Ryzin 
  (1981) kernel. All nominal and ordinal kernel functions can be found in Li and
  Racine (2007), Li and Racine (2003), Ouyan et al. (2006), and Titterington and
  Bowman (1985).
}

\value{
  \code{mscv.dkps} returns a \code{list} object, with the
  following components:
 \item{bw}{a \eqn{p}-variate vector of bandwidth values, intended to be used for 
 the \code{\link{dkps}} pairwise distance calculation }
 \item{fn_value}{a numeric value of the MSCV objective function, obtained using 
 the \code{\link{optim}} function for constrained multivariate optimization}
}
\references{

  Aitchison, J. and  C.G.G. Aitken (1976), \dQuote{Multivariate binary 
  discrimination by the kernel method,} Biometrika, 63, 413-420.
  
  Cameron, A. and P. Trivedi (2005), \dQuote{Microeconometrics: Methods and 
  Applications}, Cambridge University Press.

  Ghashti, J.S. and J.R.J Thompson (2024), \dQuote{Mixed-type Distance Shrinkage
  and Selection for Clustering via Kernel Metric Learning.} Journal of 
  Classification, Accepted.
  
  Härdle, W., and M. Müller and S. Sperlich and A. Werwatz (2004), 
  \emph{Nonparametric and Semiparametric Models}, (Vol. 1). Berlin: Springer.
  
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
  
  Titterington, D.M. and A.W. Bowman (1985), \dQuote{A comparative study of 
  smoothing procedures for ordered categorical data}, Journal of Statistical 
  Computation and Simulation, 21(3-4), 291-312.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
   for discrete distributions,}  Biometrika, 68, 301-309.

}
\author{
    John R. J. Thompson \email{john.thompson@ubc.ca}, Jesse S. Ghashti
    \email{jesse.ghashti@ubc.ca}
}


\seealso{
 \code{\link{mscv.dkss}}, \code{\link{dkss}}, \code{\link{dkps}}
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

# minimal implementation requires just the data frame, with defaults
bw <- mscv.dkps(df = df)

# specify number of starts and kernel functions
bw2 <- mscv.dkps(df = df, nstart = 5, ckernel = "c_triangle",
                  ukernel = "u_aitken", okernel = "o_liracine")
}
}

\keyword{ ~multivariate }
\keyword{ ~nonparametric }
\keyword{ ~optimize }

\concept{ ~distances }
\concept{ ~metrics }
\concept{ ~metriclearning }

