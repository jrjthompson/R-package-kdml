\name{dkss}
\alias{dkss}
\title{
Distance using Kernel Summation Similarity (DKSS) for Mixed-type Data
}
\description{
This function calculates the pairwise distances between mixed-type observations 
consisting of continuous (\code{numeric}), nominal (\code{factor}), and ordinal 
(\code{ordered}) variables using the method described in Ghashti (2024). 
This kernel metric learning methodology calculates a kernel sum similarity 
function, with a variety of options for kernel functions associated with each 
variable type and returns a distance matrix that can be used in any distance-
based algorithm.
}
\usage{
dkss(df, bw = "mscv", cFUN = "c_gaussian", uFUN = "u_aitken",
    oFUN = "o_wangvanryzin", stan = TRUE, verbose = FALSE)
}
\arguments{
  \item{df}{
    a \eqn{p}-variate data frame for which the pairwise distances between
    observations will be calculated, or a \code{kdml_mcmc} object returned by
    \code{\link{mcmc.dkss}}. Data-frame columns may be continuous
    (\code{numeric}), nominal (\code{factor}), ordinal (\code{ordered}), or any
    combination thereof.
}
  \item{bw}{
   numeric bandwidth vector of length \eqn{p}, with each element \eqn{i} 
   corresponding to the bandwidth for column \eqn{i} in \code{df}. 
   Alternatively, one of two character strings may be inputted for bandwidth 
   selection methods. \code{mscv} specifies maximum-similarity cross-validation,
   and \code{np} specifies likelihood-cross validation which is calculated using 
   \code{\link[np]{npudensbw}} in package \code{\link[np]{np}}. Defaults to 
   \code{mscv}. For joint MCMC selection of feature-specific kernels and
   bandwidths, see \code{\link{mcmc.dkss}}.
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
  \code{dkss} implements the distance using summation similarity distance (DKSS)
  as described by Ghashti (2024). This approach uses summation kernels for 
  continuous, nominal and ordinal data, which are then summed over all variable 
  types to return the pairwise distance between mixed-type data. 
  
  There are several kernels to select from. The continuous kernel functions may 
  be found in Cameron and Trivedi (2005), Härdle et al. (2004) or Silverman 
  (1986). Nominal kernels use a variation on Aitchison and Aitken's (1976) 
  kernel, while ordinal kernels use a variation of the Wang and van Ryzin (1981)
  kernel. Both nominal and ordinal kernel functions can be found in Li and 
  Racine (2007), Li and Racine (2003), Ouyan et al. (2006), and Titterington and
  Bowman (1985).
  
  Each kernel requires a bandwidth specification. Bandwidths may be supplied
  as a user-defined numeric vector of length \eqn{p}, or selected using either
  of the two built-in \code{bw} methods. The \code{mscv} method is based on
  maximum-similarity cross-validation by Ghashti and Thompson (2024), invoked
  by \code{\link{mscv.dkss}}. The \code{np} method follows the maximum-
  likelihood cross-validation techniques described by Li and Racine (2007)
  and Li and Racine (2003) for kernel density estimation of mixed-type data.
  \code{\link{mcmc.dkss}} uses Markov chain Monte Carlo to jointly select
  feature-specific kernels and bandwidths. Its fitted object can be passed
  directly to \code{dkss}, which uses the sampled MAP state by default.

  Data contained in the data frame \code{df} may constitute any combinations of 
  continuous, nominal, or ordinal data, which is to be specified in the data 
  frame \code{df} using \code{\link{factor}} for nominal data, and 
  \code{\link{ordered}} for ordinal data. Data types can be in any order and 
  will be detected automatically. Unnamed bandwidth vectors must follow the
  column order of \code{df}. Named bandwidth vectors may be supplied in any
  order but must contain every feature name exactly once. Returned bandwidths
  are ordered by variable type.

}

\value{
  \code{dkss} returns a \code{list} object, with the
  following components:
 \item{distances}{
 an \eqn{n \times n} numeric matrix containing pairwise distances between 
 observations
 }
 \item{bandwidths}{
 a \eqn{p}-variate vector of bandwidth values returned based on the \code{bw} 
 bandwidth specification method, sorted by variable type
 }
 \item{kernels}{
 a named character vector of selected feature-specific kernels. This component
 is present only when an MCMC fit is supplied
 }
 \item{state}{
 selected retained draw and chain indices. This component is present only when
 \code{df} is an MCMC fit
 }
 \item{standardized}{
 logical recording \code{stan}. This component is present only when an MCMC
 fit is supplied
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
    \email{jesse.ghashti@ubc.ca}, and Daniel Krasnov
    \email{daniel.krasnov@mail.mcgill.ca}
}


\seealso{
 \code{\link{mscv.dkss}}, \code{\link{mcmc.dkss}},
 \code{\link{dkps}},
 \code{\link{mscv.dkps}}, \code{\link{mcmc.dkps}}
}
\examples{
set.seed(1)
levels <- c("Low", "Medium", "High")
df <- data.frame(
  x = seq_len(15),
  group = factor(rep(c("A", "B", "C"), 5)),
  rating = ordered(rep(levels, 5), levels = levels)
)

# A fast distance calculation with supplied bandwidths.
d <- dkss(df, bw = c(2, 0.25, 0.25))
d$distances[1:4, 1:4]

\dontrun{
set.seed(2)
# Select bandwidths using maximum-similarity cross-validation.
selected <- dkss(df, bw = "mscv")

# Select feature-specific kernels and bandwidths using MCMC.
fit <- mcmc.dkss(df, chains = 4, warmup = 1000, draws = 2000)

# Calculate distances using the sampled MAP state.
d_map <- dkss(fit)
}
}

\keyword{ ~multivariate }
\keyword{ ~nonparametric }
\keyword{ ~optimize }

\concept{ ~distances }
\concept{ ~metrics }
\concept{ ~metriclearning }
