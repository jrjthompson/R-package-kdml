\name{kdml}
\alias{kdml}
\alias{kdml-package}
\docType{package}
\title{Kernel Metric Learning for Mixed-type Data}
\description{
This package contains nonparametric kernel methods for calculating pairwise 
distances between mixed-type observations. These methods can be used in any 
distance based algorithm, with emphasis placed on usage in clustering or 
classification applications. Descriptions of the implementation of these methods
can be found in Ghashti (2024) and Ghashti and Thompson (2025).
}

\details{
The package contains DKPS and DKSS pairwise distances, their MSCV bandwidth
selectors, and a mixed-type data generator. A compiled MCMC implementation can
jointly sample feature-specific kernels and bandwidths, with CPU and optional
CUDA backends and posterior diagnostics.
}

\author{
  John R. J. Thompson \email{john.thompson@ubc.ca}, Jesse S. Ghashti
  \email{jesse.ghashti@ubc.ca}, and Daniel Krasnov
  \email{daniel.krasnov@mail.mcgill.ca}

  Maintainer: John R. J. Thompson \email{john.thompson@ubc.ca}

  We would like to acknowledge funding support from the University of British 
  Columbia Aspire Fund (\url{https://ok.ubc.ca/}). We also acknowledge support from
  the Natural Sciences and Engineering Research Council of Canada (NSERC). 

}
\references{

     Ghashti, J.S. (2024), \dQuote{Similarity Maximization and Shrinkage Approach 
     in Kernel Metric Learning for Clustering Mixed-type Data (T)}, University 
     of British Columbia. <https://dx.doi.org/10.14288/1.044397>

     Ghashti, J.S. and J.R.J Thompson (2024), \dQuote{Mixed-type Distance 
     Shrinkage and Selection for Clustering via Kernel Metric Learning}. 
     Journal of Classification, Accepted. 


}
\keyword{package}
