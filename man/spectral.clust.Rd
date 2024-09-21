\name{spectral.clust}
\alias{spectral.clust}
\title{
Spectral Clustering using Similarity or Distance Matrices
}
\description{
This function calculates performs spectral clustering with the k-means step 
using precomputed similarity or distance matrices, and returns a vector of 
cluster assignments.
}
\usage{
spectral.clust(S, k, nstart = 10, iter.max = 1000, 
              is.sim = NULL, neighbours = 10)
}

\arguments{
  \item{S}{
    a \eqn{n \times n} numeric matrix representing either pairwise similarities
    or distances between observations. The matrix can be a similarity matrix or 
    a distance matrix, as indicated by the \code{is.sim} argument.
  }

  \item{k}{
    integer value specifying the number of clusters to form. This is passed to 
    the \code{kmeans} algorithm.
  }

  \item{nstart}{
    integer value specifying the number of random starts for the bandwidth 
    estimation. Defaults to 3 or the number of variables, whichever is larger.
  }

  \item{iter.max}{
    integer value specifying the maximum number of iterations for the 
    \code{kmeans} algorithm. Defaults to \code{1000}.
  }

  \item{is.sim}{
    logical value indicating whether the input matrix \code{S} is a similarity 
    matrix. If set to \code{TRUE}, \code{S} is treated as a similarity matrix. 
    If set to \code{FALSE}, \code{S} is treated as a distance matrix. Must be 
    specified.
  }

  \item{neighbours}{
    integer value specifying the number of nearest neighbours to consider when 
    constructing the graph Laplacian. This helps in determining the structure 
    of the graph from the similarity or distance matrix. Defaults to \code{10}.
  }
}

\details{
  \code{spectral.clust} implements spectral clustering on pairwise similarity or
  distance matrices, following the method described by Ng et al. (2001). The 
  function first constructs an adjacency matrix from the input similarity or 
  distance matrix \code{S} using the \code{neighbours} parameter to define the 
  nearest connections. If \code{S} is a similarity matrix (\code{is.sim = TRUE}),
  the function retains the largest values corresponding to the \code{neighbours}
  nearest observations. If \code{S} is a distance matrix (\code{is.sim = FALSE}), 
  it retains the smallest values for the nearest observations. The adjacency 
  matrix is symmetrized and used to compute the unnormalized Laplacian matrix. 
  The eigenvectors corresponding to the smallest eigenvalues of the Laplacian 
  are extracted and clustered using the \code{kmeans} algorithm. The number of 
  clusters, \code{k}, and parameters such as the number of random starts 
  (\code{nstart}) and maximum iterations (\code{iter.max}) for the \code{kmeans}
  step are user-specified. 
}


\value{
  \code{spectral.clust} returns a \code{list} object with the following components:
  
  \item{clusters}{
    an \eqn{n}-variate integer vector indicating the cluster assignment for each 
    observation, as determined by the \code{kmeans} algorithm.
  }

  \item{S}{
    the original \eqn{n \times n} numeric matrix used as input, representing 
    either pairwise similarities or distances between observations, depending 
    on the \code{is.sim} argument.
  }
}

\references{

  Ng, A., Jordan, M., & Weiss, Y. (2001). On spectral clustering: Analysis and 
  an algorithm. \dQuote{Advances in Neural Information processing systems}, 14.
  
}
\author{
    John R. J. Thompson \email{john.thompson@ubc.ca}, Jesse S. Ghashti
    \email{jesse.ghashti@ubc.ca}
}


\seealso{
 \code{\link{mscv.dkps}}, \code{\link{dkps}}, \code{\link{mscv.dkss}}, 
 \code{\link{dkss}}, \code{link{kss}}
}
\examples{
# load the Iris dataset
dat <- iris[,-5]

# calculate pairwise similarities using maximum likelihood cross validation
S <- kss(dat, bw = "np", npmethod = "cv.ml", cFUN = "c_gaussian", verbose = TRUE)

# cluster points using spectral clustering and compare to true class labels
cl <- spectral.clust(S$similarities, 3, is.sim = TRUE)
table(cl$clusters, iris[,5])

# try a different number of neighbours
cl2 <- spectral.clust(S$similarities, 3, is.sim = TRUE, neighbours = 4)
table(cl2$clusters, iris[,5])
}

\keyword{ ~multivariate }
\keyword{ ~nonparametric }
\keyword{ ~optimize }

\concept{ ~similaritylearning}
\concept{ ~similarity }
\concept{ ~clustering }
\concept{ ~spectralclustering }
