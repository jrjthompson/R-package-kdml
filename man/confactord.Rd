\name{confactord}
\alias{confactord}
\title{
Mixed-type Data Generation with True Membership Labels
}
\description{
This function generates a mixed-type data frame with a combination of numeric (continuous), factor (nominal), and ordered (ordinal) variables with prespecified cluster overlap for each variable type. \code{confactord} is flexible, as the user may specify the number of each variable type, the amount of variables per variable type that have cluster overlap, the amount of cluster overlap for each variable type, the number of levels for the nominal and ordinal variables, and proportion of observations per class membership. Within and across-type variables are generated independently from one another. Currently, only
}
\usage{
confactord(n = 200, 
            popProb = c(0.5,0.5), 
            numMixVar = c(1,1,1), 
            numMixVarOl = c(1,1,1),  
            olVarType = c(0.1,0.1,0.1), 
            catLevels = c(2,4))
}
\arguments{
  \item{n}{
    integer number of observations to be generated. Defaults to \code{n = 200}
}
  \item{popProb}{
   numeric vector of length two specifying the proportion of observations allocated to each class membership, which must sum to one. Defaults to \code{popProb = c(0.5, 0.5)}
}
  \item{numMixVar}{
   numeric vector of integers of length three specifying (in order) the total number of continuous (numeric), nominal (factor), and ordinal (ordered) variables to be generated. If a specific variable type is not required, set the appropriate vector indice to zero. Defaults to \code{numMixVar = c(1,1,1)}
}
 \item{numMixVarOl}{
    numeric vector of integers of length three specifying (in order) the total number of continuous (numeric), nominal (factor), and ordinal (ordered) variables that will have class membership overlap. If all variables are to be well-separated by class membership, set all indices to zero. No indice of this vector may be greater than the corresponding indice in \code{numMixVar}. Defaults to \code{numMixVarOl = c(1,1,1)}
}
 \item{olVarType}{
   numeric vector of length three specifying (in order) the percentage of class membership overlap to be applied to the continuous (numeric), nominal (factor), and ordinal (ordered) No argument required if \code{numMixVarOl = c(0,0,0)}. Permissible class membership overlap per variable type is between 0.01 and 0.99. Defaults to ten percent overlap per variable type, \code{olVarType = c(0.1,0.1,0.1)}
}
 \item{catLevels}{
  numeric vector of length two specifying (in order) the number of levels (integer values) for each of the nominal (factor) and ordinal (ordered) variable types. Defaults to \code{catLevels = c(2,4)}
}
}
\details{
  Continuous variables are generated independently from normal distributions, with means determined by true class membership. If overlap is specified, additional variance is introduced to simulate cluster overlap. Nominal variables are generated using Dirichlet distributions representing different population proportions. Ordinal variables are initially simulated as continuous variables and then discretized into ordered categories based on quantile distributions, similar to a latent class model where ordinal categories are inferred based on underlying continuous distributions and adjusted for cluster overlap parameters.
  }

\value{
  \code{confactord} returns a \code{list} object, with the
  following components:
 \item{data}{a \code{\link{data.frame}} of mixed variable types based on user-specified parameters}
 \item{class}{a numeric vector of integers specifying the true class memberships for the returned \code{data} data frame}
}

\author{
    John R. J. Thompson \email{john.thompson@ubc.ca}, Jesse S. Ghashti
    \email{jesse.ghashti@ubc.ca}
}


\seealso{
 \code{\link{mscv.dkss}}, \code{\link{mscv.dkps}}, \code{\link{dkss}}, \code{\link{dkps}}
}

\examples{
\dontrun{
# EXAMPLE1: Default implementation generates the following
# 200 observations split into two clusters of equal size (100 observations each) 
# Three variables-- one of each numeric, factor, and ordered
# Each variable has ten percent cluster overlap
# Nominal variable is binary
# Ordinal variable has four levels

df1 <- confactord()


# EXAMPLE2: 
# 500 observations; 100 observations in cluster one and 400 in cluster two 
# Three continuous variables, two nominal, one ordinal
# Only one continuous variable has cluster overlap
# All nominal and ordinal variables have cluster overlap
# Cluster overlap for continuous variable is twenty percent
# Cluster overlap for nominal variables are thirty percent
# Cluster overlap for ordinal variable is fourty percent
# Nominal variable has three levels, while ordinal has 5

df2 <- confactord(n = 500, 
                    popProb = c(0.2,0.8), 
                    numMixVar = c(3,2,1), 
                    numMixVarOl = c(1,2,1),  
                    olVarType = c(0.2,0.3,0.4), 
                    catLevels = c(3,5))
}
}


