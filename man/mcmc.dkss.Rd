\name{mcmc.dkss}
\alias{mcmc.dkss}
\title{MCMC Kernel and Bandwidth Selection for DKSS}
\description{
Sample feature-specific kernels and bandwidths for distance using kernel
summation similarity (DKSS) with a generalized Bayesian target based on maximum
similarity cross-validation (MSCV).
}
\usage{
mcmc.dkss(
  df,
  chains = 1L,
  chain_threads = 1L,
  warmup = 1000L,
  draws = 1000L,
  thin = 1L,
  beta = 1,
  proposal_sd = 0.15,
  prior_mean = 0,
  prior_sd = 4,
  adapt = TRUE,
  target_accept = 0.44,
  kernel_candidates = NULL,
  initial_kernels = NULL,
  initial_bandwidths = NULL,
  backend = c("cpu", "cuda"),
  progress = interactive()
)
}
\arguments{
\item{df}{a \eqn{p}-variate data frame. The data types may be continuous
(\code{\link{numeric}}), nominal (\code{\link{factor}}), ordinal
(\code{\link{ordered}}), or any combination thereof. Columns of \code{df}
should be set to the appropriate data type class.}
\item{chains}{positive number of independent MCMC chains.}
\item{chain_threads}{positive number of CPU chains to run concurrently,
capped at \code{chains}. Values above one require an OpenMP-enabled build and
\code{backend = "cpu"}.}
\item{warmup}{non-negative number of adaptation draws discarded from each
chain.}
\item{draws}{positive number of retained draws per chain.}
\item{thin}{positive number of sampling sweeps between retained draws.}
\item{beta}{positive learning-rate multiplier for the MSCV score in
the generalized posterior target.}
\item{proposal_sd}{positive random-walk standard deviation on the
transformed bandwidth scale.}
\item{prior_mean, prior_sd}{mean and positive standard deviation of
the independent normal prior on each transformed bandwidth.}
\item{adapt}{logical; adapt proposal standard deviations during warmup.}
\item{target_accept}{target acceptance probability between zero and one,
exclusive, used by warmup adaptation.}
\item{kernel_candidates}{optional list of non-empty character vectors giving
the allowed kernels. A named list can replace candidates for selected
features; an unnamed list must contain one element per feature in data-frame
order. The default allows every type-compatible kernel.}
\item{initial_kernels}{optional character vector of initial kernels. A named
vector can initialize selected features; an unnamed vector must contain one
kernel per feature. Every initial kernel must occur in that feature's
candidate set.}
\item{initial_bandwidths}{optional positive bandwidths. If omitted,
type- and scale-aware values are chosen. Categorical values must also lie
below the upper bound for their initial kernel.}
\item{backend}{\code{"cpu"} for the portable implementation or
\code{"cuda"} to require an available CUDA backend. Defaults to
\code{"cpu"}.}
\item{progress}{logical; report sampler progress.}
}
\details{
Let \eqn{X=(\boldsymbol{x}_1,\ldots,\boldsymbol{x}_n)} contain \eqn{p}
features, let \eqn{m=(m_1,\ldots,m_p)} specify one kernel per feature, and let
\eqn{\boldsymbol{\lambda}=(\lambda_1,\ldots,\lambda_p)} contain the associated
bandwidths. For continuous, nominal, and ordinal feature sets
\eqn{\mathcal{J}_c}, \eqn{\mathcal{J}_u}, and \eqn{\mathcal{J}_o}, the DKSS
similarity is

\eqn{s_{\mathrm{DKSS}}(\boldsymbol{x}_i,\boldsymbol{x}_l\mid m,\boldsymbol{\lambda})
= \sum_{j\in\mathcal{J}_c}K_{m_j}(x_{ij},x_{lj};\lambda_j)
+ \sum_{j\in\mathcal{J}_u}L_{m_j}(x_{ij},x_{lj};\lambda_j)
+ \sum_{j\in\mathcal{J}_o}O_{m_j}(x_{ij},x_{lj};\lambda_j).}

The sampler uses the maximum-similarity cross-validation objective

\eqn{S(m,\boldsymbol{\lambda};X)=\frac{1}{n}\sum_{i=1}^n
\log\left\{\frac{1}{n-1}\sum_{\substack{l=1\\l\ne i}}^n
s_{\mathrm{DKSS}}(\boldsymbol{x}_i,\boldsymbol{x}_l\mid
m,\boldsymbol{\lambda})\right\}.}

Bandwidths are represented by unconstrained parameters
\eqn{\boldsymbol{\theta}}. Continuous bandwidths use a log transformation
relative to the sample standard deviation, while nominal and ordinal
bandwidths use a logit transformation. Writing this transformation as
\eqn{\boldsymbol{\lambda}=g(m,\boldsymbol{\theta})}, the generalized posterior
target is

\eqn{\pi(m,\boldsymbol{\theta}\mid X)\propto
\exp\left[\beta S\{m,g(m,\boldsymbol{\theta});X\}\right]
\prod_{j=1}^p\phi(\theta_j;\mu_j,\sigma_j^2)
\prod_{j=1}^p\frac{1}{|\mathcal{M}_j|},}

where \eqn{\mathcal{M}_j} is the candidate-kernel set for feature \eqn{j},
\eqn{\phi(\cdot;\mu_j,\sigma_j^2)} is its normal prior density, and
\eqn{\beta} controls the influence of MSCV relative to the prior.

During each sweep and for each feature, the sampler first proposes a different
kernel uniformly from \eqn{\mathcal{M}_j\setminus\{m_j\}} and accepts it with
probability \eqn{\min\{1,\exp(\ell^*-\ell)\}}, where \eqn{\ell} is the log of
the target above. It then proposes
\eqn{\theta_j^*=\theta_j+\delta_j\epsilon}, with
\eqn{\epsilon\sim N(0,1)}, and applies the same Metropolis acceptance rule.
When requested, \eqn{\delta_j} is adapted during warmup toward
\code{target_accept}. After warmup, every \code{thin}-th sweep is retained.

Scalar values for \code{proposal_sd}, \code{prior_mean}, \code{prior_sd}, and
\code{initial_bandwidths} are recycled over features. Non-scalar values must
either follow feature order or be named with every feature exactly once.

Default continuous candidates are \code{c_gaussian},
\code{c_epanechnikov}, \code{c_uniform}, \code{c_triangle},
\code{c_biweight}, \code{c_triweight}, \code{c_tricube}, \code{c_cosine},
\code{c_logistic}, \code{c_sigmoid}, and \code{c_silverman}. Nominal
candidates are \code{u_aitken} and \code{u_aitchisonaitken}; ordinal
candidates are \code{o_wangvanryzin}, \code{o_habbema}, \code{o_aitken},
\code{o_aitchisonaitken}, and \code{o_liracine}.

Use \code{set.seed()} before the call for reproducible chains. The retained MAP
state is the draw with the largest sampled log target.
}
\value{
An object of class \code{kdml_mcmc}, with the following components:
\item{theta, bandwidth, kernel, kernel_code}{three-dimensional arrays of
transformed bandwidths, natural-scale bandwidths, kernel names, and kernel
codes. The dimensions represent retained draws, features, and chains.}
\item{score, log_target}{matrices of the MSCV score and generalized log
target, with one row per retained draw and one column per chain.}
\item{draws}{a draws array that follows the conventions of the
\pkg{posterior} package.}
\item{map}{the retained draw and chain indices and all values at the sampled
MAP state.}
\item{diagnostics, invalid_proposals}{proposal scales, attempt and acceptance
counts, and invalid-proposal counts.}
\item{candidates, initial, control}{normalized candidate sets, initial state,
and sampler settings.}
\item{distance, backend}{the selected metric and backend actually used.}
\item{preprocessing}{feature types, scales, levels, names, and encoded fitting
data required for reconstruction.}
}
\references{
Bissiri, P. G., Holmes, C. C., and Walker, S. G. (2016). A general framework
for updating belief distributions. \emph{Journal of the Royal Statistical
Society: Series B}, 78(5), 1103--1130.
\href{https://doi.org/10.1111/rssb.12158}{doi:10.1111/rssb.12158}.

Ghashti, J. S. and Thompson, J. R. J. (2025). Mixed-type distance shrinkage
and selection for clustering via kernel metric learning. \emph{Journal of
Classification}, 42, 311--334.
\href{https://doi.org/10.1007/s00357-024-09493-z}{doi:10.1007/s00357-024-09493-z}.
}
\seealso{
\code{\link{mcmc.dkps}}, \code{\link{dkss}},
\code{\link{kdml_mcmc_methods}}, \code{\link{cuda_info}}
}
\examples{
set.seed(43)
dat <- data.frame(
  x = seq_len(12),
  group = factor(rep(c("a", "b"), 6)),
  rank = ordered(rep(c("low", "mid", "high"), 4),
                 levels = c("low", "mid", "high"))
)

fit <- mcmc.dkss(
  dat,
  chains = 2,
  warmup = 5,
  draws = 12,
  backend = "cpu",
  progress = FALSE
)
fit$map[c("draw", "chain", "kernel", "bandwidth")]
}
\author{
Daniel Krasnov \email{daniel.krasnov@mail.mcgill.ca}
}
