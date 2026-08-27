\name{kdml_mcmc_methods}
\alias{kdml_mcmc_methods}
\alias{print.kdml_mcmc}
\alias{summary.kdml_mcmc}
\alias{print.summary.kdml_mcmc}
\alias{plot.kdml_mcmc}
\alias{as_draws.kdml_mcmc}
\alias{as_draws_array.kdml_mcmc}
\alias{kdml.diagnostics}
\title{Work with KDML MCMC Fits}
\description{
Print or summarize a KDML MCMC fit, create diagnostic and posterior plots,
convert retained samples to a \pkg{posterior} draws object, or calculate
standard Monte Carlo diagnostics.
}
\usage{
\method{print}{kdml_mcmc}(x, ...)

\method{summary}{kdml_mcmc}(
  object,
  rhat_threshold = 1.01,
  min_ess_per_chain = 100,
  ...
)

\method{print}{summary.kdml_mcmc}(
  x,
  digits = max(3L, getOption("digits") - 3L),
  ...
)

\method{plot}{kdml_mcmc}(
  x,
  type = c("trace", "acf", "diagnostics", "score", "bandwidth",
           "kernel", "acceptance"),
  features = NULL,
  chains = NULL,
  lag.max = NULL,
  ...
)

\method{as_draws}{kdml_mcmc}(
  x,
  variables = c("theta", "bandwidth", "score", "log_target"),
  ...
)

\method{as_draws_array}{kdml_mcmc}(
  x,
  variables = c("theta", "bandwidth", "score", "log_target"),
  ...
)

kdml.diagnostics(object, variables = "theta")
}
\arguments{
\item{x, object}{a \code{kdml_mcmc} fit, except that \code{x} for the
summary print method is a \code{summary.kdml_mcmc} object.}
\item{variables}{one or more retained variable groups. Choices are
\code{"theta"}, \code{"bandwidth"}, \code{"kernel_code"},
\code{"score"}, and \code{"log_target"}. The diagnostics helper defaults
to \code{"theta"}.}
\item{rhat_threshold}{finite number greater than or equal to one. Features
with larger rank-normalized split \eqn{\widehat R} values are reported in the
convergence warnings.}
\item{min_ess_per_chain}{positive target for both bulk and tail effective
sample size per chain. The summary compares each total ESS with this value
multiplied by the number of chains.}
\item{digits}{number of significant digits used by the summary print method.}
\item{type}{plot type. The default is \code{"trace"}. Available types are
transformed-bandwidth traces and autocorrelations, \eqn{\widehat R}/ESS
diagnostics, score traces, bandwidth densities, kernel occupancy, and
post-warmup acceptance rates.}
\item{features}{feature names or one-based feature indices to plot. The
default selects every feature. This selection is ignored by
\code{type = "score"}.}
\item{chains}{chain names (for example, \code{"chain1"}) or one-based chain
indices. The default selects every chain.}
\item{lag.max}{non-negative maximum lag for autocorrelation plots. The default
is the smaller of 20 and one less than the retained draw count.}
\item{...}{for \code{plot()}, named graphical parameters passed to the base
graphics function that draws each panel. Reserved for method compatibility by
the other methods.}
}
\details{
The summary reports the retained MSCV score and log-target distributions, the
retained MAP state, diagnostics for the unconstrained bandwidth parameters,
per-chain kernel occupancy and transitions, and proposal acceptance. A warning
that only one chain was run is expected when convergence cannot be assessed
across chains.

Continuous bandwidth density plots use the log-bandwidth scale. Categorical
bandwidths are plotted on their natural bounded scale. Kernel plots summarize
the sampled character-valued kernels, while draws objects represent kernels by
their numeric \code{kernel_code} values.

Feature-indexed draws have names such as \code{theta[feature]} and
\code{bandwidth[feature]}. \code{kdml.diagnostics()} delegates its calculations
to \pkg{posterior}. Rank-normalized split \eqn{\widehat R} should be interpreted
using multiple chains; effective sample sizes may be unavailable for very short
or constant chains.
}
\value{
\code{as_draws()} and \code{as_draws_array()} return a
\code{\link[posterior:draws_array]{draws_array}} with iteration, chain, and
variable dimensions.

\code{kdml.diagnostics()} returns a data frame with one row per selected
variable and columns \code{variable}, \code{mean}, \code{sd}, \code{rhat},
\code{ess_bulk}, \code{ess_tail}, and \code{mcse_mean}.

\code{summary()} returns an object of class \code{summary.kdml_mcmc}, a list
with components \code{call}, \code{metadata}, \code{score}, \code{map},
\code{theta}, \code{convergence}, \code{kernel_occupancy},
\code{kernel_transitions}, and \code{acceptance}.

The fitted-object print method invisibly returns \code{x}. The summary print
method invisibly returns its summary object. The plot method invisibly returns
the fitted \code{kdml_mcmc} object.
}
\seealso{
\code{\link{mcmc.dkps}}, \code{\link{mcmc.dkss}},
\code{\link{dkps}}, \code{\link{dkss}},
\code{\link[posterior]{summarise_draws}}
}
\examples{
set.seed(31)
dat <- data.frame(
  x = seq_len(12),
  group = factor(rep(c("a", "b"), 6))
)
fit <- mcmc.dkss(
  dat, chains = 2, warmup = 5, draws = 12,
  backend = "cpu", progress = FALSE
)

fit_summary <- suppressWarnings(summary(fit, min_ess_per_chain = 2))
fit_summary$score
posterior::as_draws(fit, variables = "theta")
posterior::as_draws_array(fit, variables = "bandwidth")
suppressWarnings(kdml.diagnostics(fit))
plot(fit, type = "trace", features = "x")
}
\author{
Daniel Krasnov \email{daniel.krasnov@mail.mcgill.ca}
}
