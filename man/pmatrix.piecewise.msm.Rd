\name{pmatrix.piecewise.msm}
\alias{pmatrix.piecewise.msm}
\title{Transition probability matrix for processes with
 piecewise-constant intensities}
\description{
  Extract the estimated transition probability matrix from a fitted
  non-time-homogeneous multi-state model for a given time interval.
  This is a generalisation of \code{\link{pmatrix.msm}} to
  models with time-dependent covariates.  Note that
  \code{\link{pmatrix.msm}} is sufficient to calculate transition
  probabilities for time-inhomogeneous
  models fitted using the \code{pci} argument to \code{\link{msm}}.
}
\usage{
pmatrix.piecewise.msm(x=NULL, t1, t2, times, covariates,
ci=c("none","normal","bootstrap"), cl=0.95, B=1000, cores=NULL,
qlist=NULL,...)
}
\arguments{

  \item{x}{A fitted multi-state model, as returned by
    \code{\link{msm}}. This should be a non-homogeneous model, whose
    transition intensity matrix depends on a time-dependent covariate. }

  \item{t1}{The start of the time interval to estimate the transition probabilities for.}

  \item{t2}{The end of the time interval to estimate the transition probabilities
    for.}

  \item{times}{Cut points at which the transition intensity matrix changes.}

  \item{covariates}{
    A list with number of components one greater than the length of
    \code{times}.  Each component of the list is specified in the same
    way as the \code{covariates} argument to \code{\link{pmatrix.msm}}.
    The components correspond to the covariate values in the intervals

    \code{(t1, times[1]],  (times[1], times[2]], ...,
      (times[length(times)], t2]}

    (assuming that all elements of \code{times} are in the interval
    \code{(t1, t2)}).
  }
  \item{ci}{If \code{"normal"}, then calculate a confidence interval for
    the transition probabilities by simulating \code{B} random vectors
    from the asymptotic multivariate normal distribution implied by the
    maximum likelihood estimates (and covariance matrix) of the log
    transition intensities and covariate effects, then calculating the
    resulting transition probability matrix for each replicate.

    If \code{"bootstrap"} then calculate a confidence interval by
    non-parametric bootstrap refitting.  This is 1-2 orders of magnitude
    slower than the \code{"normal"} method, but is expected to be more
    accurate. See \code{\link{boot.msm}} for more details of
    bootstrapping in \pkg{msm}.

    If \code{"none"} (the default) then no confidence interval is
    calculated.}

  \item{cl}{Width of the symmetric confidence interval, relative to 1.}

  \item{B}{Number of bootstrap replicates, or number of normal
    simulations from the distribution of the MLEs}

  \item{cores}{Number of cores to use for bootstrapping using parallel
    processing. See \code{\link{boot.msm}} for more details.}

  \item{qlist}{A list of transition intensity matrices, of length one
  greater than the length of \code{times}.    Either this or
  a fitted model \code{x} must be supplied.   No confidence intervals
  are available if (just) \code{qlist} is supplied.}

  \item{...}{Optional arguments to be passed to \code{\link{MatrixExp}} to
    control the method of computing the matrix exponential.}
}
\details{
    Suppose a multi-state model has been fitted, in which
  the transition intensity matrix \eqn{Q(x(t))} is modelled in terms of
  time-dependent covariates \eqn{x(t)}.  The transition probability
  matrix \eqn{P(t_1, t_n)}{P(t1, tn)} for the time interval \eqn{(t_1,
  t_n)}{(t1, tn)} cannot be calculated from the estimated intensity matrix as
  \eqn{\exp((t_n - t_1) Q)}{exp((tn - t1) Q)}, because \eqn{Q} varies within the interval
  \eqn{t_1, t_n}{t1, tn}.  However, if the covariates are piecewise-constant,
  or can be approximated as piecewise-constant, then we can calculate
  \eqn{P(t_1, t_n)}{P(t1, tn)} by multiplying together individual matrices
  \eqn{P(t_i,
  t_{i+1}) = \exp((t_{i+1} - t_i) Q)}{P(t_i,
  t_{i+1}) = exp((t_{i+1} - t_i) Q)}, calculated over intervals where Q is
constant:

\deqn{P(t_1, t_n) =  P(t_1, t_2) P(t_2, t_3)\ldots P(t_{n-1},
  t_n)}{P(t1, tn) =  P(t1, t2) P(t2, t3)\ldotsP(tn-1, tn)}
}
\value{
  The matrix of estimated transition probabilities \eqn{P(t)} for the
  time interval \code{[t1, tn]}.  That is, the probabilities of
  occupying state \eqn{s} at time \eqn{t_n}{tn}
  conditionally on occupying state \eqn{r} at time \eqn{t_1}{t1}.
  Rows correspond to "from-state" and columns to "to-state".
}
\seealso{
  \code{\link{pmatrix.msm}} }

\examples{
\dontrun{
## In a clinical study, suppose patients are given a placebo in the
## first 5 weeks, then they begin treatment 1 at 5 weeks, and
## a combination of treatments 1 and 2 from 10 weeks.
## Suppose a multi-state model x has been fitted for the patients'
## progress, with treat1 and treat2 as time dependent covariates.

## Cut points for when treatment covariate changes
times <- c(0, 5, 10)

## Indicators for which treatments are active in the four intervals
## defined by the three cut points
covariates <- list( list (treat1=0, treat2=0), list (treat1=0, treat2=0), list(treat1=1, treat2=0),
list(treat1=1, treat2=1) )

## Calculate transition probabilities from the start of the study to 15 weeks
pmatrix.piecewise.msm(x, 0, 15, times, covariates)
}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
