
\name{prevalence.msm}
\title{Tables of observed and expected prevalences}
\alias{prevalence.msm}
\concept{Goodness of fit}
\description{
  This provides a rough indication of the goodness of fit of a
  multi-state model, by estimating the observed numbers of individuals
  occupying each state at a series of times, and comparing these with
  forecasts from the fitted model.
}
\usage{
prevalence.msm(x, times=NULL, timezero=NULL, initstates=NULL, covariates="population",
               misccovariates="mean", piecewise.times=NULL, piecewise.covariates=NULL,
	       ci=c("none","normal","bootstrap"), cl=0.95, B=1000, cores=NULL,
               interp=c("start","midpoint"), censtime=Inf, subset=NULL, plot=FALSE, ...)
}
\arguments{
  \item{x}{A fitted multi-state model produced by \code{\link{msm}}.}
  \item{times}{Series of times at which to compute the observed and
    expected prevalences of states.}
  \item{timezero}{Initial time of the Markov process. Expected values
    are forecasted from here. Defaults to the minimum of the observation
    times given in the data. }
  \item{initstates}{Optional vector of the same length as the number of
    states. Gives the numbers of individuals occupying each state at the
    initial time, to be used for forecasting expected prevalences.
    The default is those observed in the data.  These should add up to
    the actual number of people in the study at the start. }
  \item{covariates}{Covariate values for which to forecast expected
    state occupancy.  With the default \code{covariates="population"},
    expected prevalences are produced by summing model predictions over
    the covariates observed in the original data, for a fair comparison
    with the observed prevalences.  This may be slow, particularly with
    continuous covariates.

    Predictions for fixed covariates can be obtained by supplying
    covariate values in the standard way, as in
    \code{\link{qmatrix.msm}}. Therefore if
    \code{covariates="population"} is too slow, using the mean observed
    values through \code{covariates="mean"} may give a reasonable
    approximation.
  }
  \item{misccovariates}{(Misclassification models only)
    Values of covariates on the misclassification probability matrix
    for converting expected true to expected misclassified states.
    Ignored if \code{covariates="population"}, otherwise
    defaults to the mean values of the covariates in the data set.}
  \item{piecewise.times}{Times at which piecewise-constant intensities
    change. See \code{\link{pmatrix.piecewise.msm}} for how
    to specify this.  Ignored if \code{covariates="population"}.
    This is only required for time-inhomogeneous
    models specified using explicit time-dependent covariates, and
    should not be used for models specified using "pci". }
  \item{piecewise.covariates}{Covariates on which the piecewise-constant
    intensities depend. See \code{\link{pmatrix.piecewise.msm}} for how
    to specify this. Ignored if \code{covariates="population"}.}

  \item{ci}{If \code{"normal"}, then calculate a confidence interval for
    the expected prevalences by simulating \code{B} random vectors
    from the asymptotic multivariate normal distribution implied by the
    maximum likelihood estimates (and covariance matrix) of the log
    transition intensities and covariate effects, then calculating the
    expected prevalences for each replicate.

    If \code{"bootstrap"} then calculate a confidence interval by
    non-parametric bootstrap refitting.  This is 1-2 orders of magnitude
    slower than the \code{"normal"} method, but is expected to be more
    accurate. See \code{\link{boot.msm}} for more details of
    bootstrapping in \pkg{msm}.

    If \code{"none"} (the default) then no confidence interval is
    calculated.}

  \item{cl}{Width of the symmetric confidence interval, relative to 1}
  \item{B}{Number of bootstrap replicates}

  \item{cores}{Number of cores to use for bootstrapping using parallel
    processing. See \code{\link{boot.msm}} for more details.}

  \item{interp}{Suppose an individual was observed in states
  \eqn{S_{r-1}} and \eqn{S_r} at two consecutive times \eqn{t_{r-1}}
  and \eqn{t_r}, and we want to estimate 'observed' prevalences at a
  time \eqn{t} between \eqn{t_{r-1}} and \eqn{t_r}.

  If \code{interp="start"}, then individuals are assumed to be in
  state \eqn{S_{r-1}} at time \eqn{t}, the same state as they were at \eqn{t_{r-1}}.

  If \code{interp="midpoint"} then if \eqn{t <= (t_{r-1} + t_r) / 2},
     the midpoint of \eqn{t_{r-1}} and \eqn{t_r}, the state at \eqn{t}
     is assumed to be \eqn{S_{r-1}}, otherwise \eqn{S_{r}}. This is
     generally more reasonable for "progressive" models.
   }

   \item{censtime}{If the time is greater than \code{censtime} and the
     patient has reached an absorbing state, then that subject will be
     removed from the risk set.  For example, if patients have died but
     would only have been observed up to this time, then this avoids
     overestimating the proportion of people who are dead at later
     times.

     This can be supplied as a single value, or as a vector with one
     element per subject (after any \code{subset} has been taken), in
     the same order as the original data.  This vector also only
     includes subjects with complete data, thus it excludes for example
     subjects with only one observation (thus no observed transitions),
     and subjects for whom every observation has missing values.
     (Note, to help construct this, the complete data used for the model
     fit can be accessed with \code{model.frame(x)}, where \code{x} is the
     fitted model object)

     This is ignored if it is less than the subject's maximum
     observation time.  }
   \item{subset}{Subset of subjects to calculate observed prevalences for.}

   \item{plot}{Generate a plot of observed against expected
     prevalences. See \code{\link{plot.prevalence.msm}}}

   \item{...}{Further arguments to pass to \code{\link{plot.prevalence.msm}}.}
}
\value{
  A list of matrices, with components:

  \item{Observed}{Table of observed numbers of individuals in each state at
    each time}

  \item{Observed percentages}{Corresponding percentage of the
    individuals at risk at each time.}

  \item{Expected}{Table of corresponding expected numbers.}

  \item{Expected percentages}{Corresponding percentage of the
    individuals at risk at each time.}

  Or if \code{ci.boot = TRUE}, the component \code{Expected} is a list
  with components \code{estimates} and \code{ci}.\cr
  \code{estimates} is a matrix of the expected prevalences, and
\code{ci} is a list of two matrices, containing the confidence limits.
The component
  \code{Expected percentages} has a similar format.

}

\details{
  The fitted transition probability matrix is used to forecast expected prevalences from the
  state occupancy at the initial time.  To
  produce the expected number in state \eqn{j} at time \eqn{t} after the
  start, the number of individuals under observation at time \eqn{t}
  (including those who have died, but not those lost to follow-up) is
  multiplied by the product of the proportion of individuals in each
  state at the initial time and the transition probability matrix in the
  time interval \eqn{t}.   The proportion of individuals in each state
  at the "initial" time is estimated, if necessary, in the same way as
  the observed prevalences.

  For misclassification models (fitted using an \code{ematrix}), this
  aims to assess the fit of the full model for the \emph{observed}
  states.  That is, the combined Markov progression model for the true
  states and the misclassification model. Thus, expected prevalences of \emph{true}
  states are estimated from the assumed proportion
  occupying each state at the initial time using the fitted transition
  probabiliy matrix. The vector of expected prevalences of true states
  is then multiplied by the fitted misclassification probability matrix
  to obtain the expected prevalences of \emph{observed} states.

  For general hidden Markov models, the observed state is taken to be the
  predicted underlying state from the Viterbi algorithm
  (\code{\link{viterbi.msm}}).  The goodness of fit of
  these states to the underlying Markov model is tested.

  For an example of this approach, see Gentleman \emph{et al.} (1994).
}

\references{
  Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P.  Multi-state
  Markov models for analysing incomplete disease history data with
  illustrations for HIV disease.  \emph{Statistics in Medicine} (1994) 13(3):
  805--821.

  Titman, A.C., Sharples, L. D.  Model diagnostics for multi-state
  models. \emph{Statistical Methods in Medical Research} (2010)
  19(6):621-651.
}

\seealso{
  \code{\link{msm}}, \code{\link{summary.msm}}
}

\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
