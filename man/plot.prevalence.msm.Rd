\name{plot.prevalence.msm}
\alias{plot.prevalence.msm}
\title{Plot of observed and expected prevalences}
\description{
  Provides a rough indication of goodness of fit of a multi-state model,
  by estimating the observed numbers of individuals occupying a state
  at a series of times, and plotting these against forecasts from the
  fitted model, for each state.  Observed prevalences are indicated as
  solid lines, expected prevalences as dashed lines. 
}
\usage{
\method{plot}{prevalence.msm}(x, mintime=NULL, maxtime=NULL, timezero=NULL,
                    initstates=NULL, interp=c("start","midpoint"),
                    censtime=Inf, subset=NULL,
                    covariates="population", misccovariates="mean",
                    piecewise.times=NULL, piecewise.covariates=NULL,
                    xlab="Times",ylab="Prevalence (\%)", lwd.obs=1,
                    lwd.exp=1, lty.obs=1, lty.exp=2, col.obs="blue",
                    col.exp="red", legend.pos=NULL,
                    ...)
}
\arguments{
  \item{x}{A fitted multi-state model produced by \code{\link{msm}}.}
  \item{mintime}{Minimum time at which to compute the observed and
    expected prevalences of states.}
  \item{maxtime}{Maximum time at which to compute the observed and
    expected prevalences of states.}
  \item{timezero}{Initial time of the Markov process. Expected values
    are forecasted from here. Defaults to the minimum of the observation
    times given in the data. }
  \item{initstates}{Optional vector of the same length as the number of
    states. Gives the numbers of individuals occupying each state at the
    initial time, to be used for forecasting expected prevalences.
    The default is those observed in the data.  These should add up to
    the actual number of people in the study at the start. }
  \item{interp}{Interpolation method for observed states, see
    \code{\link{prevalence.msm}}.}
  \item{censtime}{Subject-specific maximum follow-up times, see \code{\link{prevalence.msm}}.} 
  \item{subset}{Subset of subjects to calculated observed prevalences for.} 
  \item{covariates}{Covariate values for which to forecast expected
    state occupancy.  See \code{\link{prevalence.msm}} --- if this
    function runs too slowly, as it may if there are continuous covariates,
    replace \code{covariates="population"} with \code{covariates="mean"}. 
  }
  \item{misccovariates}{(Misclassification models only) Values of covariates on the misclassification
    probability matrix. See \code{\link{prevalence.msm}}. }  
  \item{piecewise.times}{Times at which piecewise-constant intensities
    change. See \code{\link{prevalence.msm}}.}
  \item{piecewise.covariates}{Covariates on which the piecewise-constant
    intensities depend. See \code{\link{prevalence.msm}}.}
  \item{xlab}{x axis label.}
  \item{ylab}{y axis label.}
  \item{lwd.obs}{Line width for observed prevalences. See \code{\link{par}}.}  
  \item{lwd.exp}{Line width for expected prevalences. See \code{\link{par}}.}
  \item{lty.obs}{Line type for observed prevalences. See \code{\link{par}}.} 
  \item{lty.exp}{Line type for expected prevalences. See \code{\link{par}}.} 
  \item{col.obs}{Line colour for observed prevalences. See \code{\link{par}}.}
  \item{col.exp}{Line colour for expected prevalences. See \code{\link{par}}.}
  \item{legend.pos}{Vector of the \eqn{x} and \eqn{y} position,
    respectively, of the legend.}
  \item{...}{Further arguments to be passed to the generic \code{\link{plot}}
    function.}
}
\details{
  See \code{\link{prevalence.msm}} for details of the assumptions
  underlying this method.

  Observed prevalences are plotted with a solid line, and expected
  prevalences with a dotted line.
}
\references{
  Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P.  Multi-state
  Markov models for analysing incomplete disease history data with
  illustrations for HIV disease.  \emph{Statistics in Medicine} (1994) 13(3):
  805--821.
}
\seealso{
  \code{\link{prevalence.msm}}
}
\keyword{models}
