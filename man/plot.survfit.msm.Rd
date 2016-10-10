\name{plot.survfit.msm}
\alias{plot.survfit.msm}
\title{Plot empirical and fitted survival curves}
\description{
  Plot a Kaplan-Meier estimate of the survival probability and compare
  it with the fitted survival probability from a \code{msm} model. 
}
\usage{
\method{plot}{survfit.msm}(x, from=1, to=NULL, range=NULL, covariates="mean",
                 interp=c("start","midpoint"), ci=c("none","normal","bootstrap"), B=100,
                 legend.pos=NULL, xlab="Time", ylab="Survival probability",
                 lty=1, lwd=1, col="red", lty.ci=2, lwd.ci=1, col.ci="red", 
                 mark.time=TRUE, col.surv="blue", lty.surv=2, lwd.surv=1,
                 ...)
}
\arguments{

  \item{x}{Output from \code{\link{msm}}, representing a fitted
    multi-state model object.}

  \item{from}{State from which to consider survival. Defaults to
    state 1. }

  \item{to}{Absorbing state to consider. Defaults to the
    highest-labelled absorbing state. }

  \item{range}{Vector of two elements, giving the range of times to plot
  for.}

  \item{covariates}{Covariate values for which to evaluate the expected
    probabilities.  This can either be:\cr

    the string \code{"mean"}, denoting the means of the covariates in
    the data (this is the default),\cr

    the number \code{0}, indicating that all the covariates should be
    set to zero,\cr

    or a list of values, with optional names. For example

    \code{list (60, 1)}

    where the order of the list follows the order of the covariates
    originally given in the model formula, or a named list, 

    \code{list (age = 60, sex = 1)}
  }

  \item{ci}{If \code{"none"} (the default) no confidence intervals are
    plotted.  If \code{"normal"} or \code{"bootstrap"}, confidence
    intervals are plotted based on the respective method in
    \code{\link{pmatrix.msm}}. This is very computationally-intensive,
    since intervals must be computed at a series of times.}

  \item{B}{Number of bootstrap or normal replicates for the confidence
    interval.  The default is 100 rather than the usual 1000, since
    these plots are for rough diagnostic purposes.}
  
  \item{interp}{If \code{interp="start"} (the default) then the entry
    time into the absorbing state is assumed to be the time it is first
    observed in  the data.
  
  If \code{interp="midpoint"} then the entry time into the absorbing
  state is assumed to be halfway between the time it is first observed
  and the previous observation time. This is generally more reasonable
  for "progressive" models with observations at arbitrary times. }

  \item{legend.pos}{Vector of the \eqn{x} and \eqn{y} position,
    respectively, of the legend.}
  \item{xlab}{x axis label.}
  \item{ylab}{y axis label.}
  \item{lty}{Line type for the fitted curve. See \code{\link{par}}.}  
  \item{lwd}{Line width for the fitted curve. See \code{\link{par}}.}  
  \item{col}{Colour for the fitted curve. See \code{\link{par}}.}  
  \item{lty.ci}{Line type for the fitted curve confidence limits. See \code{\link{par}}.}  
  \item{lwd.ci}{Line width for the fitted curve confidence limits. See \code{\link{par}}.}  
  \item{col.ci}{Colour for the fitted curve confidence limits. See \code{\link{par}}.}  
  \item{mark.time}{Mark the empirical survival curve at each censoring
  point, see \code{\link[survival]{lines.survfit}}.}
  \item{col.surv}{Colour for the empirical survival curve, passed to \code{\link[survival]{lines.survfit}}. See \code{\link{par}}.}
  \item{lty.surv}{Line type for the empirical survival curve, passed to \code{\link[survival]{lines.survfit}}. See \code{\link{par}}.}  
  \item{lwd.surv}{Line width for the empirical survival curve, passed to \code{\link[survival]{lines.survfit}}. See \code{\link{par}}.}  
  \item{...}{Other arguments to be passed to the
    \code{\link{plot}} function which draws the fitted curve, or the \code{\link[survival]{lines.survfit}}
    function which draws the empirical curve. }
}
\details{
  If the data represent observations of the process at arbitrary times, then
  the first occurrence of the absorbing state in the data will usually
  be greater than the actual first transition time to that state.
  Therefore the Kaplan-Meier estimate of the survival probability will
  be an overestimate.

  The method of Turnbull (1976) could be used to give a non-parametric
  estimate of the time to an interval-censored event, and compared to
  the equivalent estimate from a multi-state model.  This is implemented
  in the CRAN package \pkg{interval} (Fay and Shaw 2010).
  
  This currently only handles time-homogeneous models. 
}
\references{
  Turnbull, B. W. (1976) The empirical distribution function with arbitrarily grouped, censored and
  truncated data. J. R. Statist. Soc. B 38, 290-295.

  Fay, MP and Shaw, PA (2010). Exact and Asymptotic Weighted Logrank Tests for Interval Censored
  Data: The interval R package. Journal of Statistical Software. http://www.jstatsoft.org/v36/
  i02/. 36 (2):1-34.
}
\seealso{\code{\link[survival]{survfit}}, \code{\link[survival]{plot.survfit}}, \code{\link{plot.prevalence.msm}}}
\keyword{models}
