\name{plotprog.msm}
\alias{plotprog.msm}
\title{Kaplan Meier estimates of incidence}
\description{
  Compute and plot Kaplan-Meier estimates of the probability that each
  successive state has not occurred yet. 
}
\usage{
plotprog.msm(formula, subject, data, legend.pos=NULL, xlab="Time",
                   ylab="1 - incidence probability", lwd=1, xlim=NULL,
                   mark.time=TRUE, ...)
}
\arguments{
  \item{formula}{A formula giving the vectors containing the observed
    states and the  corresponding observation times. For example,
    
    \code{state ~ time}

    Observed states should be in the set \code{1, \dots, n},
    where \code{n} is the number of states.}
  \item{subject}{Vector of subject identification numbers for the data
    specified by \code{formula}. If missing, then all observations
    are assumed to be on the same subject. These must be sorted so that
    all observations on the same subject are adjacent.}
  \item{data}{An optional data frame in which the variables represented
    by  \code{state}, \code{time} and \code{subject} can be found.}
  \item{legend.pos}{Vector of the \eqn{x} and \eqn{y} position,
    respectively, of the legend.}

  \item{xlab}{x axis label.}
  \item{ylab}{y axis label.}
  \item{lwd}{Line width. See \code{\link{par}}.}
  \item{xlim}{x axis limits, e.g. c(0,10) for an axis ranging from 0 to 10. Default is the range of observation times.}
  \item{mark.time}{Mark the empirical survival curve at each censoring
  point, see \code{\link{lines.survfit}}.}
  
  \item{...}{Other arguments to be passed to the
    \code{\link{plot}} and \code{\link[survival]{lines.survfit}} functions.}
}
\details{
  If the data represent observations of the process at arbitrary times, then
  the first occurrence of the state in the data will usually be greater
  than the actual first transition time to that state.  Therefore the
  probabilities plotted by \code{\link{plotprog.msm}} will be
  overestimates. 
}
\seealso{\code{\link[survival]{survfit}}, \code{\link[survival]{plot.survfit}}} 
\keyword{models}
