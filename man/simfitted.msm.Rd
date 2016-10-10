\name{simfitted.msm}
\alias{simfitted.msm}
\title{Simulate from a Markov model fitted using msm}
\concept{Simulation}
\description{
  Simulate a dataset from a Markov model fitted using \code{\link{msm}}, using the
  maximum likelihood estimates as parameters, and the
  same observation times as in the original data.
}
\usage{
simfitted.msm(x, drop.absorb=TRUE, drop.pci.imp=TRUE)
}
\arguments{
  \item{x}{A fitted multi-state model object as returned by
    \code{\link{msm}}.}
  \item{drop.absorb}{Should repeated observations in an absorbing state
    be omitted.  Use the default of \code{TRUE} to avoid warnings when using the
    simulated dataset for further \code{\link{msm}} fits.  Or set to
    \code{FALSE} if exactly the same number of observations as the original
    data are needed.}
  \item{drop.pci.imp}{In time-inhomogeneous models fitted using the
    \code{pci} option to \code{\link{msm}}, censored observations
    are inserted into the data by \code{\link{msm}} at the times
    where the intensity changes, but dropped by default when simulating
    from the fitted model using this function. Set this argument to
    \code{FALSE} to keep these observations and the corresponding
    indicator variable.}
}
\value{
  A dataset with variables as described in \code{\link{simmulti.msm}}.
}
\details{
  This function is a wrapper around \code{\link{simmulti.msm}},
  and only simulates panel-observed data.  To generate datasets with
  the exact times of transition, use the lower-level
  \code{\link{sim.msm}}.

  Markov models with misclassified states fitted through the
  \code{ematrix} option to \code{\link{msm}} are supported, but not 
  general hidden Markov models with \code{hmodel}.  For
  misclassification models, this function includes misclassification in
  the simulated states. 
  
  This function is used for parametric bootstrapping to estimate the
  null distribution of the test statistic in \code{\link{pearson.msm}}.
}
\seealso{
  \code{\link{simmulti.msm}}, \code{\link{sim.msm}}, \code{\link{pearson.msm}}, \code{\link{msm}}.
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
