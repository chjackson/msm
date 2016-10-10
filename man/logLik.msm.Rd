\name{logLik.msm}
\alias{logLik.msm}
\title{Extract model log-likelihood}
\description{
  Extract the log-likelihood and the number of parameters of a model
  fitted with \code{\link{msm}}.  
}
\usage{
\method{logLik}{msm}(object, by.subject=FALSE, \dots)
}
\arguments{
  \item{object}{A fitted multi-state model object, as returned by
    \code{\link{msm}}.}
  \item{by.subject}{Return vector of subject-specific log-likelihoods,
    which should sum to the total log-likelihood.}
  \item{...}{(unused) further arguments passed to or from other methods.}
}
\value{
  The log-likelihood of the model represented by
  'object' evaluated at the maximum likelihood estimates.

  Akaike's information criterion can also be computed using
  \code{\link{AIC}(object)}.
}
\seealso{
  \code{\link{msm}},\code{\link{lrtest.msm}}.
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
