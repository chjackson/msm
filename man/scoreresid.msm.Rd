\name{scoreresid.msm}
\alias{scoreresid.msm}
\title{Score residuals}
\description{
  Score residuals for detecting outlying subjects. 
}
\usage{
scoreresid.msm(x, plot=FALSE)
}
\arguments{
  \item{x}{A fitted multi-state model, as returned by
    \code{\link{msm}}.}

  \item{plot}{If \code{TRUE}, display a simple plot of the residuals in
    subject order, labelled by subject identifiers}
}
\value{
  Vector of the residuals, named by subject identifiers. 
}
\details{
  The score residual for a single subject is

  \deqn{U(\theta)^T I(\theta)^{-1} U(\theta)}{U(theta)^T I(theta)^{-1} U(theta)}
  
  where \eqn{U(\theta)}{U(theta)} is the vector of first derivatives of the
  log-likelihood for that subject at maximum likelihood estimates
  \eqn{\theta}{theta}, and \eqn{I(\theta)}{theta} is the observed Fisher
  information matrix, that is, the matrix of second derivatives of minus
  the log-likelihood for that subject at theta.

  Subjects with a higher influence on the maximum likelihood estimates
  will have higher score residuals.

  These are only available for models with analytic derivatives (which
  includes all non-hidden and most hidden Markov models).
}
\author{Andrew Titman  \email{a.titman@lancaster.ac.uk} (theory),   
  Chris Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk} (code)}
\keyword{models}