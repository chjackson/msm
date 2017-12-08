\name{updatepars.msm}
\alias{updatepars.msm}
\title{updatepars.msm}
\description{
  Update the maximum likelihood estimates in a fitted model object.  Developer use only. 
}
\usage{
updatepars.msm(x, pars)
}
\arguments{
  \item{x}{A fitted multi-state model object, as returned by
    \code{\link{msm}}.}
  \item{pars}{Vector of new parameters, in their untransformed real-line
    parameterisations, to substitute for the maximum likelihood
    estimates corresponding to those in the \code{estimates} component
    of the fitted model object (\code{\link{msm.object}}).   The
    order of the parameters is  documented in \code{\link{msm}},
    argument \code{fixedpars}. 
}
\value{
  An updated \code{\link{msm}} model object with the updated maximum likelihood
  estimates, but with the covariances / standard errors unchanged.

  Point estimates from output functions such as \code{\link{qmatrix.msm}},
  \code{\link{pmatrix.msm}}, or any related function, can then be evaluated
  with the new parameters, and at arbitrary covariate values.

  This function is used, for example, when computing confidence
  intervals from \code{\link{pmatrix.msm}}, and related functions, using
  the \code{ci="normal"} method. 
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
