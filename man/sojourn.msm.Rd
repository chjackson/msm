\name{sojourn.msm}
\title{Mean sojourn times from a multi-state model}
\alias{sojourn.msm}
\description{
  Estimate the mean sojourn times in the transient states of a
  multi-state model and their confidence limits. 
}
\usage{
sojourn.msm(x, covariates="mean", ci=c("delta","normal","bootstrap","none"),
  cl=0.95, B=1000)
}
\arguments{

  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}.}

  \item{covariates}{
    The covariate values at which to estimate the mean sojourn times. This can either be:\cr

    the string \code{"mean"}, denoting the means of the covariates in
    the data (this is the default),\cr

    the number \code{0}, indicating that all the covariates should be
    set to zero,\cr

    a list of values, with optional names. For example,

    \code{list(60, 1)}, 
    where the order of the list follows the order of the covariates
    originally given in the model formula, or a named list, e.g.
    
    \code{list (age = 60, sex = 1)}
   
  }

  \item{ci}{
    If \code{"delta"} (the default) then confidence intervals are
    calculated by the delta method, or by simple transformation of the
    Hessian in the very simplest cases. 
    
    If \code{"normal"}, then calculate a confidence interval by simulating \code{B} random vectors
    from the asymptotic multivariate normal distribution implied by the
    maximum likelihood estimates (and covariance matrix) of the log
    transition intensities and covariate effects, then transforming.
    
    If \code{"bootstrap"} then calculate a confidence interval by
    non-parametric bootstrap refitting.  This is 1-2 orders of magnitude
    slower than the \code{"normal"} method, but is expected to be more
    accurate. See \code{\link{boot.msm}} for more details of
    bootstrapping in \pkg{msm}.}
    
  \item{cl}{Width of the symmetric confidence interval to present. 
    Defaults to 0.95.}

  \item{B}{Number of bootstrap replicates, or number of normal
    simulations from the distribution of the MLEs}
  
}
\value{
    A data frame with components:

    \item{estimates}{Estimated mean sojourn times in the transient states.}
    \item{SE}{Corresponding standard errors.}
    \item{L}{Lower confidence limits.}
    \item{U}{Upper confidence limits.}
}
\details{
  The mean sojourn time in a transient state \eqn{r} is estimated by
  \eqn{- 1 / q_{rr}}, where \eqn{q_{rr}} is the \eqn{r}th entry on the diagonal
  of the estimated transition intensity matrix.

  A continuous-time Markov model is fully specified by the mean sojourn
  times and the probability that each state is next (\code{\link{pnext.msm}}).
  This is a more intuitively meaningful description of a model than the transition intensity
  matrix (\code{\link{qmatrix.msm}}).

  Time dependent covariates, or time-inhomogeneous models, are not
  supported.  This would require the mean of a piecewise exponential
  distribution, and the package author is not aware of any general
  analytic form for that.
}
\seealso{
  \code{\link{msm}}, \code{\link{qmatrix.msm}}, \code{\link{deltamethod}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
