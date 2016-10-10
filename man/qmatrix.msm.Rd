\name{qmatrix.msm}
\alias{qmatrix.msm}
\title{Transition intensity matrix}
\description{
  Extract the estimated transition intensity matrix, and the
  corresponding standard errors, from a fitted
  multi-state model at a given set of covariate values.
}
\usage{
qmatrix.msm(x, covariates="mean", sojourn=FALSE,
            ci=c("delta","normal","bootstrap","none"), cl=0.95,
            B=1000, cores=NULL)
}
\arguments{

  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}.}

  \item{covariates}{
    The covariate values at which to estimate the intensity matrix.
    This can either be:\cr

    the string \code{"mean"}, denoting the means of the covariates in
    the data (this is the default),\cr

    the number \code{0}, indicating that all the covariates should be
    set to zero,\cr

    or a list of values, with optional names. For example

    \code{list (60, 1)}

    where the order of the list follows the order of the covariates
    originally given in the model formula. Or more clearly, a named list,

    \code{list (age = 60, sex = 1)}

    If some covariates are specified but not others, the missing ones
    default to zero.
  }
  \item{sojourn}{Set to TRUE if the estimated sojourn times and their
    standard errors should also be returned.}

  \item{ci}{
    If \code{"delta"} (the default) then confidence intervals are
    calculated by the delta method, or by simple transformation of the
    Hessian in the very simplest cases.  Normality on the log scale
    is assumed.

    If \code{"normal"}, then calculate a confidence interval by
    simulating \code{B} random vectors
    from the asymptotic multivariate normal distribution implied by the
    maximum likelihood estimates (and covariance matrix) of the log
    transition intensities and covariate effects, then transforming.

    If \code{"bootstrap"} then calculate a confidence interval by
    non-parametric bootstrap refitting.  This is 1-2 orders of magnitude
    slower than the \code{"normal"} method, but is expected to be more
    accurate. See \code{\link{boot.msm}} for more details of
    bootstrapping in \pkg{msm}.
  }

  \item{cl}{Width of the symmetric confidence interval to present.
    Defaults to 0.95.}

  \item{B}{Number of bootstrap replicates, or number of normal
    simulations from the distribution of the MLEs.}

  \item{cores}{Number of cores to use for bootstrapping using parallel
    processing. See \code{\link{boot.msm}} for more details.}
  }

\value{
  A list with components:

  \item{estimate}{Estimated transition intensity matrix.}

  \item{SE}{Corresponding approximate standard errors.}

  \item{L}{Lower confidence limits}

  \item{U}{Upper confidence limits}

  Or if \code{ci="none"}, then \code{qmatrix.msm} just returns the
  estimated transition intensity matrix.

  If \code{sojourn} is \code{TRUE}, extra components called
  \code{sojourn}, \code{sojournSE}, \code{sojournL} and \code{sojournU} are included, containing the
  estimates, standard errors and confidence limits, respectively, of the
  mean sojourn times in each transient state.

  The default print method for objects returned by
  \code{\link{qmatrix.msm}} presents estimates and confidence limits. To
  present estimates and standard errors, do something like

  \code{qmatrix.msm(x)[c("estimates","SE")]}
}

\details{
  Transition intensities and covariate effects are estimated on the log
  scale by \code{\link{msm}}. A covariance matrix is estimated from the
  Hessian of the maximised log-likelihood.

  A more practically meaningful parameterisation of a continuous-time Markov model
  with transition intensities \eqn{q_{rs}} is in terms of the mean
  sojourn times \eqn{-1 / q_{rr}} in each state \eqn{r} and the
  probabilities that the next move of the process when in state \eqn{r}
  is to state \eqn{s}, \eqn{-q_{rs} / q_{rr}}.
}

\seealso{
  \code{\link{pmatrix.msm}}, \code{\link{sojourn.msm}},
  \code{\link{deltamethod}}, \code{\link{ematrix.msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
