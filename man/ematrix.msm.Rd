\name{ematrix.msm}
\alias{ematrix.msm}
\title{Misclassification probability matrix}
\description{
  Extract the estimated misclassification probability matrix, and
  corresponding confidence intervals,  from a fitted multi-state
  model at a given set of covariate values.
}
\usage{
ematrix.msm(x, covariates="mean", ci=c("delta","normal","bootstrap","none"),
            cl=0.95, B=1000, cores=NULL)
}
\arguments{

  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}}.

  \item{covariates}{

    The covariate values for which to estimate the misclassification
    probability matrix.  This can either be:\cr

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

  \item{ci}{
    If \code{"delta"} (the default) then confidence intervals are
    calculated by the delta method, or by simple transformation of the
    Hessian in the very simplest cases.

    If \code{"normal"}, then calculate a confidence interval by simulating \code{B} random vectors
    from the asymptotic multivariate normal distribution implied by the
    maximum likelihood estimates (and covariance matrix) of the
    multinomial-logit-transformed misclassification probabilities and
    covariate effects, then transforming back.

    If \code{"bootstrap"} then calculate a confidence interval by
    non-parametric bootstrap refitting.  This is 1-2 orders of magnitude
    slower than the \code{"normal"} method, but is expected to be more
    accurate. See \code{\link{boot.msm}} for more details of
    bootstrapping in \pkg{msm}.}

  \item{cl}{Width of the symmetric confidence interval to present.
    Defaults to 0.95.}

  \item{B}{Number of bootstrap replicates, or number of normal
    simulations from the distribution of the MLEs}

  \item{cores}{Number of cores to use for bootstrapping using parallel
    processing. See \code{\link{boot.msm}} for more details.}

}
\value{
  A list with components:

  \item{estimate}{Estimated misclassification probability matrix. The
    rows correspond to true states, and columns observed states.}
  \item{SE}{Corresponding approximate standard errors.}
  \item{L}{Lower confidence limits.}
  \item{U}{Upper confidence limits.}

  Or if \code{ci="none"}, then \code{ematrix.msm} just returns the
  estimated misclassification probability matrix.

  The default print method for objects returned by
  \code{\link{ematrix.msm}} presents estimates and confidence limits. To
  present estimates and standard errors, do something like

  \code{ematrix.msm(x)[c("estimates","SE")]}

}
\details{
  Misclassification probabilities and covariate effects are estimated on
  the multinomial-logit scale by \code{\link{msm}}. A covariance matrix
  is estimated from the Hessian of the maximised log-likelihood.  From
  these, the delta method can be used to obtain standard errors of the
  probabilities on the natural scale at arbitrary covariate values.
  Confidence intervals are estimated by assuming normality on the
  multinomial-logit scale.
}
\seealso{
  \code{\link{qmatrix.msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
