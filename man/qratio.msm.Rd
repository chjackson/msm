\name{qratio.msm}
\alias{qratio.msm}
\title{Estimated ratio of transition intensities}
\description{
  Compute the estimate and approximate standard error of the ratio of
  two estimated transition intensities from a fitted multi-state model
  at a given set of covariate values.
}
\usage{
qratio.msm(x, ind1, ind2, covariates = "mean",
           ci=c("delta","normal","bootstrap","none"), cl = 0.95,
           B=1000, cores=NULL)
}
\arguments{
  \item{x}{A fitted multi-state model, as returned by
    \code{\link{msm}}.}

  \item{ind1}{Pair of numbers giving the indices in the intensity matrix
    of the numerator of the ratio, for example, \code{c(1,2)}. }

  \item{ind2}{Pair of numbers giving the indices in the intensity matrix
    of the denominator of the ratio, for example, \code{c(2,1)}. }

  \item{covariates}{
    The covariate values at which to estimate the intensities.
    This can either be:\cr

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
    calculated by the delta method.

    If \code{"normal"}, then calculate a confidence interval by
    simulating \code{B} random vectors
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

  \item{cores}{Number of cores to use for bootstrapping using parallel
    processing. See \code{\link{boot.msm}} for more details.}
}
\value{
  A named vector with elements \code{estimate}, \code{se}, \code{L} and \code{U}
  containing the estimate, standard error, lower and upper confidence
  limits, respectively, of the ratio of intensities.
}
\details{

  For example, we might want to compute the ratio of the progression
  rate and recovery rate for a fitted model \code{disease.msm} with a
  health state (state 1) and a disease state
  (state 2).  In this case, the progression rate is the (1,2) entry of
  the intensity matrix, and the recovery rate is the (2,1) entry.
  Thus to compute this ratio with covariates set to their means, we
  call

  \code{qratio.msm(disease.msm, c(1,2), c(2,1))} .

  Standard errors are estimated by the delta method. Confidence limits
  are estimated by assuming normality on the log scale.

}
\seealso{
  \code{\link{qmatrix.msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
