\name{pnext.msm}
\alias{pnext.msm}
\title{Probability of each state being next}
\description{
  Compute a matrix of the probability of each state \eqn{s} being the
  next state of the process after each state \eqn{r}.  Together with
  the mean sojourn times in each state (\code{\link{sojourn.msm}}),
  these fully define a continuous-time Markov model.
}
\usage{
pnext.msm(x, covariates = "mean",
           ci=c("normal","bootstrap","delta","none"), cl = 0.95,
           B=1000, cores=NULL)
}
\arguments{
  \item{x}{A fitted multi-state model, as returned by
    \code{\link{msm}}.}

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
    If \code{"normal"} (the default) then calculate a confidence interval by
    simulating \code{B} random vectors
    from the asymptotic multivariate normal distribution implied by the
    maximum likelihood estimates (and covariance matrix) of the log
    transition intensities and covariate effects, then transforming.

    If \code{"bootstrap"} then calculate a confidence interval by
    non-parametric bootstrap refitting.  This is 1-2 orders of magnitude
    slower than the \code{"normal"} method, but is expected to be more
    accurate. See \code{\link{boot.msm}} for more details of
    bootstrapping in \pkg{msm}.

    If \code{"delta"} then confidence intervals are
    calculated based on the delta method SEs of the log rates, but this is
    not recommended since it may not respect the constraint 
    that probabilities are less than one.
  }

  \item{cl}{Width of the symmetric confidence interval to present.
    Defaults to 0.95.}

  \item{B}{Number of bootstrap replicates, or number of normal
    simulations from the distribution of the MLEs.}

  \item{cores}{Number of cores to use for bootstrapping using parallel
    processing. See \code{\link{boot.msm}} for more details.}
}
\value{
  The matrix of probabilities that the next move of a process in state
    \eqn{r} (rows) is to state \eqn{s} (columns).
}
\details{
  For a continuous-time Markov process in state \eqn{r}, the probability
  that the next state is \eqn{s} is \eqn{-q_{rs} / q_{rr}}, where
  \eqn{q_{rs}} is the transition intensity (\code{\link{qmatrix.msm}}).

  A continuous-time Markov model is fully specified by these probabilities together with
  the mean sojourn times \eqn{-1/q_{rr}} in each state \eqn{r}. This
  gives a more intuitively meaningful description of a model than the
  intensity matrix.

  Remember that \pkg{msm} deals with continuous-time, not discrete-time
  models, so these are \emph{not} the same as the probability of observing
  state \eqn{s} at a fixed time in the future.  Those probabilities are
  given by \code{\link{pmatrix.msm}}.
}
\seealso{
  \code{\link{qmatrix.msm}},\code{\link{pmatrix.msm}},\code{\link{qratio.msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
