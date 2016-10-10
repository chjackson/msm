\name{pmatrix.msm}
\alias{pmatrix.msm}
\title{Transition probability matrix}
\description{
  Extract the estimated transition probability matrix from a fitted multi-state
  model for a given time interval, at a given set of covariate values.
}
\usage{
pmatrix.msm(x=NULL, t=1, t1=0, covariates="mean",
            ci=c("none","normal","bootstrap"), cl=0.95, B=1000,
            cores=NULL, qmatrix=NULL, 
            ...)
}
\arguments{

  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}.}

  \item{t}{The time interval to estimate the transition probabilities
  for, by default one unit. }

\item{t1}{The starting time of
  the interval. Used for models \code{x} with piecewise-constant intensities fitted
  using the \code{pci} option to \code{\link{msm}}. The probabilities will be computed on the interval [t1, t1+t].}

  \item{covariates}{
    The covariate values at which to estimate the transition
    probabilities.  This can either be:\cr

    the string \code{"mean"}, denoting the means of the covariates in
    the data (this is the default),\cr

    the number \code{0}, indicating that all the covariates should be
    set to zero,\cr

    or a list of values, with optional names. For example

    \code{list (60, 1)}

    where the order of the list follows the order of the covariates
    originally given in the model formula, or a named list,

    \code{list (age = 60, sex = 1)}

    For time-inhomogeneous models fitted using the \code{pci} option to
    \code{\link{msm}}, "covariates" here include only those specified using the
    \code{covariates} argument to \code{\link{msm}}, and exclude the
    artificial covariates representing the time period.

    For time-inhomogeneous models fitted "by hand" by using a
    time-dependent covariate in the \code{covariates} argument to
    \code{\link{msm}}, the function \code{\link{pmatrix.piecewise.msm}}
    should be used to to calculate transition probabilities.
  }

  \item{ci}{If \code{"normal"}, then calculate a confidence interval for
    the transition probabilities by simulating \code{B} random vectors
    from the asymptotic multivariate normal distribution implied by the
    maximum likelihood estimates (and covariance matrix) of the log
    transition intensities and covariate effects, then calculating the
    resulting transition probability matrix for each replicate. See,
    e.g. Mandel (2013) for a discussion of this approach.

    If \code{"bootstrap"} then calculate a confidence interval by
    non-parametric bootstrap refitting.  This is 1-2 orders of magnitude
    slower than the \code{"normal"} method, but is expected to be more
    accurate. See \code{\link{boot.msm}} for more details of
    bootstrapping in \pkg{msm}.

    If \code{"none"} (the default) then no confidence interval is
    calculated.}

  \item{cl}{Width of the symmetric confidence interval, relative to 1.}

  \item{B}{Number of bootstrap replicates, or number of normal
    simulations from the distribution of the MLEs}

  \item{cores}{Number of cores to use for bootstrapping using parallel
    processing. See \code{\link{boot.msm}} for more details.}

  \item{qmatrix}{A transition intensity matrix.   Either this or
  a fitted model \code{x} must be supplied.   No confidence intervals
  are available if \code{qmatrix} is supplied.}

  \item{...}{Optional arguments to be passed to \code{\link{MatrixExp}} to
    control the method of computing the matrix exponential.}
}
\value{
  The matrix of estimated transition probabilities \eqn{P(t)} in the given time.
  Rows correspond to "from-state" and columns to "to-state".

  Or if \code{ci="normal"} or \code{ci="bootstrap"}, \code{pmatrix.msm}
returns a list with
  components \code{estimates} and \code{ci}, where \code{estimates} is
  the matrix of estimated transition probabilities, and \code{ci} is a
  list of two matrices containing the upper and lower confidence
  limits.
}
\details{
  For a continuous-time homogeneous Markov process with transition
  intensity matrix
  \eqn{Q},  the probability of occupying state \eqn{s} at time \eqn{u + t}
  conditionally on occupying state \eqn{r} at time \eqn{u} is given by the
  \eqn{(r,s)} entry of the matrix \eqn{P(t) = \exp(tQ)}{P(t) = exp(tQ)},
  where \eqn{\exp()}{exp()} is the matrix exponential.

  For non-homogeneous processes, where covariates and hence the
  transition intensity matrix \eqn{Q} are piecewise-constant in time,
  the transition probability matrix is calculated as
  a product of matrices over a series of intervals, as explained in
  \code{\link{pmatrix.piecewise.msm}}.

  The \code{\link{pmatrix.piecewise.msm}}
  function is only necessary for models fitted using a
  time-dependent covariate in the \code{covariates} argument to
  \code{\link{msm}}. For time-inhomogeneous models fitted using "pci",
  \code{pmatrix.msm} can be used, with arguments \code{t} and \code{t1},
  to calculate transition probabilities over any time period.
}
\references{
  Mandel, M. (2013). "Simulation based confidence intervals for
  functions with complicated derivatives." The American Statistician 67(2):76-81
}
\seealso{
  \code{\link{qmatrix.msm}},  \code{\link{pmatrix.piecewise.msm}}, \code{\link{boot.msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}.}
\keyword{models}
