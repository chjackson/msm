\name{phasemeans.msm}
\alias{phasemeans.msm}
\title{Parameters of phase-type models in mixture form}
\description{
  Parameters of fitted two-phase models, in mixture model parameterisation.
}
\usage{
phasemeans.msm(x, covariates="mean", ci=c("none","normal","bootstrap"),
               cl=0.95, B=1000, cores=NULL)
}
\arguments{
  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}.}
  \item{covariates}{Covariate values, see \code{\link{qmatrix.msm}}.}
  \item{ci}{
    If \code{"none"} (the default) no confidence intervals are
    calculated.   Otherwise \code{"normal"}, or \code{"boot"}
    as described by \code{\link{qmatrix.msm}}.}

  \item{cl}{Width of the symmetric confidence interval, relative to 1.}

  \item{B}{Number of bootstrap replicates, or number of normal
    simulations from the distribution of the MLEs.}

  \item{cores}{Number of cores to use for bootstrapping using parallel
    processing. See \code{\link{boot.msm}} for more details.}
}
\value{
  Matrix with one row for each state that has a two-phase distribution,
  and three columns: the short-stay mean, long-stay mean and long-stay
  probability.  These are functions of the transition intensities of the
  expanded hidden Markov model, defined in \code{\link{d2phase}}.
}
\seealso{
  \code{\link{d2phase}}.
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
