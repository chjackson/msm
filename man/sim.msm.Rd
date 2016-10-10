\name{sim.msm}
\title{Simulate one individual trajectory from a continuous-time Markov model}
\alias{sim.msm}
\concept{Simulation}
\description{Simulate one realisation from a continuous-time
  Markov process up to a given time.
}
\usage{
sim.msm(qmatrix, maxtime, covs=NULL, beta=NULL, obstimes=0, start=1,
mintime=0)
}
\arguments{

  \item{qmatrix}{ The transition intensity matrix of the
    Markov process. The diagonal of \code{qmatrix} is ignored,
    and computed as appropriate so that the rows sum to zero. For
    example, a possible \code{qmatrix} for a three state illness-death
    model with recovery is:

    \code{rbind(
    c( 0, 0.1, 0.02 ),
    c( 0.1, 0, 0.01 ),
    c( 0, 0, 0 )
    )}

  }

  \item{maxtime}{ Maximum time for the simulated process. }

  \item{covs}{Matrix of time-dependent covariates, with one row for each
  observation time and one column for each covariate.}

  \item{beta}{Matrix of linear covariate effects on log transition
    intensities. The rows correspond to different covariates, and the
    columns to the transition intensities. The intensities are ordered
    by reading across rows of the intensity matrix, starting with the
    first, counting the positive off-diagonal elements of the matrix.}

  \item{obstimes}{Vector of times at which the covariates are observed.}

  \item{start}{Starting state of the process. Defaults to 1.}

  \item{mintime}{Starting time of the process. Defaults to 0.}

}
\value{ A list with components,

  \item{states}{Simulated states through which the process moves.  This
    ends with either an absorption before \code{obstime}, or a transient state
    at \code{obstime}. }

  \item{times}{Exact times at which the process changes to the corresponding
    states}

  \item{qmatrix}{The given transition intensity matrix}
}
\details{
  The effect of time-dependent covariates on the transition intensity
  matrix for an individual is determined by assuming that the covariate is a step function
  which remains constant in between the individual's observation times.
}
\seealso{ \code{\link{simmulti.msm}}
}
\examples{

qmatrix <- rbind(
                 c(-0.2,   0.1,  0.1 ),
                 c(0.5,   -0.6,  0.1 ),
                 c(0,  0,  0)
                 )
sim.msm(qmatrix, 30)
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{datagen}
