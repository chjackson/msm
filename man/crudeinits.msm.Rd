\name{crudeinits.msm}
\alias{crudeinits.msm}
\title{Calculate crude initial values for transition intensities}
\description{
  Calculates crude initial values for transition intensities by assuming
  that the data represent the exact transition times of the Markov
  process.
}
\usage{
crudeinits.msm(formula, subject, qmatrix, data=NULL, censor=NULL, censor.states=NULL)
}
\arguments{
  \item{formula}{ A formula giving the vectors containing the observed
    states and the  corresponding observation times. For example,
    
    \code{state ~ time}

    Observed states should be in the set \code{1, \dots, n},
    where \code{n} is the number of states.  Note hidden Markov models
    are not supported by this function.
  }
  \item{subject}{Vector of subject identification numbers for the data
    specified by \code{formula}. If missing, then all observations
    are assumed to be on the same subject. These must be sorted so that
    all observations on the same subject are adjacent.}

  \item{qmatrix}{Matrix of indicators for the allowed
    transitions.  An initial value will be estimated for each
    value of qmatrix that is greater than zero.  Transitions are
    taken as disallowed for each entry of \code{qmatrix} that is 0. }

  \item{data}{An optional data frame in which the variables represented
    by \code{subject} and \code{state} can be found.}

  \item{censor}{ A state, or vector of states, which indicates
    censoring.  See \code{\link{msm}}. 
  }

  \item{censor.states}{
    Specifies the underlying states which censored observations can
    represent.  See \code{\link{msm}}. 
  }
}
\value{
  The estimated transition intensity matrix.  This can be used as
  the \code{qmatrix} argument to \code{\link{msm}}.  
}
\details{

  Suppose we want a crude estimate of the transition intensity \eqn{q_{rs}}{q_rs} from state \eqn{r} to
  state \eqn{s}.   If we observe \eqn{n_{rs}}{n_rs} transitions from state \eqn{r} to state \eqn{s}, and a
  total of \eqn{n_r} transitions from state \eqn{r}, then \eqn{q_{rs} /
    q_{rr}}{q_rs / q_rr} can be estimated by \eqn{n_{rs} / n_r}{n_rs / n_r}. Then, given a total of \eqn{T_r} years
  spent in state \eqn{r},  the mean sojourn time \eqn{1 / q_{rr}}{1 / q_rr} can be estimated as \eqn{T_r / n_r}.
  Thus, \eqn{n_{rs} / T_r}{n_rs / T_r} is a crude estimate of
  \eqn{q_{rs}}{q_rs}.

  If the data do represent the exact transition times of the Markov
  process, then these are the exact maximum likelihood estimates.

  Observed transitions which are incompatible with the given
  \code{qmatrix} are ignored.  Censored states are ignored. 

}
\seealso{
  \code{\link{statetable.msm}}
}
\examples{
data(cav)
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
statetable.msm(state, PTNUM, data=cav)
crudeinits.msm(state ~ years, PTNUM, data=cav, qmatrix=twoway4.q)
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
