\name{statetable.msm}
\alias{statetable.msm}
\title{Table of transitions}
\description{
  Calculates a frequency table counting the number of times each pair of
  states were observed in successive observation times.  This can be a
  useful way of summarising multi-state data.
}
\usage{
statetable.msm(state, subject, data=NULL)
}
\arguments{
  \item{state}{Observed states, assumed to be ordered by time
  within each subject.}
  \item{subject}{Subject identification numbers corresponding to
    \code{state}. If not given, all observations are assumed to be on
    the same subject.}
  \item{data}{An optional data frame in which the variables represented
    by \code{subject} and \code{state} can be found.}
}
\value{
  A frequency table with starting states as rows and finishing states as
  columns.
}
\details{
  If the data are intermittently observed (panel data) this table should
  not be used to decide what transitions should be allowed in the
  \eqn{Q} matrix, which works in continuous time.  This function counts
  the transitions between states over a time interval, not in real time.
  There can be observed transitions between state \eqn{r} and \eqn{s}
  over an interval even if \eqn{q_{rs}=0}, because the process may have
  passed through one or more intermediate states in the middle of the
  interval.  
}
\seealso{
  \code{\link{crudeinits.msm}}
}
\examples{
## Heart transplant data
data(cav)

## 148 deaths from state 1, 48 from state 2 and 55 from state 3.
statetable.msm(state, PTNUM, data=cav)

}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
