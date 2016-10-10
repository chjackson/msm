\name{msm2Surv}
\alias{msm2Surv}
\title{Convert data for `msm' to data for `survival', `mstate' or `flexsurv' analysis}
\description{
  Converts longitudinal data for a \code{\link{msm}} model fit, where observations
  represent the exact transition times of the process, to counting process
  data.
  This enables, for example, flexible parametric multi-state models to
  be fitted with \code{\link[flexsurv]{flexsurvreg}} from the
  \pkg{flexsurv} package, or semiparametric models to be implemented with
  \code{\link[survival]{coxph}} and the \pkg{mstate} package.
}
\usage{
msm2Surv(data, subject, time, state, covs, Q)
}
\arguments{
  \item{data}{Data frame in the format expected by a \code{\link{msm}}
    model fit with \code{exacttimes=TRUE} or all \code{obstype=2}.  
    Each row represents an observation of a state, and the time variable
    contains the exact and complete transition times of the underlying process. 
    This is explained in more detail
    in the help page for \code{\link{msm}}, section \code{obstype=2}. 
  }
  \item{subject}{Name of the subject ID in the data (character format,
    i.e. quoted).}
  \item{time}{Name of the time variable in the data (character).}
  \item{state}{Name of the state variable in the data (character).}
  \item{covs}{Vector of covariate names to carry through (character).
    If not supplied, this is taken to be all remaining variables in the
    data.
  }
  \item{Q}{Transition intensity matrix.  This should have number of rows
    and number of columns both equal to the number of states.  If an
    instantaneous transition is not allowed from state \eqn{r} to state
    \eqn{s}, then \code{Q} should have \eqn{(r,s)} entry 0, otherwise it
    should be non-zero.  The diagonal entries are ignored.}
}
\value{
  A data frame of class \code{"msdata"}, with rows representing observed or censored transitions. 
  There will be one row for each observed transition in the original
    data, and additional rows for every potential transition that could have
    occurred out of each observed state.

    The data frame will have columns called:

    \item{id}{Subject ID}
    \item{from}{Starting state of the transition}
    \item{to}{Finishing state of the transition}
    \item{Tstart}{The starting time of the transition}
    \item{Tstop}{The finishing time of the transition}
    \item{time}{The time difference = \code{Tstop} - \code{Tstart}}
    \item{status}{Event or censoring indicator, with 1 indicating an
    observed transition, and 0 indicating censoring}
    \item{trans}{Transition number} 

    and any remaining columns will represent covariates.  Any covariates
    whose names clash with the standard variables in the returned data
    (\code{"id"}, \code{"from"}, \code{"to"}, \code{"Tstart"}, \code{"Tstop"}, \code{"time"}, \code{"status"} or \code{"trans"}) have \code{".2"} appended to their names.

     The transition matrix in \pkg{mstate} format is stored in the \code{trans}
     attribute of the returned object.  See the example code below. 
   }
\details{
  For example, if the data supplied to \code{\link{msm}} look like this:

  \tabular{lllll}{
    \code{subj} \tab \code{days} \tab \code{status}  \tab \code{age} \tab \code{treat} \cr
       1\tab       0\tab     1  \tab 66\tab  1\cr
       1\tab      27\tab     2  \tab 66\tab  1\cr
       1\tab      75\tab     3  \tab 66\tab  1\cr
       1\tab      97\tab     4  \tab 66\tab  1\cr
       1\tab    1106\tab     4  \tab 69\tab  1\cr
       2\tab       0\tab     1  \tab 49\tab  0\cr
       2\tab      90\tab     2  \tab 49\tab  0\cr
       2\tab    1037\tab     2  \tab 51\tab  0\cr
  }

  then the output of \code{\link{msm2Surv}} will be a data frame looking like this:

  \tabular{lllllllllll}{
\code{id}  \tab  \code{from}  \tab  \code{to}  \tab  \code{Tstart}  \tab  \code{Tstop}  \tab  \code{time}  \tab  \code{status}  \tab  \code{age}  \tab  \code{treat}  \tab  \code{trans}\cr
1  \tab  1  \tab  2  \tab  0  \tab  27  \tab  27  \tab  1  \tab  66  \tab  1  \tab  1\cr
1  \tab  1  \tab  4  \tab  0  \tab  27  \tab  27  \tab  0  \tab  66  \tab  1  \tab  2\cr
1  \tab  2  \tab  3  \tab  27  \tab  75  \tab  48  \tab  1  \tab  66  \tab  1  \tab  3\cr
1  \tab  2  \tab  4  \tab  27  \tab  75  \tab  48  \tab  0  \tab  66  \tab  1  \tab  4\cr
1  \tab  3  \tab  4  \tab  75  \tab  97  \tab  22  \tab  1  \tab  69  \tab  1  \tab  5\cr
2  \tab  1  \tab  2  \tab  0  \tab  90  \tab  90  \tab  1  \tab  49  \tab  0  \tab  1\cr
2  \tab  1  \tab  4  \tab  0  \tab  90  \tab  90  \tab  0  \tab  49  \tab  0  \tab  2\cr
2  \tab  2  \tab  3  \tab  90  \tab  1037  \tab  947  \tab  0  \tab  49  \tab  0  \tab  3\cr
2  \tab  2  \tab  4  \tab  90  \tab  1037  \tab  947  \tab  0\tab  49  \tab  0  \tab  4\cr
}

At 27 days, subject 1 is observed to move from state 1 to
state 2 (first row, status 1), which means that their potential
transition from state 1 to state 4 is censored (second row, status 0).

See the \pkg{mstate} package and the references below for more details
of this data format and using it for semi-parametric multi-state
modelling. 
}
\references{
   Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
   Competing risks and multi-state models. \emph{Statistics in Medicine}
   26: 2389-2430.

   Liesbeth C. de Wreede, Marta Fiocco, Hein Putter (2011). \pkg{mstate}: An R
  Package for the Analysis of Competing Risks and Multi-State
  Models. \emph{Journal of Statistical Software}, 38(7), 1-30. 
  \url{http://www.jstatsoft.org/v38/i07}  

  Jackson, C. H. (2014). flexsurv: Flexible parametric survival and
  multi-state models.  R package version 0.5.

}
\seealso{
  \code{\link[mstate]{msprep}}, in \pkg{mstate}, which produces data in a similar format,
  given data in "wide" format with one row per subject.
}
\examples{
msmdat <- data.frame(
 subj = c(1, 1, 1, 1, 1, 2, 2, 2),
 days = c(0, 27, 75, 97, 1106, 0, 90, 1037),
 status = c(1, 2, 3, 4, 4, 1, 2, 2),
 age = c(66, 66, 66, 66, 69, 49, 49, 51),
 treat = c(1, 1, 1, 1, 1, 0, 0, 0)
)
# transitions only allowed to next state up or state 4
Q <- rbind(c(1, 1, 0, 1), 
           c(0, 1, 1, 1),
           c(0, 0, 1, 1),
           c(0, 0, 0, 0))
dat <- msm2Surv(data=msmdat, subject="subj", time="days", state="status", 
         Q=Q)
dat
attr(dat, "trans")
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
