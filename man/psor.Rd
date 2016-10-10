\name{psor}
\docType{data}
\alias{psor}
\title{Psoriatic arthritis data}
\description{
  A series of observations of grades of psoriatic arthritis, as
  indicated by numbers of damaged joints. 
}
\usage{psor}

\format{
  A data frame containing 806 observations, representing visits to a
  psoriatic arthritis (PsA) clinic from 305 patients.  The rows
  are grouped by patient number and ordered by examination time. Each
  row represents an examination and contains additional covariates.
  
  \tabular{rll}{
    \code{ptnum}   \tab (numeric) \tab Patient identification number \cr
    \code{months}     \tab (numeric) \tab  Examination time in months \cr
    \code{state}  \tab (numeric) \tab Clinical state of PsA.  Patients in
    states 1, 2, 3 and 4 \cr \tab \tab have 0, 1 to 4, 5 to 9 and 10 or more damaged
    joints, \cr \tab \tab respectively.  \cr
    \code{hieffusn}   \tab (numeric) \tab  Presence of five or more effusions \cr
    \code{ollwsdrt}   \tab (character) \tab Erythrocyte sedimentation rate of
    less than 15 mm/h \cr
  }}
\references{
  Gladman, D. D. and Farewell, V.T. (1999) Progression in psoriatic arthritis:
  role of time-varying clinical indicators.  J. Rheumatol. 26(11):2409-13
}
\examples{
## Four-state progression-only model with high effusion and low
## sedimentation rate as covariates on the progression rates.  High
## effusion is assumed to have the same effect on the 1-2, 2-3, and 3-4
## progression rates, while low sedimentation rate has the same effect
## on the 1-2 and 2-3 intensities, but a different effect on the 3-4. 

data(psor)
psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
psor.msm <- msm(state ~ months, subject=ptnum, data=psor, 
                qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)),
                fixedpars=FALSE, control = list(REPORT=1,trace=2), method="BFGS")
qmatrix.msm(psor.msm)
sojourn.msm(psor.msm)
hazard.msm(psor.msm)
}
\keyword{datasets}
