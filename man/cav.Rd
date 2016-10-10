\name{cav}
\docType{data}
\alias{cav}
\title{Heart transplant monitoring data}
\description{
  A series of approximately yearly angiographic examinations of heart transplant
  recipients.  The state at each time is a grade of cardiac allograft 
  vasculopathy (CAV),  a deterioration of the arterial walls.  
}
\usage{cav}

\format{
  A data frame containing 2846 rows.  There are 622 patients, the rows
  are grouped by patient number and ordered by years after transplant,
  with each row representing an examination and containing additional
  covariates.
  
  \tabular{rll}{
    \code{PTNUM}   \tab (numeric) \tab Patient identification number \cr
    \code{age}     \tab (numeric) \tab Recipient age at examination (years) \cr
    \code{years}   \tab (numeric) \tab Examination time (years after transplant)\cr
    \code{dage}    \tab (numeric) \tab Age of heart donor (years) \cr
    \code{sex}   \tab (numeric) \tab sex (0=male, 1=female) \cr
    \code{pdiag} \tab (factor) \tab Primary diagnosis (reason for transplant)  \cr \tab \tab IHD=ischaemic heart disease, IDC=idiopathic dilated  cardiomyopathy. \cr
    \code{cumrej}     \tab (numeric) \tab Cumulative number of acute rejection  episodes \cr
    \code{state} \tab (numeric) \tab State at the examination. \cr \tab
    \tab State 1 represents no CAV, state 2 is mild/moderate CAV \cr \tab \tab and state 3 is severe CAV.  State 4 indicates death.   \cr
    \code{firstobs}   \tab (numeric) \tab 0 = record represents an angiogram or date of death.\cr \tab \tab
    1 = record represents transplant (patient's first observation) \cr
    \code{statemax}  \tab (numeric) \tab Maximum observed state so far for this patient (added in version 1.5.1) 
  }}
\source{
  Papworth Hospital, U.K.
}
\references{
  Sharples, L.D. and Jackson, C.H. and Parameshwar, J. and Wallwork,
  J. and Large, S.R. (2003). Diagnostic accuracy of coronary angiopathy
  and risk factors for post-heart-transplant cardiac allograft
  vasculopathy. Transplantation 76(4):679-82
}
\keyword{datasets}
