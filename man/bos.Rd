\name{bos}
\docType{data}
\alias{bos}
\alias{bos3}
\alias{bos4}
\title{Bronchiolitis obliterans syndrome after lung transplants}
\description{
  A dataset containing histories of bronchiolitis obliterans syndrome
  (BOS) from lung transplant recipients. BOS is a chronic decline
  in lung function, often observed after lung transplantation.
  The condition is classified into four stages of severity: none, mild, moderate
  and severe.
}
\usage{bos}
\format{
  A data frame containing 638 rows, grouped by patient, including histories of 204 patients.
  The first observation for each patient is defined to be stage 1, no
  BOS, at six months after transplant.  Subsequent observations denote
  the entry times into stages 2, 3, 4, representing mild, moderate and
  severe BOS respectively, and stage 5, representing death.   
  \tabular{rll}{
    \code{ptnum}   \tab (numeric) \tab Patient identification number \cr
    \code{time}     \tab (numeric) \tab Months after transplant \cr
    \code{state}  \tab (numeric) \tab BOS state entered at this time \cr
}}
\details{
  The entry time of each patient into each stage of BOS was estimated
  by clinicians, based on their history of lung function measurements
  and acute rejection and infection episodes.  BOS is only assumed to occur
  beyond six months after transplant.   In the first six months
  the function of each patient's new lung stabilises.  Subsequently BOS
  is diagnosed by comparing the lung function against the "baseline" value. 

  The objects \code{bos3} and \code{bos4} contain the same data, but
  with mild/moderate/severe combined, and moderate/severe combined,
  to give 3 and 4-state representations respectively.
}
\source{
  Papworth Hospital, U.K.
}
\references{
  Heng. D. et al. (1998).  Bronchiolitis Obliterans Syndrome: Incidence,
  Natural History, Prognosis, and Risk Factors.  Journal of Heart and Lung
  Transplantation 17(12)1255--1263.
}
\keyword{datasets}
