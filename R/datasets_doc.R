#' Aortic aneurysm progression data
#' 
#' This dataset contains longitudinal measurements of grades of aortic
#' aneurysms, measured by ultrasound examination of the diameter of the aorta.
#' 
#' 
#' @name aneur
#' @docType data
#' @format A data frame containing 4337 rows, with each row corresponding to an
#' ultrasound scan from one of 838 men over 65 years of age.
#' 
#' \tabular{rll}{ \code{ptnum} \tab (numeric) \tab Patient identification
#' number \cr \code{age} \tab (numeric) \tab Recipient age at examination
#' (years) \cr \code{diam} \tab (numeric) \tab Aortic diameter\cr \code{state}
#' \tab (numeric) \tab State of aneurysm. \cr }
#' 
#' The states represent successive degrees of aneurysm severity, as indicated
#' by the aortic diameter.
#' 
#' \tabular{rll}{ State 1 \tab Aneurysm-free \tab < 30 cm \cr State 2 \tab Mild
#' aneurysm \tab 30-44 cm \cr State 3 \tab Moderate aneurysm \tab 45-54 cm \cr
#' State 4 \tab Severe aneurysm \tab > 55 cm \cr }
#' 
#' 683 of these men were aneurysm-free at age 65 and were re-screened every two
#' years.  The remaining men were aneurysmal at entry and had successive
#' screens with frequency depending on the state of the aneurysm.  Severe
#' aneurysms are repaired by surgery.
#' @references Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W.
#' and Couto, E.  Multi-state Markov models for disease progression with
#' classification error. \emph{The Statistician}, 52(2): 193--209 (2003)
#' 
#' Couto, E. and Duffy, S. W. and Ashton, H. A. and Walker, N. M.  and Myles,
#' J. P. and Scott, R. A. P. and Thompson, S. G. (2002) \emph{Probabilities of
#' progression of aortic aneurysms: estimates and implications for screening
#' policy} Journal of Medical Screening 9(1):40--42
#' @source The Chichester, U.K. randomised controlled trial of screening for
#' abdominal aortic aneurysms by ultrasonography.
#' @keywords datasets
NULL





#' Bronchiolitis obliterans syndrome after lung transplants
#' 
#' A dataset containing histories of bronchiolitis obliterans syndrome (BOS)
#' from lung transplant recipients. BOS is a chronic decline in lung function,
#' often observed after lung transplantation.  The condition is classified into
#' four stages of severity: none, mild, moderate and severe.
#' 
#' The entry time of each patient into each stage of BOS was estimated by
#' clinicians, based on their history of lung function measurements and acute
#' rejection and infection episodes.  BOS is only assumed to occur beyond six
#' months after transplant.  In the first six months the function of each
#' patient's new lung stabilises.  Subsequently BOS is diagnosed by comparing
#' the lung function against the "baseline" value.
#' 
#' The objects \code{bos3} and \code{bos4} contain the same data, but with
#' mild/moderate/severe combined, and moderate/severe combined, to give 3 and
#' 4-state representations respectively.
#' 
#' @name bos
#' @aliases bos bos3 bos4
#' @docType data
#' @format A data frame containing 638 rows, grouped by patient, including
#' histories of 204 patients.  The first observation for each patient is
#' defined to be stage 1, no BOS, at six months after transplant.  Subsequent
#' observations denote the entry times into stages 2, 3, 4, representing mild,
#' moderate and severe BOS respectively, and stage 5, representing death.
#' \tabular{rll}{ \code{ptnum} \tab (numeric) \tab Patient identification
#' number \cr \code{time} \tab (numeric) \tab Months after transplant \cr
#' \code{state} \tab (numeric) \tab BOS state entered at this time \cr }
#' @references Heng. D. et al. (1998).  Bronchiolitis Obliterans Syndrome:
#' Incidence, Natural History, Prognosis, and Risk Factors.  Journal of Heart
#' and Lung Transplantation 17(12)1255--1263.
#' @source Papworth Hospital, U.K.
#' @keywords datasets
NULL





#' Heart transplant monitoring data
#' 
#' A series of approximately yearly angiographic examinations of heart
#' transplant recipients.  The state at each time is a grade of cardiac
#' allograft vasculopathy (CAV), a deterioration of the arterial walls.
#' 
#' 
#' @name cav
#' @docType data
#' @format A data frame containing 2846 rows.  There are 622 patients, the rows
#' are grouped by patient number and ordered by years after transplant, with
#' each row representing an examination and containing additional covariates.
#' 
#' \tabular{rll}{ \code{PTNUM} \tab (numeric) \tab Patient identification
#' number \cr \code{age} \tab (numeric) \tab Recipient age at examination
#' (years) \cr \code{years} \tab (numeric) \tab Examination time (years after
#' transplant)\cr \code{dage} \tab (numeric) \tab Age of heart donor (years)
#' \cr \code{sex} \tab (numeric) \tab sex (0=male, 1=female) \cr \code{pdiag}
#' \tab (factor) \tab Primary diagnosis (reason for transplant) \cr \tab \tab
#' IHD=ischaemic heart disease, IDC=idiopathic dilated cardiomyopathy. \cr
#' \code{cumrej} \tab (numeric) \tab Cumulative number of acute rejection
#' episodes \cr \code{state} \tab (numeric) \tab State at the examination. \cr
#' \tab \tab State 1 represents no CAV, state 2 is mild/moderate CAV \cr \tab
#' \tab and state 3 is severe CAV.  State 4 indicates death.  \cr
#' \code{firstobs} \tab (numeric) \tab 0 = record represents an angiogram or
#' date of death.\cr \tab \tab 1 = record represents transplant (patient's
#' first observation) \cr \code{statemax} \tab (numeric) \tab Maximum observed
#' state so far for this patient (added in version 1.5.1) }
#' @references Sharples, L.D. and Jackson, C.H. and Parameshwar, J. and
#' Wallwork, J. and Large, S.R. (2003). Diagnostic accuracy of coronary
#' angiopathy and risk factors for post-heart-transplant cardiac allograft
#' vasculopathy. Transplantation 76(4):679-82
#' @source Papworth Hospital, U.K.
#' @keywords datasets
NULL







#' FEV1 measurements from lung transplant recipients
#' 
#' A series of measurements of the forced expiratory volume in one second
#' (FEV1) from lung transplant recipients, from six months onwards after their
#' transplant.
#' 
#' A baseline "normal" FEV1 for each individual is calculated using
#' measurements from the first six months after transplant. After six months,
#' as presented in this dataset, FEV1 is expressed as a percentage of the
#' baseline value.
#' 
#' FEV1 is monitored to diagnose bronchiolitis obliterans syndrome (BOS), a
#' long-term lung function decline, thought to be a form of chronic rejection.
#' Acute rejections and infections also affect the lung function in the short
#' term.
#' 
#' @name fev
#' @docType data
#' @format A data frame containing 5896 rows.  There are 204 patients, the rows
#' are grouped by patient number and ordered by days after transplant.  Each
#' row represents an examination and containing an additional covariate.
#' 
#' \tabular{rll}{ \code{ptnum} \tab (numeric) \tab Patient identification
#' number. \cr \code{days} \tab (numeric) \tab Examination time (days after
#' transplant). \cr \code{fev} \tab (numeric) \tab Percentage of baseline FEV1.
#' A code of 999 indicates the patient's date of death. \cr \code{acute} \tab
#' (numeric) \tab 0/1 indicator for whether the patient suffered an acute
#' infection or rejection \cr \tab \tab within 14 days of the visit.  \cr }
#' @references Jackson, C.H. and Sharples, L.D. Hidden Markov models for the
#' onset and progression of bronchiolitis obliterans syndrome in lung
#' transplant recipients \emph{Statistics in Medicine}, 21(1): 113--128 (2002).
#' @source Papworth Hospital, U.K.
#' @keywords datasets
NULL

#' Psoriatic arthritis data
#' 
#' A series of observations of grades of psoriatic arthritis, as indicated by
#' numbers of damaged joints.
#' 
#' 
#' @name psor
#' @docType data
#' @format A data frame containing 806 observations, representing visits to a
#' psoriatic arthritis (PsA) clinic from 305 patients.  The rows are grouped by
#' patient number and ordered by examination time. Each row represents an
#' examination and contains additional covariates.
#' 
#' \tabular{rll}{ \code{ptnum} \tab (numeric) \tab Patient identification
#' number \cr \code{months} \tab (numeric) \tab Examination time in months \cr
#' \code{state} \tab (numeric) \tab Clinical state of PsA.  Patients in states
#' 1, 2, 3 and 4 \cr \tab \tab have 0, 1 to 4, 5 to 9 and 10 or more damaged
#' joints, \cr \tab \tab respectively.  \cr \code{hieffusn} \tab (numeric) \tab
#' Presence of five or more effusions \cr \code{ollwsdrt} \tab (character) \tab
#' Erythrocyte sedimentation rate of less than 15 mm/h \cr }
#' @references Gladman, D. D. and Farewell, V.T. (1999) Progression in
#' psoriatic arthritis: role of time-varying clinical indicators.  J.
#' Rheumatol. 26(11):2409-13
#' @keywords datasets
#' @examples
#' 
#' ## Four-state progression-only model with high effusion and low
#' ## sedimentation rate as covariates on the progression rates.  High
#' ## effusion is assumed to have the same effect on the 1-2, 2-3, and 3-4
#' ## progression rates, while low sedimentation rate has the same effect
#' ## on the 1-2 and 2-3 intensities, but a different effect on the 3-4. 
#' 
#' data(psor)
#' psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
#' psor.msm <- msm(state ~ months, subject=ptnum, data=psor, 
#'                 qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
#'                 constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)),
#'                 fixedpars=FALSE, control = list(REPORT=1,trace=2), method="BFGS")
#' qmatrix.msm(psor.msm)
#' sojourn.msm(psor.msm)
#' hazard.msm(psor.msm)
#' 
NULL
