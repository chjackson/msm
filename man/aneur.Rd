\name{aneur}
\docType{data}
\alias{aneur}
\title{Aortic aneurysm progression data}
\description{
  This dataset contains longitudinal measurements of grades of aortic
  aneurysms,  measured by ultrasound examination of the
  diameter of the aorta. 
  
}
\usage{aneur}

\format{
  A data frame containing 4337 rows, with each row corresponding to an
  ultrasound scan from one of 838 men over 65 years of age.

  \tabular{rll}{
    \code{ptnum}   \tab (numeric) \tab Patient identification number \cr
    \code{age}     \tab (numeric) \tab Recipient age at examination (years) \cr
    \code{diam}   \tab (numeric) \tab Aortic diameter\cr
    \code{state} \tab (numeric) \tab State of aneurysm. \cr
  }

  The states represent successive degrees of aneurysm severity, as
  indicated by the aortic diameter. 

  \tabular{rll}{
    State 1 \tab Aneurysm-free \tab < 30 cm \cr
    State 2 \tab Mild aneurysm \tab 30-44 cm \cr
    State 3 \tab Moderate aneurysm \tab 45-54 cm \cr
    State 4 \tab Severe aneurysm \tab > 55 cm \cr
  }

  683 of these men were aneurysm-free at age 65 and were re-screened every two
  years.  The remaining men were aneurysmal at entry and had successive
  screens with frequency depending on the state of the aneurysm.
  Severe aneurysms are repaired by surgery.
}
\source{
  The Chichester, U.K. randomised controlled trial of screening for
  abdominal aortic aneurysms by ultrasonography.
}
\references{
  Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W. and
  Couto, E.  Multi-state Markov models for disease progression with
  classification error. \emph{The Statistician}, 52(2): 193--209 (2003)

  Couto, E. and Duffy, S. W. and Ashton, H. A. and Walker, N. M.
  and Myles, J. P. and Scott, R. A. P. and Thompson, S. G. (2002)
  \emph{Probabilities of progression of aortic aneurysms: estimates
    and implications for screening policy}
  Journal of Medical Screening 9(1):40--42
}
\keyword{datasets}
