\name{msm.summary}
\title{Summarise a fitted multi-state model}
\alias{summary.msm}
\alias{print.summary.msm}
\description{
  Summary method for fitted \code{\link{msm}} models. This
  is simply a wrapper around \code{\link{prevalence.msm}} which 
  produces a table of observed and
  expected state prevalences for each time, and for models with
  covariates, \code{\link{hazard.msm}} to print
  hazard ratios with 95\% confidence intervals for covariate effects.
}
\usage{
\method{summary}{msm}(object, hazard.scale=1, \dots)
}
\arguments{
  \item{object}{A fitted multi-state model object, as returned by
    \code{\link{msm}}.}
  \item{hazard.scale}{Vector with same elements as number of covariates
    on transition rates. Corresponds to the increase in each covariate
    used to calculate its hazard ratio. Defaults to all 1.}  
  \item{...}{Further arguments passed to \code{\link{prevalence.msm}}. }
}
\value{
  A list of class \code{summary.msm}, with components:
  
  \item{prevalences}{Output from \code{\link{prevalence.msm}}.}
  
  \item{hazard}{Output from \code{\link{hazard.msm}}.}
  
  \item{hazard.scale}{Value of the \code{hazard.scale} argument.}
}
\seealso{
  \code{\link{msm}},\code{\link{prevalence.msm}}, \code{\link{hazard.msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
