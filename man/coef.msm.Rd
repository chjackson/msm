\name{coef.msm}
\alias{coef.msm}
\title{Extract model coefficients}
\description{
  Extract the estimated log transition intensities and the corresponding
  linear effects of each covariate. 
}
\usage{
\method{coef}{msm}(object, \dots)
}
\arguments{
  \item{object}{A fitted multi-state model object, as returned by
    \code{\link{msm}}.}
  \item{...}{(unused) further arguments passed to or from other methods.}
}
\value{
  If there is no misclassification, \code{coef.msm} returns a list of matrices.
  The first component, labelled \code{logbaseline}, is a matrix containing the estimated
  transition intensities on the log scale with any covariates fixed at
  their means in the data. Each remaining component is a matrix giving the linear
  effects of the labelled covariate on the matrix of log
  intensities. \cr

  
  For misclassification models, \code{coef.msm} returns a list of
  lists. The first component, \code{Qmatrices}, is a list of matrices as described in the
  previous paragraph.  The additional component \code{Ematrices} is a
  list of similar format containing the logit-misclassification
  probabilities and any estimated covariate effects. 
}
\seealso{
  \code{\link{msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
