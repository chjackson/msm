\name{transient.msm}
\alias{transient.msm}
\alias{absorbing.msm}
\title{Transient and absorbing states}
\description{
  Returns the transient and absorbing states of either a fitted model or a
  transition intensity matrix. 
}
\usage{
transient.msm(x=NULL, qmatrix=NULL)
absorbing.msm(x=NULL, qmatrix=NULL)
}
\arguments{
  \item{x}{A fitted multi-state model as returned by \code{\link{msm}}.}
  \item{qmatrix}{A transition intensity matrix. The diagonal is ignored
    and taken to be minus the sum of the rest of the row.}
}
\value{
  A vector of the ordinal indices of the transient or absorbing states.  
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
