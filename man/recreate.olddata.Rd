\name{recreate.olddata}
\alias{recreate.olddata}
\title{Convert data stored in msm object to old format}
\description{
  Converts the \code{data} element of msm objects to the old format.
}
\usage{
recreate.olddata(x)
}
\arguments{
  \item{x}{Object returned by the \code{\link{msm}} function,
    representing a fitted multi-state model.} 
}
\value{
  A list of vectors and matrices in the undocumented ad-hoc format used
  for the \code{data} component of \code{msm} objects in \pkg{msm}
  versions 1.3.1 and earlier.  
}
\details{
  This is just provided for convenience and to illustrate the changes.
  It is not guaranteed to be complete, and is liable to be withdrawn.
  Users who were relying on the previous undocumented format are advised
  to upgrade their code to use the new format, which uses model frames
  and model design matrices in the standard format used in version 1.4,
  based on \code{\link{model.frame}} and \code{\link{model.matrix}}.
}

