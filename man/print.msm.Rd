\name{print.msm}
\title{Print a fitted msm model object}
\alias{print.msm}
\alias{printnew.msm}
\description{
  Print a fitted msm model object
}
\usage{
\method{print}{msm}(x,covariates=NULL, digits=4, \dots)
printnew.msm(x, covariates=NULL, digits=4, \dots)
}
\arguments{

  \item{x}{Output from \code{\link{msm}}, representing a fitted
    multi-state model object.}

  \item{covariates}{Covariates for which to print ``baseline''
  transition intensities or misclassification probabilities. See
  \code{\link{qmatrix.msm}} for more details. }

  \item{digits}{Minimum number of significant digits, passed to
    \code{\link{format}}. Defaults to 4.}

  \item{...}{Other arguments to be passed to \code{\link{format}}.}
}
\details{
  This is the new method of formatting msm objects for printing.
  The old method was based on printing lists of matrices.
That produced a lot of wasted space for parameters which were zero,
and it was difficult to match corresponding numbers between matrices.
The new method presents all the transition intensities and covariate
effects as a single compact table, and likewise for misclassification
matrices.

Also in the old method, covariate effects were presented as log hazard ratios or
log odds ratios.  The log scale is more convenient mathematically, but
unnatural to interpret.  The new method presents hazard
ratios for covariates on transition intensities and odds ratios for
misclassification probabilities.

\code{printnew.msm} is an alias for \code{print.msm}.
  }
\value{
The object returned by \code{print.msm} is a numeric matrix with one column
for each estimate or confidence limit for intensities and their
covariates, in the same arrangement as
printed, but with the underlying numbers in full precision.  The results
formatted for printing are stored in the \code{"formatted"} attribute of
the object, as a character matrix.   These can alternatively be produced
by \code{\link{msm.form.qoutput}}, which has no printing side-effect.
\code{\link{msm.form.eoutput}} produces the same arrangement for
misclassification probabilities instead of intensities.
}
\seealso{
  \code{\link{msm}}, \code{\link{printold.msm}}, \code{\link{msm.form.qoutput}}.
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
