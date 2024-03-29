% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals_doc.R
\name{ecmodel.object}
\alias{ecmodel.object}
\title{Developer documentation: model for covariates on misclassification
probabilities}
\value{
\item{npars}{Number of covariate effect parameters.  This is defined
as the number of covariates on misclassification (with factors expanded as
contrasts) multiplied by the number of allowed misclassifications in the
model.  } \item{ndpars}{Number of distinct covariate effect parameters, as
\code{npars}, but after any equality constraints have been applied.}
\item{ncovs}{Number of covariates on misclassification, with factors
expanded as contrasts.} \item{constr}{List of equality constraints on these
covariate effects, as supplied in the \code{miscconstraint} argument to
\code{\link{msm}}.} \item{covlabels}{Names / labels of these covariates in
the model matrix (see \code{\link{model.matrix.msm}}).} \item{inits}{Initial
values for these covariate effects, as a vector formed from the
\code{misccovinits} list supplied to \code{\link{msm}}.}
\item{covmeans}{Means of these covariates in the data (excluding data not
required to fit the model, such as observations with missing data in other
elements or subjects' last observations).  This includes means of 0/1 factor
contrasts as well as continuous covariates (for historic reasons, which may
not be sensible).}
}
\description{
A list representing the model for covariates on misclassification
probabilities.
}
\seealso{
\code{\link{msm.object}}.
}
