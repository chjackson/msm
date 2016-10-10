\name{model.frame.msm}
\alias{model.frame.msm}
\alias{model.matrix.msm}
\title{Extract original data from \code{msm} objects.}
\description{
Extract the data from a multi-state model fitted with \code{msm}.
}
\usage{
\method{model.frame}{msm}(formula, agg=FALSE, \dots)
\method{model.matrix}{msm}(object, model="intens", state=1, \dots)
}
\arguments{
  \item{formula}{A fitted multi-state model object, as returned by
    \code{\link{msm}}.}
  \item{agg}{Return the model frame in the efficient aggregated form
    used to calculate the likelihood internally for non-hidden Markov
    models.  This has one row for each unique combination of from-state,
    to-state, time lag, covariate value and observation type.  The
    variable named \code{"(nocc)"} counts how many observations of that
    combination there are in the original data. }
  \item{object}{A fitted multi-state model object, as returned by
    \code{\link{msm}}.}
  \item{model}{\code{"intens"} to return the design matrix for covariates on
    intensities, \code{"misc"} for misclassification probabilities, \code{"hmm"} for a
    general hidden Markov model, and \code{"inits"} for initial state
    probabilities in hidden Markov models.}
  \item{state}{State corresponding to the required covariate design
    matrix in a hidden Markov model.}
  \item{...}{Further arguments (not used).}
}
\value{
  \code{model.frame} returns a data frame with all the original
  variables used for the model fit, with any missing data removed (see
  \code{na.action} in \code{\link{msm}}).  The state, time, subject,
  \code{obstype} and \code{obstrue} variables are named
  \code{"(state)"}, \code{"(time)"}, \code{"(subject)"},
  \code{"(obstype)"} and \code{"(obstrue)"} respectively (note the
  brackets).  A variable called \code{"(obs)"} is the observation number
  from the original data before any missing data were dropped.
  The variable \code{"(pcomb)"} is used for computing the likelihood for hidden
  Markov models, and identifies
  which distinct time difference, \code{obstype} and covariate values
  (thus which distinct interval transition probability matrix) each observation
  corresponds to. 

  The model frame object has some other useful attributes, including
  \code{"usernames"} giving the user's original names for these
  variables (used for model refitting, e.g. in bootstrapping or cross
  validation) and \code{"covnames"} identifying which ones are covariates.  
  
  \code{model.matrix} returns a design matrix for a part of the
  model that includes covariates.  The required part is indicated by the
  \code{"model"} argument.

  For time-inhomogeneous models fitted with \code{"pci"}, these datasets
  will have imputed observations at each time change point, indicated
  where the variable \code{"(pci.imp)"} in the model frame is 1.  The
  model matrix for intensities will have factor contrasts for the
  \code{timeperiod} covariate.
}
\seealso{
  \code{\link{msm}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
