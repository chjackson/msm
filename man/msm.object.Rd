\name{msm.object}
\alias{msm.object}
\title{Fitted msm model objects}
\description{
  The \code{\link{msm}} function returns a list with the following
  components.  These are intended for developers and confident
  users.  To extract results from fitted model objects, functions such
  as \code{\link{qmatrix.msm}} or \code{\link{print.msm}} should be used
  instead.
}
\value{
\item{call}{The original call to \code{\link{msm}}, as returned by \code{\link{match.call}}.}
    \item{Qmatrices}{A list of matrices. The first component, labelled
    \code{logbaseline}, is a matrix containing the estimated
    transition intensities on the log scale with any covariates fixed at
    their means in the data (or at zero, if \code{center=FALSE}). The
    component labelled \code{baseline} is the equivalent on the
    untransformed scale. Each
    remaining component is a matrix giving the linear
    effects of the labelled covariate on the matrix of log
    intensities.  To extract an estimated intensity matrix on the natural
    scale, at an arbitrary combination of covariate values, use the
    function \code{\link{qmatrix.msm}}.
  }
  \item{QmatricesSE}{The standard error matrices corresponding to
    \code{Qmatrices}.
  }
  \item{QmatricesL,QmatricesU}{Corresponding lower and upper symmetric
    confidence limits, of width 0.95 unless specified otherwise by the
    \code{cl} argument.
  }
  \item{Ematrices}{A list of matrices. The first component, labelled
    \code{logitbaseline}, is the estimated misclassification probability
    matrix (expressed as as log odds relative to the probability of the
    true state) with any covariates fixed at their means in the data (or
    at zero, if \code{center=FALSE}). The
    component labelled \code{baseline} is the equivalent on the
    untransformed scale. Each
    remaining component is a matrix giving the linear
    effects of the labelled covariate on the matrix of logit
    misclassification probabilities.   To extract an estimated misclassification
    probability matrix on the natural scale, at an arbitrary combination
    of covariate values, use the function \code{\link{ematrix.msm}}.}
  \item{EmatricesSE}{The standard error matrices corresponding to \code{Ematrices}.}
  \item{EmatricesL,EmatricesU}{Corresponding lower and upper symmetric
    confidence limits, of width 0.95 unless specified otherwise by the
    \code{cl} argument.
  }

  \item{minus2loglik}{Minus twice the maximised log-likelihood.}
  \item{deriv}{Derivatives of the minus twice log-likelihood at its maximum.}

    \item{estimates}{Vector of untransformed maximum likelihood estimates
    returned from \code{\link{optim}}.  Transition intensities are on
    the log scale and misclassification probabilities are given as log
    odds relative to the probability of the true state.}
  \item{estimates.t}{Vector of transformed maximum likelihood estimates
    with intensities and probabilities on their natural scales.}

  \item{fixedpars}{Indices of \code{estimates} which were fixed during
    the maximum likelihood estimation.}
  \item{center}{Indicator for whether the estimation was performed with
    covariates centered on their means in the data.}
  \item{covmat}{Covariance matrix corresponding to \code{estimates}.}
  \item{ci}{Matrix of confidence intervals corresponding to \code{estimates.t}}

  \item{opt}{Return value from the optimisation routine (such as
    \code{\link{optim}} or \code{\link{nlm}}), giving information about
    the results of the optimisation.}
  \item{foundse}{Logical value indicating whether the Hessian was positive-definite
    at the supposed maximum of the likelihood.  If not, the covariance matrix of the
    parameters is unavailable. In these cases the optimisation has
    probably not converged to a maximum.
  }
  \item{data}{A list giving the data used for the model fit, for use in
    post-processing.  To extract it, use the methods
    \code{\link{model.frame.msm}} or \code{\link{model.matrix.msm}}.

    The format of this element changed in version
    1.4 of \pkg{msm}, so that it now contains a
    \code{\link{model.frame}} object \code{mf} with all the variables
    used in the model.  The previous format (an ad-hoc list of vectors and matrices) can be
    obtained with the function \code{recreate.olddata(msmobject)}, where
    \code{msmobject} is the object returned by \code{msm}.
  }
\item{qmodel}{A list of objects representing the transition matrix
  structure and options for likelihood calculation.  See
  \code{\link{qmodel.object}} for documentation of the components.}
\item{emodel}{A list of objects representing the misclassification model
structure, for models specified using the \code{ematrix} argument to
\code{\link{msm}}. See \code{\link{emodel.object}}.}
\item{qcmodel}{A list of objects representing the model for covariates
    on transition intensities.  See \code{\link{qcmodel.object}}.}
\item{ecmodel}{A list of objects representing the model for covariates
    on transition intensities.  See \code{\link{ecmodel.object}}.}
\item{hmodel}{A list of objects representing the hidden Markov model
    structure. See \code{\link{hmodel.object}}.}
\item{cmodel}{A list giving information about censored states.  See
  \code{\link{cmodel.object}}. }
\item{pci}{Cut points for time-varying intensities, as supplied to
  \code{\link{msm}}, but excluding any that are outside the times
  observed in the data.}
\item{paramdata}{A list giving information about the parameters of the multi-state
  model.  See \code{\link{paramdata.object}}.}
\item{cl}{Confidence interval width, as supplied to \code{\link{msm}}.}
\item{covariates}{Formula for covariates on intensities, as supplied to \code{\link{msm}}.}
\item{misccovariates}{Formula for covariates on misclassification probabilities, as supplied to \code{\link{msm}}.}
\item{hcovariates}{Formula for covariates on hidden Markov model outcomes, as supplied to \code{\link{msm}}.}
\item{initcovariates}{Formula for covariates on initial state occupancy
  probabilities in hidden Markov models, as supplied to
  \code{\link{msm}}.}
\item{sojourn}{
    A list as returned by \code{\link{sojourn.msm}}, with components:

    \code{mean} = estimated mean sojourn times in the transient states,
    with covariates fixed at their means (if center=TRUE) or at zero
    (if center=FALSE).

    \code{se} = corresponding standard errors.
}
}
