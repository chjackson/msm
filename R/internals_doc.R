

#' Developer documentation: censoring model object
#' 
#' A list giving information about censored states, their labels in the data
#' and what true states they represent.
#' 
#' @name cmodel.object
#' @return \item{ncens}{The number of distinct values used for censored
#' observations in the \code{state} data supplied to \code{\link{msm}}.}
#' \item{censor}{A vector of length \code{ncens}, giving the labels used for
#' censored states in the data.} \item{states}{A vector obtained by
#' \code{unlist()}ing a list with \code{ncens} elements, each giving the set of
#' true states that an observation with this label could be.}
#' \item{index}{Index into \code{states} for the first state corresponding to
#' each \code{censor}, plus an extra \code{length(states)+1}.}
#' @seealso \code{\link{msm.object}}.
NULL





#' Developer documentation: model for covariates on misclassification
#' probabilities
#' 
#' A list representing the model for covariates on misclassification
#' probabilities.
#' 
#' @name ecmodel.object
#' @return \item{npars}{Number of covariate effect parameters.  This is defined
#' as the number of covariates on misclassification (with factors expanded as
#' contrasts) multiplied by the number of allowed misclassifications in the
#' model.  } \item{ndpars}{Number of distinct covariate effect parameters, as
#' \code{npars}, but after any equality constraints have been applied.}
#' \item{ncovs}{Number of covariates on misclassification, with factors
#' expanded as contrasts.} \item{constr}{List of equality constraints on these
#' covariate effects, as supplied in the \code{miscconstraint} argument to
#' \code{\link{msm}}.} \item{covlabels}{Names / labels of these covariates in
#' the model matrix (see \code{\link{model.matrix.msm}}).} \item{inits}{Initial
#' values for these covariate effects, as a vector formed from the
#' \code{misccovinits} list supplied to \code{\link{msm}}.}
#' \item{covmeans}{Means of these covariates in the data (excluding data not
#' required to fit the model, such as observations with missing data in other
#' elements or subjects' last observations).  This includes means of 0/1 factor
#' contrasts as well as continuous covariates (for historic reasons, which may
#' not be sensible).}
#' @seealso \code{\link{msm.object}}.
NULL





#' Developer documentation: misclassification model structure object
#' 
#' A list giving information about the misclassifications assumed in a
#' multi-state model fitted with the \code{ematrix} argument of
#' \code{\link{msm}}.  Returned in a fitted \code{\link{msm}} model object.
#' This information is converted internally to a \code{hmodel} object (see
#' \code{\link{hmodel.object}}) for use in likelihood computations.
#' 
#' @name emodel.object
#' @return \item{nstates}{Number of states (same as \code{qmodel$nstates}).}
#' \item{npars}{Number of allowed misclassifications, equal to
#' \code{sum(imatrix)}.} \item{imatrix}{Indicator matrix for allowed
#' misclassifications.  This has \eqn{(r,s)} entry 1 if misclassification of
#' true state \eqn{r} as observed state \eqn{s} is possible.  diagonal entries
#' are arbitrarily set to 0.} \item{ematrix}{Matrix of initial values for the
#' misclassification probabilities, supplied as the \code{ematrix} argument of
#' \code{\link{msm}}.} \item{inits}{Vector of these initial values, reading
#' across rows of \code{qmatrix} and excluding the diagonal and disallowed
#' transitions.} \item{constr}{Indicators for equality constraints on baseline
#' misclassification probabilities, taken from the \code{econstraint} argument
#' to \code{\link{msm}}, and mapped if necessary to the set (1,2,3,...)}
#' \item{ndpars}{Number of distinct misclassification probabilities, after
#' applying equality constraints.} \item{nipars}{Number of initial state
#' occupancy probabilities being estimated.  This is zero if
#' \code{est.initprobs=FALSE}, otherwise equal to the number of states.}
#' \item{initprobs}{Initial state occupancy probabilities, as supplied to
#' \code{\link{msm}} (initial values before estimation, if
#' \code{est.initprobs=TRUE}.)} \item{est.initprobs}{Are initial state
#' occupancy probabilities estimated (\code{TRUE} or \code{FALSE}), as supplied
#' in the \code{est.initprobs} argument of \code{\link{msm}}.}
#' @seealso \code{\link{msm.object}},\code{\link{qmodel.object}},
#' \code{\link{hmodel.object}}.
NULL


#' Developer documentation: hidden Markov model structure object
#' 
#' A list giving information about the models for the outcome data
#' conditionally on the states of a hidden Markov model.  Used in internal
#' computations, and returned in a fitted \code{\link{msm}} model object.
#' 
#' @name hmodel.object
#' @return \item{hidden}{\code{TRUE} for hidden Markov models, \code{FALSE}
#' otherwise.} \item{nstates}{Number of states, the same as
#' \code{qmodel$nstates}.} \item{fitted}{\code{TRUE} if the parameter values in
#' \code{pars} are the maximum likelihood estimates, \code{FALSE} if they are
#' the initial values.} \item{models}{The outcome distribution for each hidden
#' state.  A vector of length \code{nstates} whose \eqn{r}th entry is the index
#' of the state \eqn{r} outcome distributions in the vector of supported
#' distributions.  The vector of supported distributions is given in full by
#' \code{msm:::.msm.HMODELS}: the first few are 1 for categorical outcome, 2
#' for identity, 3 for uniform and 4 for normal. } \item{labels}{String
#' identifying each distribution in \code{models}.} \item{npars}{Vector of
#' length \code{nstates} giving the number of parameters in each outcome
#' distribution, excluding covariate effects.} \item{nipars}{Number of initial
#' state occupancy probabilities being estimated.  This is zero if
#' \code{est.initprobs=FALSE}, otherwise equal to the number of states.}
#' \item{totpars}{Total number of parameters, equal to \code{sum(npars)}. }
#' \item{pars}{A vector of length \code{totpars}, made from concatenating a
#' list of length \code{nstates} whose \eqn{r}th component is vector of the
#' parameters for the state \eqn{r} outcome distribution.  } \item{plabs}{List
#' with the names of the parameters in \code{pars}.} \item{parstate}{A vector
#' of length \code{totpars}, whose \eqn{i}th element is the state corresponding
#' to the \eqn{i}th parameter.} \item{firstpar}{A vector of length
#' \code{nstates} giving the index in \code{pars} of the first parameter for
#' each state.} \item{locpars}{Index in \code{pars} of parameters which can
#' have covariates on them. } \item{initprobs}{Initial state occupancy
#' probabilities, as supplied to \code{\link{msm}} (initial values before
#' estimation, if \code{est.initprobs=TRUE}.)} \item{est.initprobs}{Are initial
#' state occupancy probabilities estimated (\code{TRUE} or \code{FALSE}), as
#' supplied in the \code{est.initprobs} argument of \code{\link{msm}}.}
#' \item{ncovs}{Number of covariate effects per parameter in \code{pars}, with,
#' e.g. factor contrasts expanded.} \item{coveffect}{Vector of covariate
#' effects, of length \code{sum(ncovs)}.} \item{covlabels}{Labels of these
#' effects.} \item{coveffstate}{Vector indicating state corresponding to each
#' element of \code{coveffect}.} \item{ncoveffs}{Number of covariate effects on
#' HMM outcomes, equal to \code{sum(ncovs)}.} \item{nicovs}{Vector of length
#' \code{nstates-1} giving the number of covariate effects on each initial
#' state occupancy probability (log relative to the baseline probability).}
#' \item{icoveffect}{Vector of length \code{sum(nicovs)} giving covariate
#' effects on initial state occupancy probabilities.} \item{nicoveffs}{Number
#' of covariate effects on initial state occupancy probabilities, equal to
#' \code{sum(nicovs)}.} \item{constr}{Constraints on (baseline) hidden Markov
#' model outcome parameters, as supplied in the \code{hconstraint} argument of
#' \code{\link{msm}}, excluding covariate effects, converted to a vector and
#' mapped to the set 1,2,3,\ldots{} if necessary.} \item{covconstr}{Vector of
#' constraints on covariate effects in hidden Markov outcome models, as
#' supplied in the \code{hconstraint} argument of \code{\link{msm}}, excluding
#' baseline parameters, converted to a vector and mapped to the set
#' 1,2,3,\ldots{} if necessary.} \item{ranges}{Matrix of range restrictions for
#' HMM parameters, including those given to the \code{hranges} argument to
#' \code{\link{msm}}.} \item{foundse}{\code{TRUE} if standard errors are
#' available for the estimates.} \item{initpmat}{Matrix of initial state
#' occupancy probabilities with one row for each subject (estimated if
#' \code{est.initprobs=TRUE}).} \item{ci}{Confidence intervals for baseline HMM
#' outcome parameters.} \item{covci}{Confidence intervals for covariate effects
#' in HMM outcome models.}
#' @seealso \code{\link{msm.object}},\code{\link{qmodel.object}},
#' \code{\link{emodel.object}}.
NULL



#' Fitted msm model objects
#' 
#' The \code{\link{msm}} function returns a list with the following components.
#' These are intended for developers and confident users.  To extract results
#' from fitted model objects, functions such as \code{\link{qmatrix.msm}} or
#' \code{\link{print.msm}} should be used instead.
#' 
#' @name msm.object
#' @return \item{call}{The original call to \code{\link{msm}}, as returned by
#' \code{\link{match.call}}.} \item{Qmatrices}{A list of matrices. The first
#' component, labelled \code{logbaseline}, is a matrix containing the estimated
#' transition intensities on the log scale with any covariates fixed at their
#' means in the data (or at zero, if \code{center=FALSE}). The component
#' labelled \code{baseline} is the equivalent on the untransformed scale. Each
#' remaining component is a matrix giving the linear effects of the labelled
#' covariate on the matrix of log intensities.  To extract an estimated
#' intensity matrix on the natural scale, at an arbitrary combination of
#' covariate values, use the function \code{\link{qmatrix.msm}}.  }
#' \item{QmatricesSE}{The standard error matrices corresponding to
#' \code{Qmatrices}.  } \item{QmatricesL,QmatricesU}{Corresponding lower and
#' upper symmetric confidence limits, of width 0.95 unless specified otherwise
#' by the \code{cl} argument.  } \item{Ematrices}{A list of matrices. The first
#' component, labelled \code{logitbaseline}, is the estimated misclassification
#' probability matrix (expressed as as log odds relative to the probability of
#' the true state) with any covariates fixed at their means in the data (or at
#' zero, if \code{center=FALSE}). The component labelled \code{baseline} is the
#' equivalent on the untransformed scale. Each remaining component is a matrix
#' giving the linear effects of the labelled covariate on the matrix of logit
#' misclassification probabilities.  To extract an estimated misclassification
#' probability matrix on the natural scale, at an arbitrary combination of
#' covariate values, use the function \code{\link{ematrix.msm}}.}
#' \item{EmatricesSE}{The standard error matrices corresponding to
#' \code{Ematrices}.} \item{EmatricesL,EmatricesU}{Corresponding lower and
#' upper symmetric confidence limits, of width 0.95 unless specified otherwise
#' by the \code{cl} argument.  }
#' 
#' \item{minus2loglik}{Minus twice the maximised log-likelihood.}
#' \item{deriv}{Derivatives of the minus twice log-likelihood at its maximum.}
#' 
#' \item{estimates}{Vector of untransformed maximum likelihood estimates
#' returned from \code{\link{optim}}.  Transition intensities are on the log
#' scale and misclassification probabilities are given as log odds relative to
#' the probability of the true state.} \item{estimates.t}{Vector of transformed
#' maximum likelihood estimates with intensities and probabilities on their
#' natural scales.}
#' 
#' \item{fixedpars}{Indices of \code{estimates} which were fixed during the
#' maximum likelihood estimation.} \item{center}{Indicator for whether the
#' estimation was performed with covariates centered on their means in the
#' data.} \item{covmat}{Covariance matrix corresponding to \code{estimates}.}
#' \item{ci}{Matrix of confidence intervals corresponding to
#' \code{estimates.t}}
#' 
#' \item{opt}{Return value from the optimisation routine (such as
#' \code{\link{optim}} or \code{\link{nlm}}), giving information about the
#' results of the optimisation.} \item{foundse}{Logical value indicating
#' whether the Hessian was positive-definite at the supposed maximum of the
#' likelihood.  If not, the covariance matrix of the parameters is unavailable.
#' In these cases the optimisation has probably not converged to a maximum.  }
#' \item{data}{A list giving the data used for the model fit, for use in
#' post-processing.  To extract it, use the methods
#' \code{\link{model.frame.msm}} or \code{\link{model.matrix.msm}}.
#' 
#' The format of this element changed in version 1.4 of \pkg{msm}, so that it
#' now contains a \code{\link{model.frame}} object \code{mf} with all the
#' variables used in the model.  The previous format (an ad-hoc list of vectors
#' and matrices) can be obtained with the function
#' \code{recreate.olddata(msmobject)}, where \code{msmobject} is the object
#' returned by \code{msm}.  } \item{qmodel}{A list of objects representing the
#' transition matrix structure and options for likelihood calculation.  See
#' \code{\link{qmodel.object}} for documentation of the components.}
#' \item{emodel}{A list of objects representing the misclassification model
#' structure, for models specified using the \code{ematrix} argument to
#' \code{\link{msm}}. See \code{\link{emodel.object}}.} \item{qcmodel}{A list
#' of objects representing the model for covariates on transition intensities.
#' See \code{\link{qcmodel.object}}.} \item{ecmodel}{A list of objects
#' representing the model for covariates on transition intensities.  See
#' \code{\link{ecmodel.object}}.} \item{hmodel}{A list of objects representing
#' the hidden Markov model structure. See \code{\link{hmodel.object}}.}
#' \item{cmodel}{A list giving information about censored states.  See
#' \code{\link{cmodel.object}}. } \item{pci}{Cut points for time-varying
#' intensities, as supplied to \code{\link{msm}}, but excluding any that are
#' outside the times observed in the data.} \item{paramdata}{A list giving
#' information about the parameters of the multi-state model.  See
#' \code{\link{paramdata.object}}.} \item{cl}{Confidence interval width, as
#' supplied to \code{\link{msm}}.} \item{covariates}{Formula for covariates on
#' intensities, as supplied to \code{\link{msm}}.}
#' \item{misccovariates}{Formula for covariates on misclassification
#' probabilities, as supplied to \code{\link{msm}}.} \item{hcovariates}{Formula
#' for covariates on hidden Markov model outcomes, as supplied to
#' \code{\link{msm}}.} \item{initcovariates}{Formula for covariates on initial
#' state occupancy probabilities in hidden Markov models, as supplied to
#' \code{\link{msm}}.} \item{sojourn}{ A list as returned by
#' \code{\link{sojourn.msm}}, with components:
#' 
#' \code{mean} = estimated mean sojourn times in the transient states, with
#' covariates fixed at their means (if center=TRUE) or at zero (if
#' center=FALSE).
#' 
#' \code{se} = corresponding standard errors. }
NULL





#' Developer documentation: internal msm parameters object
#' 
#' An object giving information about the parameters of the multi-state model.
#' Used internally during maximum likelihood estimation and arranging results.
#' Returned as the \code{paramdata} component of a fitted \code{\link{msm}}
#' model object.
#' 
#' @name paramdata.object
#' @return \item{inits}{Vector of initial values for distinct parameters which
#' are being estimated.  These have been transformed to the real line (e.g. by
#' log), and exclude parameters being fixed at their initial values, parameters
#' defined to be always fixed (e.g. binomial denominators) and parameters
#' constrained to equal previous ones.} \item{plabs}{Names of parameters in
#' \code{allinits}.} \item{allinits}{Vector of parameter values before
#' estimation, including those which are fixed or constrained to equal other
#' parameters, and transformed to the real line.} \item{hmmpars}{Indices of
#' \code{allinits} which represent baseline parameters of hidden Markov outcome
#' models (thus excluding covariate effects in HMMs and initial state occupancy
#' probabilities). } \item{fixed}{\code{TRUE} if all parameters are fixed,
#' \code{FALSE} otherwise.} \item{fixedpars}{Indices of parameters in
#' \code{allinits} which are fixed, either by definition or as requested by the
#' user in the \code{fixedpars} argument to \code{\link{msm}}. Excludes
#' parameters fixed by constraining to equal other parameters.}
#' \item{notfixed}{Indices of parameters which are not fixed by the definition
#' of \code{fixedpars}.} \item{optpars}{Indices of parameters in
#' \code{allinits} being estimated, thus those included in \code{inits}.}
#' \item{auxpars}{Indices of "auxiliary" parameters which are always fixed, for
#' example, binomial denominators (\code{\link{hmmBinom}}) and the \code{which}
#' parameter in \code{\link{hmmIdent}}.} \item{constr}{Vector of integers, of
#' length \code{npars}, indicating which sets of parameters are constrained to
#' be equal to each other.  If two of these integers are equal the
#' corresponding parameters are equal.  A negative element indicates that
#' parameter is defined to be minus some other parameter (this is used for
#' covariate effects on transition intensities).} \item{npars}{Total number of
#' parameters, equal to \code{length(allinits)}.} \item{nfix}{Number of fixed
#' parameters, equal to \code{length(fixedpars)}. } \item{nopt}{Number of
#' parameters being estimated, equal to \code{length(inits)} and
#' \code{length(optpars)}.} \item{ndup}{Number of parameters defined as
#' duplicates of previous parameters by equality constraints (currently
#' unused).} \item{ranges}{Matrix of defined ranges for each parameter on the
#' natural scale (e.g. 0 to infinity for rate parameters). } \item{opt}{Object
#' returned by the optimisation routine (such as \code{\link{optim}}).}
#' \item{foundse}{\code{TRUE} if standard errors are available after
#' optimisation. If \code{FALSE} the optimisation probably hasn't converged. }
#' \item{lik}{Minus twice the log likelihood at the parameter estimates.}
#' \item{deriv}{Derivatives of the minus twice log likelihood at the parameter
#' estimates, if available. } \item{information}{Corresponding expected
#' information matrix at the parameter estimates, if available.}
#' \item{params}{Vector of parameter values after maximum likelihood
#' estimation, corresponding to \code{allinits}, still on the real-line
#' transformed scale.} \item{covmat}{Covariance matrix corresponding to
#' \code{params}.} \item{ci}{Matrix of confidence intervals corresponding to
#' \code{params}, with nominal coverage (default 0.95) defined by the \code{cl}
#' argument of \code{\link{msm}}. } \item{estimates.t}{Vector of parameter
#' estimates, as \code{params} but with parameters on their natural scales.}
#' @seealso \code{\link{msm.object}}
NULL


#' Developer documentation: model for covariates on transition intensities
#' 
#' A list representing the model for covariates on transition intensities
#' 
#' @name qcmodel.object
#' @return \item{npars}{Number of covariate effect parameters.  This is defined
#' as the number of covariates on intensities (with factors expanded as
#' contrasts) multiplied by the number of allowed transitions in the model.
#' 
#' Note if \code{\link{msm}} was called with \code{covariates} set to a list of
#' different covariates for different intensities, then this will include
#' covariate effects that are implicitly defined as zero by this list.  The
#' information in \code{\link[=paramdata.object]{paramdata}} objects can be
#' used to identify wich ones are fixed at zero.
#' 
#' This also includes any \code{timeperiod} covariates in a time-inhomogeneous
#' model defined by the \code{pci} option to \code{\link{msm}}.  }
#' \item{ndpars}{Number of distinct covariate effect parameters, as
#' \code{npars}, but after any equality constraints have been applied.}
#' \item{ncovs}{Number of covariates on intensities, with factors expanded as
#' contrasts.} \item{constr}{List of equality constraints on these covariate
#' effects, as supplied in the \code{constraint} argument to
#' \code{\link{msm}}.} \item{covlabels}{Names / labels of these covariates in
#' the model matrix (see \code{\link{model.matrix.msm}}).} \item{inits}{Initial
#' values for these covariate effects, as a vector formed from the
#' \code{covinits} list supplied to \code{\link{msm}}.} \item{covmeans}{Means
#' of these covariates in the data (excluding data not required to fit the
#' model, such as observations with missing data in other elements or subjects'
#' last observations).  This includes means of 0/1 factor contrasts as well as
#' continuous covariates (for historic reasons, which may not be sensible).}
#' @seealso \code{\link{msm.object}}.
NULL





#' Developer documentation: transition model structure object
#' 
#' A list giving information about the structure of states and allowed
#' transitions in a multi-state model, and options for likelihood calculation.
#' Used in internal computations, and returned in a fitted \code{\link{msm}}
#' model object.
#' 
#' @name qmodel.object
#' @return \item{nstates}{Number of states} \item{iso}{Label for which basic
#' structure the model is isomorphic to in the list of structures for which
#' analytic formulae for the transition probabilities are implemented in the
#' source file \code{src/analyticp.c}.  This list is given by the internal
#' object \code{msm:::.msm.graphs} which is defined and documented in the
#' source file \code{R/constants.R}.
#' 
#' \code{iso} is 0 if the analytic P matrix is not implemented for this
#' structure, or if analytic P matrix calculations are disabled using
#' \code{use.analyticp=FALSE} in the call to \code{\link{msm}}.  }
#' \item{perm}{Permutation required to convert the base isomorphism into the
#' structure of this model. A vector of integers whose \eqn{r}th element is the
#' state number in the base structure representing state \eqn{r} in the current
#' structure.  } \item{qperm}{Inverse permutation: vector whose \eqn{r}th
#' element is the state number in the current structure representing the
#' \eqn{r}th state in the base structure.} \item{npars}{Number of allowed
#' instantaneous transitions, equal to \code{sum(imatrix)}.}
#' \item{imatrix}{Indicator matrix for allowed instantaneous transitions.  This
#' has \eqn{(r,s)} entry 1 if the transition from \eqn{r} to \eqn{s} is
#' permitted in continuous time, and 0 otherwise.  The diagonal entries are
#' arbitrarily set to 0.} \item{qmatrix}{Matrix of initial values for the
#' transition intensities, supplied as the \code{qmatrix} argument of
#' \code{\link{msm}}.} \item{inits}{Vector of these initial values, reading
#' across rows of \code{qmatrix} and excluding the diagonal and disallowed
#' transitions.} \item{constr}{Indicators for equality constraints on baseline
#' intensities, taken from the \code{qconstraint} argument to
#' \code{\link{msm}}, and mapped if necessary to the set (1,2,3,...).}
#' \item{ndpars}{Number of distinct allowed instantaneous transitions, after
#' applying equality constraints.} \item{expm}{Use \pkg{expm} package to
#' calculate matrix exponentials for likelihoods, as supplied to the
#' \code{use.expm} argument of \code{\link{msm}}.  \code{TRUE} or
#' \code{FALSE}.}
#' @seealso \code{\link{msm.object}},\code{\link{emodel.object}},
#' \code{\link{hmodel.object}}.
NULL


