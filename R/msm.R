#' Multi-state Markov and hidden Markov models in continuous time
#' 
#' Fit a continuous-time Markov or hidden Markov multi-state model by maximum
#' likelihood. Observations of the process can be made at arbitrary times, or
#' the exact times of transition between states can be known.  Covariates can
#' be fitted to the Markov chain transition intensities or to the hidden Markov
#' observation process.
#' 
#' For full details about the methodology behind the \pkg{msm} package, refer
#' to the PDF manual \file{msm-manual.pdf} in the \file{doc} subdirectory of
#' the package. This includes a tutorial in the typical use of \pkg{msm}.  The
#' paper by Jackson (2011) in Journal of Statistical Software presents the
#' material in this manual in a more concise form.
#' 
#' \pkg{msm} was designed for fitting \emph{continuous-time} Markov models,
#' processes where transitions can occur at any time.  These models are defined
#' by \emph{intensities}, which govern both the time spent in the current state
#' and the probabilities of the next state.  In \emph{discrete-time models},
#' transitions are known in advance to only occur at multiples of some time
#' unit, and the model is purely governed by the probability distributions of
#' the state at the next time point, conditionally on the state at the current
#' time. These can also be fitted in \pkg{msm}, assuming that there is a
#' continuous-time process underlying the data.  Then the fitted transition
#' probability matrix over one time period, as returned by
#' \code{pmatrix.msm(...,t=1)} is equivalent to the matrix that governs the
#' discrete-time model. However, these can be fitted more efficiently using
#' multinomial logistic regression, for example, using \code{multinom} from the
#' R package \pkg{nnet} (Venables and Ripley, 2002).
#' 
#' For simple continuous-time multi-state Markov models, the likelihood is
#' calculated in terms of the transition intensity matrix \eqn{Q}. When the
#' data consist of observations of the Markov process at arbitrary times, the
#' exact transition times are not known.  Then the likelihood is calculated
#' using the transition probability matrix \eqn{P(t) = \exp(tQ)}{P(t) =
#' exp(tQ)}, where \eqn{\exp}{exp} is the matrix exponential.  If state \eqn{i}
#' is observed at time \eqn{t} and state \eqn{j} is observed at time \eqn{u},
#' then the contribution to the likelihood from this pair of observations is
#' the \eqn{i,j} element of \eqn{P(u - t)}. See, for example, Kalbfleisch and
#' Lawless (1985), Kay (1986), or Gentleman \emph{et al.} (1994).
#' 
#' For hidden Markov models, the likelihood for an individual with \eqn{k}
#' observations is calculated directly by summing over the unknown state at
#' each time, producing a product of \eqn{k} matrices. The calculation is a
#' generalisation of the method described by Satten and Longini (1996), and
#' also by Jackson and Sharples (2002), and Jackson \emph{et al.} (2003).
#' 
#' There must be enough information in the data on each state to estimate each
#' transition rate, otherwise the likelihood will be flat and the maximum will
#' not be found.  It may be appropriate to reduce the number of states in the
#' model, the number of allowed transitions, or the number of covariate
#' effects, to ensure convergence.  Hidden Markov models, and situations where
#' the value of the process is only known at a series of snapshots, are
#' particularly susceptible to non-identifiability, especially when combined
#' with a complex transition matrix.  Choosing an appropriate set of initial
#' values for the optimisation can also be important.  For flat likelihoods,
#' 'informative' initial values will often be required.  See the PDF manual for
#' other tips.
#' 
#' @param formula A formula giving the vectors containing the observed states
#' and the corresponding observation times. For example,
#' 
#' \code{state ~ time}
#' 
#' Observed states should be numeric variables in the set \code{1, \dots{}, n},
#' where \code{n} is the number of states.  Factors are allowed only if their
#' levels are called \code{"1", \dots{}, "n"}.
#' 
#' The times can indicate different types of observation scheme, so be careful
#' to choose the correct \code{obstype}.
#' 
#' For hidden Markov models, \code{state} refers to the outcome variable, which
#' need not be a discrete state.  It may also be a matrix, giving multiple
#' observations at each time (see \code{\link{hmmMV}}).
#' @param subject Vector of subject identification numbers for the data
#' specified by \code{formula}. If missing, then all observations are assumed
#' to be on the same subject. These must be sorted so that all observations on
#' the same subject are adjacent.
#' @param data Optional data frame in which to interpret the variables supplied
#' in \code{formula}, \code{subject}, \code{covariates}, \code{misccovariates},
#' \code{hcovariates}, \code{obstype} and \code{obstrue}.
#' @param qmatrix Matrix which indicates the allowed transitions in the
#' continuous-time Markov chain, and optionally also the initial values of
#' those transitions. If an instantaneous transition is not allowed from state
#' \eqn{r} to state \eqn{s}, then \code{qmatrix} should have \eqn{(r,s)} entry
#' 0, otherwise it should be non-zero.
#' 
#' If supplying initial values yourself, then the non-zero entries should be
#' those values.  If using \code{gen.inits=TRUE} then the non-zero entries can
#' be anything you like (conventionally 1).  Any diagonal entry of
#' \code{qmatrix} is ignored, as it is constrained to be equal to minus the sum
#' of the rest of the row.
#' 
#' For example,\cr
#' 
#' \code{ rbind( c( 0, 0.1, 0.01 ), c( 0.1, 0, 0.2 ), c( 0, 0, 0 ) ) }\cr
#' 
#' represents a 'health - disease - death' model, with initial transition
#' intensities 0.1 from health to disease, 0.01 from health to death, 0.1 from
#' disease to health, and 0.2 from disease to death.
#' 
#' If the states represent ordered levels of severity of a disease, then this
#' matrix should usually only allow transitions between adjacent states.  For
#' example, if someone was observed in state 1 ("mild") at their first
#' observation, followed by state 3 ("severe") at their second observation,
#' they are assumed to have passed through state 2 ("moderate") in between, and
#' the 1,3 entry of \code{qmatrix} should be zero.
#' 
#' The initial intensities given here are with any covariates set to their
#' means in the data (or set to zero, if \code{center = FALSE}). If any
#' intensities are constrained to be equal using \code{qconstraint}, then the
#' initial value is taken from the first of these (reading across rows).
#' @param gen.inits If \code{TRUE}, then initial values for the transition
#' intensities are generated automatically using the method in
#' \code{\link{crudeinits.msm}}. The non-zero entries of the supplied
#' \code{qmatrix} are assumed to indicate the allowed transitions of the model.
#' This is not available for hidden Markov models, including models with
#' misclassified states.
#' @param ematrix If misclassification between states is to be modelled, this
#' should be a matrix of initial values for the misclassification
#' probabilities.  The rows represent underlying states, and the columns
#' represent observed states.  If an observation of state \eqn{s} is not
#' possible when the subject occupies underlying state \eqn{r}, then
#' \code{ematrix} should have \eqn{(r,s)} entry 0.  Otherwise \code{ematrix}
#' should have \eqn{(r,s)} entry corresponding to the probability of observing
#' \eqn{s} conditionally on occupying true state \eqn{r}. The diagonal of
#' \code{ematrix} is ignored, as rows are constrained to sum to 1.  For
#' example, \cr
#' 
#' \code{ rbind( c( 0, 0.1, 0 ), c( 0.1, 0, 0.1 ), c( 0, 0.1, 0 ) ) }\cr
#' 
#' represents a model in which misclassifications are only permitted between
#' adjacent states.
#' 
#' If any probabilities are constrained to be equal using \code{econstraint},
#' then the initial value is taken from the first of these (reading across
#' rows).
#' 
#' For an alternative way of specifying misclassification models, see
#' \code{hmodel}.
#' @param hmodel Specification of the hidden Markov model (HMM).  This should
#' be a list of return values from HMM constructor functions.  Each element of
#' the list corresponds to the outcome model conditionally on the corresponding
#' underlying state.  Univariate constructors are described in
#' the\code{\link{hmm-dists}} help page.  These may also be grouped together to
#' specify a multivariate HMM with a set of conditionally independent
#' univariate outcomes at each time, as described in \code{\link{hmmMV}}.
#' 
#' For example, consider a three-state hidden Markov model.  Suppose the
#' observations in underlying state 1 are generated from a Normal distribution
#' with mean 100 and standard deviation 16, while observations in underlying
#' state 2 are Normal with mean 54 and standard deviation 18. Observations in
#' state 3, representing death, are exactly observed, and coded as 999 in the
#' data.  This model is specified as
#' 
#' \code{hmodel = list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18),
#' hmmIdent(999))}
#' 
#' The mean and standard deviation parameters are estimated starting from these
#' initial values. If multiple parameters are constrained to be equal using
#' \code{hconstraint}, then the initial value is taken from the value given on
#' the first occasion that parameter appears in \code{hmodel}.
#' 
#' See the \code{\link{hmm-dists}} help page for details of the constructor
#' functions for each univariate distribution.
#' 
#' A misclassification model, that is, a hidden Markov model where the outcomes
#' are misclassified observations of the underlying states, can either be
#' specified using a list of \code{\link{hmmCat}} or \code{\link{hmmIdent}}
#' objects, or by using an \code{ematrix}.
#' 
#' For example, \cr
#' 
#' \code{ ematrix = rbind( c( 0, 0.1, 0, 0 ), c( 0.1, 0, 0.1, 0 ), c( 0, 0.1,
#' 0, 0), c( 0, 0, 0, 0) ) }\cr
#' 
#' is equivalent to \cr
#' 
#' \code{hmodel = list( hmmCat(prob=c(0.9, 0.1, 0, 0)), hmmCat(prob=c(0.1, 0.8,
#' 0.1, 0)), hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()) }\cr
#' 
#' @param obstype A vector specifying the observation scheme for each row of
#' the data. This can be included in the data frame \code{data} along with the
#' state, time, subject IDs and covariates.  Its elements should be either 1, 2
#' or 3, meaning as follows:
#' 
#' \describe{ \item{1}{An observation of the process at an arbitrary time (a
#' "snapshot" of the process, or "panel-observed" data). The states are unknown
#' between observation times.} \item{2}{An exact transition time, with the
#' state at the previous observation retained until the current observation.
#' An observation may represent a transition to a different state or a repeated
#' observation of the same state (e.g. at the end of follow-up).  Note that if
#' all transition times are known, more flexible models could be fitted with
#' packages other than \pkg{msm} - see the note under \code{exacttimes}.
#' 
#' Note also that if the previous state was censored using \code{censor}, for
#' example known only to be state 1 or state 2, then \code{obstype} 2 means
#' that either state 1 is retained or state 2 is retained until the current
#' observation - this does not allow for a change of state in the middle of the
#' observation interval.  } \item{3}{An exact transition time, but the state at
#' the instant before entering this state is unknown. A common example is death
#' times in studies of chronic diseases.} } If \code{obstype} is not specified,
#' this defaults to all 1. If \code{obstype} is a single number, all
#' observations are assumed to be of this type.  The obstype value for the
#' first observation from each subject is not used.
#' 
#' This is a generalisation of the \code{deathexact} and \code{exacttimes}
#' arguments to allow different schemes per observation.  \code{obstype}
#' overrides both \code{deathexact} and \code{exacttimes}.
#' 
#' \code{exacttimes=TRUE} specifies that all observations are of obstype 2.
#' 
#' \code{deathexact = death.states} specifies that all observations of
#' \code{death.states} are of type 3.  \code{deathexact = TRUE} specifies that
#' all observations in the final absorbing state are of type 3.
#' @param obstrue In misclassification models specified with \code{ematrix},
#' \code{obstrue} is a vector of logicals (\code{TRUE} or \code{FALSE}) or
#' numerics (1 or 0) specifying which observations (\code{TRUE}, 1) are
#' observations of the underlying state without error, and which (\code{FALSE},
#' 0) are realisations of a hidden Markov model.
#' 
#' In HMMs specified with \code{hmodel}, where the hidden state is known at
#' some times, if \code{obstrue} is supplied it is assumed to contain the
#' actual true state data.  Elements of \code{obstrue} at times when the hidden
#' state is unknown are set to \code{NA}.  This allows the information from HMM
#' outcomes generated conditionally on the known state to be included in the
#' model, thus improving the estimation of the HMM outcome distributions.
#' 
#' HMMs where the true state is known to be within a specific set at specific
#' times can be defined with a combination of \code{censor} and \code{obstrue}.
#' In these models, a code is defined for the \code{state} outcome (see
#' \code{censor}), and \code{obstrue} is set to 1 for observations where the
#' true state is known to be one of the elements of \code{censor.states} at the
#' corresponding time.
#' @param covariates A formula or a list of formulae representing the
#' covariates on the transition intensities via a log-linear model. If a single
#' formula is supplied, like
#' 
#' \code{covariates = ~ age + sex + treatment}
#' 
#' then these covariates are assumed to apply to all intensities.  If a named
#' list is supplied, then this defines a potentially different model for each
#' named intensity.  For example,
#' 
#' \code{covariates = list("1-2" = ~ age, "2-3" = ~ age + treatment)}
#' 
#' specifies an age effect on the state 1 - state 2 transition, additive age
#' and treatment effects on the state 2 - state 3 transition, but no covariates
#' on any other transitions that are allowed by the \code{qmatrix}.
#' 
#' If covariates are time dependent, they are assumed to be constant in between
#' the times they are observed, and the transition probability between a pair
#' of times \eqn{(t1, t2)} is assumed to depend on the covariate value at
#' \eqn{t1}.
#' @param covinits Initial values for log-linear effects of covariates on the
#' transition intensities. This should be a named list with each element
#' corresponding to a covariate.  A single element contains the initial values
#' for that covariate on each transition intensity, reading across the rows in
#' order.  For a pair of effects constrained to be equal, the initial value for
#' the first of the two effects is used.
#' 
#' For example, for a model with the above \code{qmatrix} and age and sex
#' covariates, the following initialises all covariate effects to zero apart
#' from the age effect on the 2-1 transition, and the sex effect on the 1-3
#' transition.  \code{ covinits = list(sex=c(0, 0, 0.1, 0), age=c(0, 0.1, 0,
#' 0))}
#' 
#' For factor covariates, name each level by concatenating the name of the
#' covariate with the level name, quoting if necessary. For example, for a
#' covariate \code{agegroup} with three levels \code{0-15, 15-60, 60-}, use
#' something like
#' 
#' \code{ covinits = list("agegroup15-60"=c(0, 0.1, 0, 0), "agegroup60-"=c(0.1,
#' 0.1, 0, 0))}
#' 
#' If not specified or wrongly specified, initial values are assumed to be
#' zero.
#' @param constraint A list of one numeric vector for each named covariate. The
#' vector indicates which covariate effects on intensities are constrained to
#' be equal. Take, for example, a model with five transition intensities and
#' two covariates. Specifying\cr
#' 
#' \code{constraint = list (age = c(1,1,1,2,2), treatment = c(1,2,3,4,5))}\cr
#' 
#' constrains the effect of age to be equal for the first three intensities,
#' and equal for the fourth and fifth. The effect of treatment is assumed to be
#' different for each intensity. Any vector of increasing numbers can be used
#' as indicators. The intensity parameters are assumed to be ordered by reading
#' across the rows of the transition matrix, starting at the first row,
#' ignoring the diagonals.
#' 
#' Negative elements of the vector can be used to indicate that particular
#' covariate effects are constrained to be equal to minus some other effects.
#' For example:
#' 
#' \code{constraint = list (age = c(-1,1,1,2,-2), treatment = c(1,2,3,4,5))
#' }\cr
#' 
#' constrains the second and third age effects to be equal, the first effect to
#' be minus the second, and the fifth age effect to be minus the fourth.  For
#' example, it may be realisitic that the effect of a covariate on the
#' "reverse" transition rate from state 2 to state 1 is minus the effect on the
#' "forward" transition rate, state 1 to state 2.  Note that it is not possible
#' to specify exactly which of the covariate effects are constrained to be
#' positive and which negative.  The maximum likelihood estimation chooses the
#' combination of signs which has the higher likelihood.
#' 
#' For categorical covariates, defined as factors, specify constraints as
#' follows:\cr
#' 
#' \code{list(..., covnameVALUE1 = c(...), covnameVALUE2 = c(...), ...)}\cr
#' 
#' where \code{covname} is the name of the factor, and \code{VALUE1},
#' \code{VALUE2}, ... are the labels of the factor levels (usually excluding
#' the baseline, if using the default contrasts).
#' 
#' Make sure the \code{contrasts} option is set appropriately, for example, the
#' default
#' 
#' \code{options(contrasts=c(contr.treatment, contr.poly))}
#' 
#' sets the first (baseline) level of unordered factors to zero, then the
#' baseline level is ignored in this specification.
#' 
#' To assume no covariate effect on a certain transition, use the
#' \code{fixedpars} argument to fix it at its initial value (which is zero by
#' default) during the optimisation.
#' 
#' @param misccovariates A formula representing the covariates on the
#' misclassification probabilities, analogously to \code{covariates}, via
#' multinomial logistic regression. Only used if the model is specified using
#' \code{ematrix}, rather than \code{hmodel}.
#' 
#' This must be a single formula - lists are not supported, unlike
#' \code{covariates}.  If a different model on each probability is required,
#' include all covariates in this formula, and use \code{fixedpars} to fix some
#' of their effects (for particular probabilities) at their default initial
#' values of zero.
#' @param misccovinits Initial values for the covariates on the
#' misclassification probabilities, defined in the same way as \code{covinits}.
#' Only used if the model is specified using \code{ematrix}.
#' @param miscconstraint A list of one vector for each named covariate on
#' misclassification probabilities. The vector indicates which covariate
#' effects on misclassification probabilities are constrained to be equal,
#' analogously to \code{constraint}.  Only used if the model is specified using
#' \code{ematrix}.
#' @param hcovariates List of formulae the same length as \code{hmodel},
#' defining any covariates governing the hidden Markov outcome models.  The
#' covariates operate on a suitably link-transformed linear scale, for example,
#' log scale for a Poisson outcome model. If there are no covariates for a
#' certain hidden state, then insert a NULL in the corresponding place in the
#' list.  For example, \code{hcovariates = list(~acute + age, ~acute, NULL).}
#' @param hcovinits Initial values for the hidden Markov model covariate
#' effects. A list of the same length as \code{hcovariates}. Each element is a
#' vector with initial values for the effect of each covariate on that state.
#' For example, the above \code{hcovariates} can be initialised with
#' \code{hcovariates = list(c(-8, 0), -8, NULL)}. Initial values must be given
#' for all or no covariates, if none are given these are all set to zero.  The
#' initial value given in the \code{hmodel} constructor function for the
#' corresponding baseline parameter is interpreted as the value of that
#' parameter with any covariates fixed to their means in the data.  If multiple
#' effects are constrained to be equal using \code{hconstraint}, then the
#' initial value is taken from the first of the multiple initial values
#' supplied.
#' @param hconstraint A named list. Each element is a vector of constraints on
#' the named hidden Markov model parameter. The vector has length equal to the
#' number of times that class of parameter appears in the whole model.
#' 
#' For example consider the three-state hidden Markov model described above,
#' with normally-distributed outcomes for states 1 and 2.  To constrain the
#' outcome variance to be equal for states 1 and 2, and to also constrain the
#' effect of \code{acute} on the outcome mean to be equal for states 1 and 2,
#' specify
#' 
#' \code{hconstraint = list(sd = c(1,1), acute=c(1,1))}
#' 
#' Note this excludes initial state occupancy probabilities and covariate
#' effects on those probabilities, which cannot be constrained.
#' @param hranges Range constraints for hidden Markov model parameters.
#' Supplied as a named list, with each element corresponding to the named
#' hidden Markov model parameter.  This element is itself a list with two
#' elements, vectors named "lower" and "upper".  These vectors each have length
#' equal to the number of times that class of parameter appears in the whole
#' model, and give the corresponding mininum amd maximum allowable values for
#' that parameter.  Maximum likelihood estimation is performed with these
#' parameters constrained in these ranges (through a log or logit-type
#' transformation).  Lower bounds of \code{-Inf} and upper bounds of \code{Inf}
#' can be given if the parameter is unbounded above or below.
#' 
#' For example, in the three-state model above, to constrain the mean for state
#' 1 to be between 0 and 6, and the mean of state 2 to be between 7 and 12,
#' supply
#' 
#' \code{hranges=list(mean=list(lower=c(0, 7), upper=c(6, 12)))}
#' 
#' These default to the natural ranges, e.g. the positive real line for
#' variance parameters, and [0,1] for probabilities.  Therefore \code{hranges}
#' need not be specified for such parameters unless an even stricter constraint
#' is desired.  If only one limit is supplied for a parameter, only the first
#' occurrence of that parameter is constrained.
#' 
#' Initial values should be strictly within any ranges, and not on the range
#' boundary, otherwise optimisation will fail with a "non-finite value" error.
#' @param qconstraint A vector of indicators specifying which baseline
#' transition intensities are equal. For example,
#' 
#' \code{qconstraint = c(1,2,3,3)}
#' 
#' constrains the third and fourth intensities to be equal, in a model with
#' four allowed instantaneous transitions.  When there are covariates on the
#' intensities and \code{center=TRUE} (the default), \code{qconstraint} is
#' applied to the intensities with covariates taking the values of the means in
#' the data.  When \code{center=FALSE}, \code{qconstraint} is applied to the
#' intensities with covariates set to zero.
#' @param econstraint A similar vector of indicators specifying which baseline
#' misclassification probabilities are constrained to be equal.  Only used if
#' the model is specified using \code{ematrix}, rather than \code{hmodel}.
#' @param initprobs Only used in hidden Markov models.  Underlying state
#' occupancy probabilities at each subject's first observation.  Can either be
#' a vector of \eqn{nstates} elements with common probabilities to all
#' subjects, or a \eqn{nsubjects} by \eqn{nstates} matrix of subject-specific
#' probabilities.  This refers to observations after missing data and subjects
#' with only one observation have been excluded.
#' 
#' If these are estimated (see \code{est.initprobs}), then this represents an
#' initial value, and defaults to equal probability for each state.  Otherwise
#' this defaults to \code{c(1, rep(0, nstates-1))}, that is, in state 1 with a
#' probability of 1.  Scaled to sum to 1 if necessary.  The state 1 occupancy
#' probability should be non-zero.
#' @param est.initprobs Only used in hidden Markov models.  If \code{TRUE},
#' then the underlying state occupancy probabilities at the first observation
#' will be estimated, starting from a vector of initial values supplied in the
#' \code{initprobs} argument.  Structural zeroes are allowed: if any of these
#' initial values are zero they will be fixed at zero during optimisation, even
#' if \code{est.initprobs=TRUE}, and no covariate effects on them are
#' estimated.  The exception is state 1, which should have non-zero occupancy
#' probability.
#' 
#' Note that the free parameters during this estimation exclude the state 1
#' occupancy probability, which is fixed at one minus the sum of the other
#' probabilities.
#' @param initcovariates Formula representing covariates on the initial state
#' occupancy probabilities, via multinomial logistic regression.  The linear
#' effects of these covariates, observed at the individual's first observation
#' time, operate on the log ratio of the state \eqn{r} occupancy probability to
#' the state 1 occupancy probability, for each \eqn{r = 2} to the number of
#' states.  Thus the state 1 occupancy probability should be non-zero. If
#' \code{est.initprobs} is \code{TRUE}, these effects are estimated starting
#' from their initial values.  If \code{est.initprobs} is \code{FALSE}, these
#' effects are fixed at theit initial values.
#' @param initcovinits Initial values for the covariate effects
#' \code{initcovariates}.  A named list with each element corresponding to a
#' covariate, as in \code{covinits}. Each element is a vector with (1 - number
#' of states) elements, containing the initial values for the linear effect of
#' that covariate on the log odds of that state relative to state 1, from state
#' 2 to the final state.  If \code{initcovinits} is not specified, all
#' covariate effects are initialised to zero.
#' @param deathexact Vector of indices of absorbing states whose time of entry
#' is known exactly, but the individual is assumed to be in an unknown
#' transient state ("alive") at the previous instant.  This is the usual
#' situation for times of death in chronic disease monitoring data.  For
#' example, if you specify \code{deathexact = c(4, 5)} then states 4 and 5 are
#' assumed to be exactly-observed death states.
#' 
#' See the \code{obstype} argument.  States of this kind correspond to
#' \code{obstype=3}. \code{deathexact = TRUE} indicates that the final
#' absorbing state is of this kind, and \code{deathexact = FALSE} or
#' \code{deathexact = NULL} (the default) indicates that there is no state of
#' this kind.
#' 
#' The \code{deathexact} argument is overridden by \code{obstype} or
#' \code{exacttimes}.
#' 
#' Note that you do not always supply a \code{deathexact} argument, even if
#' there are states that correspond to deaths, because they do not necessarily
#' have \code{obstype=3}.  If the state is known between the time of death and
#' the previous observation, then you should specify \code{obstype=2} for the
#' death times, or \code{exacttimes=TRUE} if the state is known at all times,
#' and the \code{deathexact} argument is ignored.
#' @param death Old name for the \code{deathexact} argument.  Overridden by
#' \code{deathexact} if both are supplied. Deprecated.
#' @param censor A state, or vector of states, which indicates censoring.
#' Censoring means that the observed state is known only to be one of a
#' particular set of states. For example, \code{censor=999} indicates that all
#' observations of \code{999} in the vector of observed states are censored
#' states.  By default, this means that the true state could have been any of
#' the transient (non-absorbing) states. To specify corresponding true states
#' explicitly, use a \code{censor.states} argument.
#' 
#' Note that in contrast to the usual terminology of survival analysis, here it
#' is the \emph{state} which is considered to be censored, rather than the
#' \emph{event time}.  If at the end of a study, an individual has not died,
#' but their true state is \emph{known}, then \code{censor} is unnecessary,
#' since the standard multi-state model likelihood is applicable.  Also a
#' "censored" state here can be at any time, not just at the end.
#' 
#' For hidden Markov models, censoring may indicate either a set of possible
#' observed states, or a set of (hidden) true states. The later case is
#' specified by setting the relevant elements of \code{obstrue} to 1 (and
#' \code{NA} otherwise).
#' 
#' Note in particular that general time-inhomogeneous Markov models with
#' piecewise constant transition intensities can be constructed using the
#' \code{censor} facility. If the true state is unknown on occasions when a
#' piecewise constant covariate is known to change, then censored states can be
#' inserted in the data on those occasions.  The covariate may represent time
#' itself, in which case the \code{pci} option to msm can be used to perform
#' this trick automatically, or some other time-dependent variable.
#' 
#' Not supported for multivariate hidden Markov models specified with
#' \code{\link{hmmMV}}.
#' @param censor.states Specifies the underlying states which censored
#' observations can represent. If \code{censor} is a single number (the
#' default) this can be a vector, or a list with one element.  If \code{censor}
#' is a vector with more than one element, this should be a list, with each
#' element a vector corresponding to the equivalent element of \code{censor}.
#' For example
#' 
#' \code{censor = c(99, 999), censor.states = list(c(2,3), c(3,4))}
#' 
#' means that observations coded 99 represent either state 2 or state 3, while
#' observations coded 999 are really either state 3 or state 4.
#' @param pci Model for piecewise-constant intensities.  Vector of cut points
#' defining the times, since the start of the process, at which intensities
#' change for all subjects.  For example
#' 
#' \code{pci = c(5, 10)}
#' 
#' specifies that the intensity changes at time points 5 and 10.  This will
#' automatically construct a model with a categorical (factor) covariate called
#' \code{timeperiod}, with levels \code{"[-Inf,5)"}, \code{"[5,10)"} and
#' \code{"[10,Inf)"}, where the first level is the baseline.  This covariate
#' defines the time period in which the observation was made.  Initial values
#' and constraints on covariate effects are specified the same way as for a
#' model with a covariate of this name, for example,
#' 
#' \code{covinits = list("timeperiod[5,10)"=c(0.1,0.1),
#' "timeperiod[10,Inf)"=c(0.1,0.1))}
#' 
#' Thus if \code{pci} is supplied, you cannot have a previously-existing
#' variable called \code{timeperiod} as a covariate in any part of a \code{msm}
#' model.
#' 
#' To assume piecewise constant intensities for some transitions but not others
#' with \code{pci}, use the \code{fixedpars} argument to fix the appropriate
#' covariate effects at their default initial values of zero.
#' 
#' Internally, this works by inserting censored observations in the data at
#' times when the intensity changes but the state is not observed.
#' 
#' If the supplied times are outside the range of the time variable in the
#' data, \code{pci} is ignored and a time-homogeneous model is fitted.
#' 
#' After fitting a time-inhomogeneous model, \code{\link{qmatrix.msm}} can be
#' used to obtain the fitted intensity matrices for each time period, for
#' example,
#' 
#' \code{qmatrix.msm(example.msm, covariates=list(timeperiod="[5,Inf)"))}
#' 
#' This facility does not support interactions between time and other
#' covariates.  Such models need to be specified "by hand", using a state
#' variable with censored observations inserted.  Note that the \code{data}
#' component of the \code{msm} object returned from a call to \code{msm} with
#' \code{pci} supplied contains the states with inserted censored observations
#' and time period indicators.  These can be used to construct such models.
#' 
#' Note that you do not need to use \code{pci} in order to model the effect of
#' a time-dependent covariate in the data.  \code{msm} will automatically
#' assume that covariates are piecewise-constant and change at the times when
#' they are observed.  \code{pci} is for when you want all intensities to
#' change at the same pre-specified times for all subjects.
#'
#' \code{pci} is not supported for multivariate hidden Markov models specified with
#' \code{\link{hmmMV}}.   An approximate equivalent can be constructed by
#' creating a variable in the data to represent the time period, and treating
#' that as a covariate using the \code{covariates} argument to \code{msm}.
#' This will assume that the value of this variable is constant between
#' observations. 
#' 
#' @param phase.states Indices of states which have a two-phase sojourn
#' distribution.  This defines a semi-Markov model, in which the hazard of an
#' onward transition depends on the time spent in the state.
#' 
#' This uses the technique described by Titman and Sharples (2009).  A hidden
#' Markov model is automatically constructed on an expanded state space, where
#' the phases correspond to the hidden states.  The "tau" proportionality
#' constraint described in this paper is currently not supported.
#' 
#' Covariates, constraints, \code{deathexact} and \code{censor} are expressed
#' with respect to the expanded state space.  If not supplied by hand,
#' \code{initprobs} is defined automatically so that subjects are assumed to
#' begin in the first of the two phases.
#' 
#' Hidden Markov models can additionally be given phased states.  The user
#' supplies an outcome distribution for each original state using
#' \code{hmodel}, which is expanded internally so that it is assumed to be the
#' same within each of the phased states.  \code{initprobs} is interpreted on
#' the expanded state space.  Misclassification models defined using
#' \code{ematrix} are not supported, and these must be defined using
#' \code{hmmCat} or \code{hmmIdent} constructors, as described in the
#' \code{hmodel} section of this help page.  Or the HMM on the expanded state
#' space can be defined by hand.
#' 
#' Output functions are presented as it were a hidden Markov model on the
#' expanded state space, for example, transition probabilities between states,
#' covariate effects on transition rates, or prevalence counts, are not
#' aggregated over the hidden phases.
#' 
#' Numerical estimation will be unstable when there is weak evidence for a
#' two-phase sojourn distribution, that is, if the model is close to Markov.
#' 
#' See \code{\link{d2phase}} for the definition of the two-phase distribution
#' and the interpretation of its parameters.
#' 
#' This is an experimental feature, and some functions are not implemented.
#' Please report any experiences of using this feature to the author!
#' @param phase.inits Initial values for phase-type models.  A list with one
#' component for each "two-phased" state.  Each component is itself a list of
#' two elements.  The first of these elements is a scalar defining the
#' transition intensity from phase 1 to phase 2.  The second element is a
#' matrix, with one row for each potential destination state from the
#' two-phased state, and two columns.  The first column is the transition rate
#' from phase 1 to the destination state, and the second column is the
#' transition rate from phase 2 to the destination state.  If there is only one
#' destination state, then this may be supplied as a vector.
#' 
#' In phase type models, the initial values for transition rates out of
#' non-phased states are taken from the \code{qmatrix} supplied to msm, and
#' entries of this matrix corresponding to transitions out of phased states are
#' ignored.
#' @param exacttimes By default, the transitions of the Markov process are
#' assumed to take place at unknown occasions in between the observation times.
#' If \code{exacttimes} is set to \code{TRUE}, then the observation times are
#' assumed to represent the exact times of transition of the process.  The
#' subject is assumed to be in the same state between these times.  An
#' observation may represent a transition to a different state or a repeated
#' observation of the same state (e.g. at the end of follow-up).  This is
#' equivalent to every row of the data having \code{obstype = 2}.  See the
#' \code{obstype} argument.  If both \code{obstype} and \code{exacttimes} are
#' specified then \code{exacttimes} is ignored.
#' 
#' Note that the complete history of the multi-state process is known with this
#' type of data.  The models which \pkg{msm} fits have the strong assumption of
#' constant (or piecewise-constant) transition rates.  Knowing the exact
#' transition times allows more realistic models to be fitted with other
#' packages.  For example parametric models with sojourn distributions more
#' flexible than the exponential can be fitted with the \pkg{flexsurv} package,
#' or semi-parametric models can be implemented with \pkg{survival} in
#' conjunction with \pkg{mstate}.
#' @param cl Width of symmetric confidence intervals for maximum likelihood
#' estimates, by default 0.95.
#' @param fixedpars Vector of indices of parameters whose values will be fixed
#' at their initial values during the optimisation. These are given in the
#' order: transition intensities (reading across rows of the transition
#' matrix), covariates on intensities (ordered by intensities within
#' covariates), hidden Markov model parameters, including misclassification
#' probabilities or parameters of HMM outcome distributions (ordered by
#' parameters within states), hidden Markov model covariate parameters (ordered
#' by covariates within parameters within states), initial state occupancy
#' probabilities (excluding the first probability, which is fixed at one minus
#' the sum of the others).
#' 
#' If there are equality constraints on certain parameters, then
#' \code{fixedpars} indexes the set of unique parameters, excluding those which
#' are constrained to be equal to previous parameters.
#' 
#' To fix all parameters, specify \code{fixedpars = TRUE}.
#' 
#' This can be useful for profiling likelihoods, and building complex models
#' stage by stage.
#' @param center If \code{TRUE} (the default, unless \code{fixedpars=TRUE})
#' then covariates are centered at their means during the maximum likelihood
#' estimation. This usually improves stability of the numerical optimisation.
#' @param opt.method If "optim", "nlm" or "bobyqa", then the corresponding R
#' function will be used for maximum likelihood estimation.
#' \code{\link{optim}} is the default.  "bobyqa" requires the package
#' \pkg{minqa} to be installed.  See the help of these functions for further
#' details.  Advanced users can also add their own optimisation methods, see
#' the source for \code{optim.R} in msm for some examples.
#' 
#' If "fisher", then a specialised Fisher scoring method is used (Kalbfleisch
#' and Lawless, 1985) which can be faster than the generic methods, though less
#' robust.  This is only available for Markov models with panel data
#' (\code{obstype=1}), that is, not for models with censored states, hidden
#' Markov models, exact observation or exact death times (\code{obstype=2,3}).
#' @param hessian If \code{TRUE} then standard errors and confidence intervals
#' are obtained from a numerical estimate of the Hessian (the observed
#' information matrix).  This is the default when maximum likelihood estimation
#' is performed.  If all parameters are fixed at their initial values and no
#' optimisation is performed, then this defaults to \code{FALSE}.  If
#' requested, the actual Hessian is returned in \code{x$paramdata$opt$hessian},
#' where \code{x} is the fitted model object.
#' 
#' If \code{hessian} is set to \code{FALSE}, then standard errors and
#' confidence intervals are obtained from the Fisher (expected) information
#' matrix, if this is available.  This may be preferable if the numerical
#' estimation of the Hessian is computationally intensive, or if the resulting
#' estimate is non-invertible or not positive definite.
#' @param use.deriv If \code{TRUE} then analytic first derivatives are used in
#' the optimisation of the likelihood, where available and an appropriate
#' quasi-Newton optimisation method, such as BFGS, is being used.  Analytic
#' derivatives are not available for all models.
#' @param use.expm If \code{TRUE} then any matrix exponentiation needed to
#' calculate the likelihood is done using the \pkg{expm} package.  Otherwise
#' the original routines used in \pkg{msm} 1.2.4 and earlier are used.  Set to
#' \code{FALSE} for backward compatibility, and let the package maintainer know
#' if this gives any substantive differences.
#' @param analyticp By default, the likelihood for certain simpler 3, 4 and 5
#' state models is calculated using an analytic expression for the transition
#' probability (P) matrix. For all other models, matrix exponentiation is used
#' to obtain P. To revert to the original method of using the matrix
#' exponential for all models, specify \code{analyticp=FALSE}. See the PDF
#' manual for a list of the models for which analytic P matrices are
#' implemented.
#' @param na.action What to do with missing data: either \code{na.omit} to drop
#' it and carry on, or \code{na.fail} to stop with an error.  Missing data
#' includes all NAs in the states, times, \code{subject} or \code{obstrue}, all
#' NAs at the first observation for a subject for covariates in
#' \code{initcovariates}, all NAs in other covariates (excluding the last
#' observation for a subject), all NAs in \code{obstype} (excluding the first
#' observation for a subject), and any subjects with only one observation (thus
#' no observed transitions).
#' @param ... Optional arguments to the general-purpose optimisation routine,
#' \code{\link{optim}} by default. For example \code{method="Nelder-Mead"} to
#' change the optimisation algorithm from the \code{"BFGS"} method that msm
#' calls by default.
#' 
#' It is often worthwhile to normalize the optimisation using
#' \code{control=list(fnscale = a)}, where \code{a} is the a number of the
#' order of magnitude of the -2 log likelihood.
#' 
#' If 'false' convergence is reported and the standard errors cannot be
#' calculated due to a non-positive-definite Hessian, then consider tightening
#' the tolerance criteria for convergence. If the optimisation takes a long
#' time, intermediate steps can be printed using the \code{trace} argument of
#' the control list. See \code{\link{optim}} for details.
#' 
#' For the Fisher scoring method, a \code{control} list can be supplied in the
#' same way, but the only supported options are \code{reltol}, \code{trace} and
#' \code{damp}.  The first two are used in the same way as for
#' \code{\link{optim}}.  If the algorithm fails with a singular information
#' matrix, adjust \code{damp} from the default of zero (to, e.g. 1).  This adds
#' a constant identity matrix multiplied by \code{damp} to the information
#' matrix during optimisation.
#' @return To obtain summary information from models fitted by the
#' \code{\link{msm}} function, it is recommended to use extractor functions
#' such as \code{\link{qmatrix.msm}}, \code{\link{pmatrix.msm}},
#' \code{\link{sojourn.msm}}, \code{\link{msm.form.qoutput}}.  These provide
#' estimates and confidence intervals for quantities such as transition
#' probabilities for given covariate values.
#' 
#' For advanced use, it may be necessary to directly use information stored in
#' the object returned by \code{\link{msm}}.  This is documented in the help
#' page \code{\link{msm.object}}.
#' 
#' Printing a \code{msm} object by typing the object's name at the command line
#' implicitly invokes \code{\link{print.msm}}.  This formats and prints the
#' important information in the model fit, and also returns that information in
#' an R object.  This includes estimates and confidence intervals for the
#' transition intensities and (log) hazard ratios for the corresponding
#' covariates. When there is a hidden Markov model, the chief information in
#' the \code{hmodel} component is also formatted and printed. This includes
#' estimates and confidence intervals for each parameter.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{simmulti.msm}}, \code{\link{plot.msm}},
#' \code{\link{summary.msm}}, \code{\link{qmatrix.msm}},
#' \code{\link{pmatrix.msm}}, \code{\link{sojourn.msm}}.
#' @references Jackson, C.H. (2011). Multi-State Models for Panel Data: The msm
#' Package for R., Journal of Statistical Software, 38(8), 1-29. URL
#' http://www.jstatsoft.org/v38/i08/.
#' 
#' Kalbfleisch, J., Lawless, J.F., The analysis of panel data under a Markov
#' assumption \emph{Journal of the Americal Statistical Association} (1985)
#' 80(392): 863--871.
#' 
#' Kay, R.  A Markov model for analysing cancer markers and disease states in
#' survival studies.  \emph{Biometrics} (1986) 42: 855--865.
#' 
#' Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P.  Multi-state
#' Markov models for analysing incomplete disease history data with
#' illustrations for HIV disease.  \emph{Statistics in Medicine} (1994) 13(3):
#' 805--821.
#' 
#' Satten, G.A. and Longini, I.M.  Markov chains with measurement error:
#' estimating the 'true' course of a marker of the progression of human
#' immunodeficiency virus disease (with discussion) \emph{Applied Statistics}
#' 45(3): 275-309 (1996)
#' 
#' Jackson, C.H. and Sharples, L.D. Hidden Markov models for the onset and
#' progression of bronchiolitis obliterans syndrome in lung transplant
#' recipients \emph{Statistics in Medicine}, 21(1): 113--128 (2002).
#' 
#' Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W. and Couto, E.
#' Multi-state Markov models for disease progression with classification error.
#' \emph{The Statistician}, 52(2): 193--209 (2003)
#' 
#' Titman, A.C. and Sharples, L.D.  Semi-Markov models with phase-type sojourn
#' distributions.  \emph{Biometrics} 66, 742-752 (2009).
#' 
#' Venables, W.N. and Ripley, B.D. (2002) \emph{Modern Applied Statistics with
#' S}, second edition. Springer.
#' @keywords models
#' @examples
#' 
#' ### Heart transplant data
#' ### For further details and background to this example, see
#' ### Jackson (2011) or the PDF manual in the doc directory.
#' print(cav[1:10,])
#' twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
#' c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
#' statetable.msm(state, PTNUM, data=cav)
#' crudeinits.msm(state ~ years, PTNUM, data=cav, qmatrix=twoway4.q)
#' cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
#'                  qmatrix = twoway4.q, deathexact = 4, 
#'                  control = list ( trace = 2, REPORT = 1 )  )
#' cav.msm
#' qmatrix.msm(cav.msm)
#' pmatrix.msm(cav.msm, t=10)
#' sojourn.msm(cav.msm)
#' 
#' @export msm
msm <- function(formula, subject=NULL, data=list(), qmatrix, gen.inits=FALSE,
                ematrix=NULL,  hmodel=NULL, obstype=NULL, obstrue=NULL,
                covariates = NULL, covinits = NULL, constraint = NULL,
                misccovariates = NULL, misccovinits = NULL, miscconstraint = NULL,
                hcovariates = NULL, hcovinits = NULL, hconstraint = NULL, hranges=NULL,
                qconstraint=NULL, econstraint=NULL, initprobs = NULL,
                est.initprobs=FALSE, initcovariates = NULL, initcovinits = NULL,
                deathexact = NULL, death = NULL, exacttimes = FALSE, censor=NULL,
                censor.states=NULL, pci=NULL, phase.states=NULL,
                phase.inits = NULL, # TODO merge with inits eventually
                cl = 0.95, fixedpars = NULL, center=TRUE,
                opt.method="optim", hessian=NULL, use.deriv=TRUE,
                use.expm=TRUE, analyticp=TRUE, na.action=na.omit, ...)
{
    call <- match.call()
    if (missing(formula)) stop("state ~ time formula not given")
    if (missing(data)) data <- environment(formula)

### MODEL FOR TRANSITION INTENSITIES
    if (gen.inits) {
        if (is.null(hmodel) && is.null(ematrix)) {
            subj <- eval(substitute(subject), data, parent.frame())
            qmatrix <- crudeinits.msm(formula, subj, qmatrix, data, censor, censor.states)
        }
        else warning("gen.inits not supported for hidden Markov models, ignoring")
    }
    qmodel <- qmodel.orig <- msm.form.qmodel(qmatrix, qconstraint, analyticp, use.expm, phase.states)

    if (!is.null(phase.states)) {
        qmodel <- msm.phase2qmodel(qmodel, phase.states, phase.inits, qconstraint, analyticp, use.expm)
    } 
### MISCLASSIFICATION MODEL
    if (!is.null(ematrix)) {
        msm.check.ematrix(ematrix, qmodel.orig$nstates)
        if (!is.null(phase.states)){
            stop("phase-type models with additional misclassification must be specified through \"hmodel\" with hmmCat() or hmmIdent() constructors, or as HMMs by hand") 
        }
        emodel <- msm.form.emodel(ematrix, econstraint, initprobs, est.initprobs, qmodel)
    }
    else emodel <- list(misc=FALSE, npars=0, ndpars=0)

### GENERAL HIDDEN MARKOV MODEL
    if (!is.null(hmodel)) {
        msm.check.hmodel(hmodel, qmodel.orig$nstates)
        if (!is.null(phase.states)){
            hmodel.orig <- hmodel
            hmodel <- rep(hmodel, qmodel$phase.reps)
        }
        hmodel <- msm.form.hmodel(hmodel, hconstraint, initprobs, est.initprobs)
    }
    else {
        if (!is.null(hcovariates)) stop("hcovariates have been specified, but no hmodel")
        if (!is.null(phase.states)){
            hmodel <- msm.phase2hmodel(qmodel, hmodel)
        }
        else hmodel <- list(hidden=FALSE, models=rep(0, qmodel$nstates), nipars=0, nicoveffs=0, totpars=0, ncoveffs=0) # might change later if misc
    }
### CONVERT OLD STYLE MISCLASSIFICATION MODEL TO NEW GENERAL HIDDEN MARKOV MODEL
    if (emodel$misc) {
        hmodel <- msm.emodel2hmodel(emodel, qmodel)
    }
    else {
        emodel <- list(misc=FALSE, npars=0, ndpars=0, nipars=0, nicoveffs=0)
        hmodel$ematrix <- FALSE
    }
    
### EXACT DEATH TIMES. Logical values allowed for backwards compatibility (TRUE means final state has exact death time, FALSE means no states with exact death times)
    if (!is.null(deathexact)) death <- deathexact
    dmodel <- msm.form.dmodel(death, qmodel, hmodel)  # returns death, ndeath,
    if (dmodel$ndeath > 0 && exacttimes) warning("Ignoring death argument, as all states have exact entry times")

### CENSORING MODEL
    cmodel <- msm.form.cmodel(censor, censor.states, qmodel$qmatrix, hmodel)

### SOME CHECKS    
    if (!inherits(formula, "formula")) stop("formula is not a formula")
    if (!is.null(covariates) && (!(is.list(covariates) || inherits(covariates, "formula"))))
        stop(deparse(substitute(covariates)), " should be a formula or list of formulae")
    if (!is.null(misccovariates) && (!inherits(misccovariates, "formula")))
        stop(deparse(substitute(misccovariates)), " should be a formula")
    if (is.list(covariates)) {  # different covariates on each intensity
        covlist <- covariates
        msm.check.covlist(covlist, qmodel)
        ter <- lapply(covlist, function(x)attr(terms(x),"term.labels"))
        covariates <- reformulate(unique(unlist(ter))) # merge the formulae into one
    } else covlist <- NULL # parameters will be constrained later, see msm.form.cri
    if (is.null(covariates)) covariates <- ~1
    if (emodel$misc && is.null(misccovariates)) misccovariates <- ~1
    if (hmodel$hidden && !is.null(hcovariates)) msm.check.hcovariates(hcovariates, qmodel)
    
### BUILD MODEL FRAME containing all data required for model fit,
### using all variables found in formulae.
    ## Names include factor() around covariate names, and interactions,
    ## if specified.  Need to build and evaluate a call, instead of
    ## running model.frame() directly, to find subject and other
    ## extras. Not sure why.
    indx <- match(c("data", "subject", "obstrue"), names(call), nomatch = 0)
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    temp[["state"]] <- as.name(all.vars(formula[[2]]))
    temp[["time"]] <- as.name(all.vars(formula[[3]]))
    varnames <- function(x){if(is.null(x)) NULL else attr(terms(x), "term.labels")}
    forms <- c(covariates, misccovariates, hcovariates, initcovariates)
    covnames <- unique(unlist(lapply(forms, varnames)))
    temp[["formula"]] <- if (length(covnames) > 0) reformulate(covnames) else ~1
    temp[["na.action"]] <- na.pass # run na.action later so we can pass aux info to it
    temp[["data"]] <- data
    mf <- eval(temp, parent.frame())
    
    ## remember user-specified names for later (e.g. bootstrap/cross validation)
    usernames <- c(state=all.vars(formula[[2]]), time=all.vars(formula[[3]]), subject=as.character(temp$subject), obstype=as.character(substitute(obstype)), obstrue=as.character(temp$obstrue))
    attr(mf, "usernames") <- usernames

    ## handle matrices in state outcome constructed in formula with cbind()
    indx <- match(c("formula", "data"), names(call), nomatch = 0)
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    mfst <- eval(temp, parent.frame())
    if (is.matrix(mfst[[1]]) && !is.matrix(mf$"(state)"))
        mf$"(state)" <- mfst[[1]]

    if (is.factor(mf$"(state)")){
        if (!all(grepl("^[[:digit:]]+$", as.character(mf$"(state)"))))
            stop("state variable should be numeric or a factor with ordinal numbers as levels")
        else mf$"(state)" <- as.numeric(as.character(mf$"(state)"))
    } else if (is.character(mf$"(state)")) stop("state variable is character, should be numeric") 
    msm.check.state(qmodel$nstates, mf$"(state)", cmodel$censor, hmodel)
    if (is.null(mf$"(subject)")) mf$"(subject)" <- rep(1, nrow(mf))
    msm.check.times(mf$"(time)", mf$"(subject)", mf$"(state)")
    obstype <- if (missing(obstype)) NULL else eval(substitute(obstype), data, parent.frame()) # handle separately to allow passing a scalar (1, 2 or 3)
    mf$"(obstype)" <- msm.form.obstype(mf, obstype, dmodel, exacttimes)
    mf$"(obstrue)" <- msm.form.obstrue(mf, hmodel, cmodel)
    mf$"(obs)" <- seq_len(nrow(mf)) # row numbers before NAs omitted, for reporting in msm.check.*
    basenames <- c("(state)","(time)","(subject)","(obstype)","(obstrue)","(obs)")
    attr(mf, "covnames") <- setdiff(names(mf), basenames)
    attr(mf, "covnames.q") <- colnames(attr(terms(covariates), "factors")) # names as in data, plus factor() if user
    if (emodel$misc) attr(mf, "covnames.e") <- colnames(attr(terms(misccovariates), "factors"))
    attr(mf, "ncovs") <- length(attr(mf, "covnames"))


### HANDLE MISSING DATA
    ## Find which cols of model frame correspond to covs which appear only on initprobs
    ## Pass through to na.action as attribute
    ic <- all.vars(initcovariates)
    others <- c(covariates,misccovariates,hcovariates)
    oic <- ic[!ic %in% unlist(lapply(others, all.vars))]
    attr(mf, "icovi") <- match(oic, colnames(mf))
    if (missing(na.action) || identical(na.action, na.omit) || (identical(na.action,"na.omit")))
        mf <- na.omit.msmdata(mf, hidden=hmodel$hidden, misc=emodel$misc)
    else if (identical(na.action, na.fail) || (identical(na.action,"na.fail")))
        mf <- na.fail.msmdata(mf, hidden=hmodel$hidden, misc=emodel$misc)
    else stop ("na.action should be \"na.omit\" or \"na.fail\"")
    attr(mf, "npts") <- length(unique(mf$"(subject)"))
    attr(mf, "ntrans") <- nrow(mf) - attr(mf, "npts")

### UTILITY FOR PIECEWISE-CONSTANT INTENSITIES.  Insert censored
    ## observations into the model frame, add a factor "timeperiod" to
    ## the data, and add the corresponding term to the formula.
    ## NOTE pci kept pre-imputed data as msmdata.obs.orig
    if (!is.null(pci)) {
        if (isTRUE(hmodel$mv)) stop("`pci` not supported for multivariate hidden Markov models.  As an approximation, create a variable in your data representing the time period, and treat it as a covariate")
        tdmodel <- msm.pci(pci, mf, qmodel, cmodel, covariates)
        if (!is.null(tdmodel)) {
            mf <- tdmodel$mf; covariates <- tdmodel$covariates; cmodel <- tdmodel$cmodel
            pci <- tdmodel$tcut # was returned in msm object
        } else {pci <- NULL; mf$"(pci.imp)" <- 0}
    } else tdmodel <- NULL

### CALCULATE COVARIATE MEANS from data by transition, including means of 0/1 contrasts.
    forms <- c(covariates, misccovariates, hcovariates, initcovariates) # covariates may have been updated to include timeperiod for pci
    covnames <- unique(unlist(lapply(forms, varnames)))
    if (length(covnames) > 0) {
        mm.mean <- model.matrix(reformulate(covnames), mf)
        cm <- colMeans(mm.mean[duplicated(mf$"(subject)",fromLast=TRUE),,drop=FALSE])
        cm["(Intercept)"] <- 0
    } else cm <- NULL
    attr(mf, "covmeans") <- cm

### MAKE AGGREGATE DATA for computing likelihood of non-hidden models efficiently
### This can also be used externally to extract aggregate data from model frame in msm object
    mf.agg <- msm.form.mf.agg(list(data=list(mf=mf), qmodel=qmodel, hmodel=hmodel, cmodel=cmodel))
### Make indicator for which distinct P matrix each observation corresponds to. Only used in HMMs.
    mf$"(pcomb)" <- msm.form.hmm.agg(mf);

### CONVERT misclassification covariate formula to hidden covariate formula
    if (inherits(misccovariates, "formula")){
        if (!emodel$misc) stop("misccovariates supplied but no ematrix")
        hcovariates <- lapply(ifelse(rowSums(emodel$imatrix)>0, deparse(misccovariates), deparse(~1)), as.formula)
    }

### FORM DESIGN MATRICES FOR COVARIATE MODELS
### These can also be used externally to extract design matrices from model frame in msm object
    mm.cov <- msm.form.mm.cov(list(data=list(mf=mf), covariates=covariates, center=center))
    mm.cov.agg <- msm.form.mm.cov.agg(list(data=list(mf.agg=mf.agg), covariates=covariates, hmodel=hmodel, cmodel=cmodel, center=center))
    mm.mcov <- msm.form.mm.mcov(list(data=list(mf=mf), misccovariates=misccovariates, emodel=emodel, center=center))
    mm.hcov <- msm.form.mm.hcov(list(data=list(mf=mf), hcovariates=hcovariates, qmodel=qmodel, hmodel=hmodel, center=center))
    mm.icov <- msm.form.mm.icov(list(data=list(mf=mf), initcovariates=initcovariates, hmodel=hmodel, center=center))

### MODEL FOR COVARIATES ON INTENSITIES
    if (!is.null(covlist)) {
        cri <- msm.form.cri(covlist, qmodel, mf, mm.cov, tdmodel)
    } else cri <- NULL
    qcmodel <-
        if (ncol(mm.cov) > 1)
            msm.form.covmodel(mf, mm.cov, constraint, covinits, cm, qmodel$npars, cri)
        else {
            if (!is.null(constraint)) warning("constraint specified but no covariates")
            list(npars=0, ncovs=0, ndpars=0)
        }
### MODEL FOR COVARIATES ON MISCLASSIFICATION PROBABILITIES
    if (!emodel$misc || is.null(misccovariates))
        ecmodel <- list(npars=0, ncovs=0)
    if (!is.null(misccovariates)) {
        if (!emodel$misc) {
            warning("misccovariates have been specified, but misc is FALSE. Ignoring misccovariates.")
        }
        else {
            ecmodel <- msm.form.covmodel(mf, mm.mcov, miscconstraint, misccovinits, cm, emodel$npars, cri=NULL)
            hcovariates <- msm.misccov2hcov(misccovariates, emodel)
            hcovinits <- msm.misccovinits2hcovinits(misccovinits, hcovariates, emodel, ecmodel)
        }
    }
### MODEL FOR COVARIATES ON GENERAL HIDDEN PARAMETERS
    if (!is.null(hcovariates)) {
        if (hmodel$mv) stop("hcovariates not supported for multivariate hidden Markov models")
        hmodel <- msm.form.hcmodel(hmodel, mm.hcov, hcovinits, hconstraint)
        if (emodel$misc)
            hmodel$covconstr <- msm.form.hcovconstraint(miscconstraint, hmodel)
    }
    else if (hmodel$hidden) {
        npars <- if(hmodel$mv) colSums(hmodel$npars) else hmodel$npars
        hmodel <- c(hmodel, list(ncovs=rep(rep(0, hmodel$nstates), npars), ncoveffs=0))
        class(hmodel) <- "hmodel"
    }
    if (!is.null(initcovariates)) {
        if (hmodel$hidden)
            hmodel <- msm.form.icmodel(hmodel, mm.icov, initcovinits)
        else warning("initprobs and initcovariates ignored for non-hidden Markov models")
    }
    else if (hmodel$hidden) {
        hmodel <- c(hmodel, list(nicovs=rep(0, hmodel$nstates-1), nicoveffs=0, cri=ecmodel$cri))
        class(hmodel) <- "hmodel"
    }
    if (hmodel$hidden && !emodel$misc) {
        hmodel$constr <- msm.form.hconstraint(hconstraint, hmodel)
        hmodel$covconstr <- msm.form.hcovconstraint(hconstraint, hmodel)
    }
    if (hmodel$hidden) hmodel$ranges <- msm.form.hranges(hranges, hmodel)
### INITIAL STATE OCCUPANCY PROBABILITIES IN HMMS
    if (hmodel$hidden) hmodel <- msm.form.initprobs(hmodel, initprobs, mf)
### FORM LIST OF INITIAL PARAMETERS, MATCHING PROVIDED INITS WITH SPECIFIED MODEL, FIXING SOME PARS IF REQUIRED
    p <- msm.form.params(qmodel, qcmodel, emodel, hmodel, fixedpars)

    msmdata <- list(mf=mf, mf.agg=mf.agg, mm.cov=mm.cov, mm.cov.agg=mm.cov.agg,
                    mm.mcov=mm.mcov, mm.hcov=mm.hcov, mm.icov=mm.icov)

### CALCULATE LIKELIHOOD AT INITIAL VALUES OR DO OPTIMISATION (see optim.R)
    if (p$fixed) opt.method <- "fixed"
    if (is.null(hessian)) hessian <- !p$fixed
    p <- msm.optim(opt.method, p, hessian, use.deriv, msmdata, qmodel, qcmodel, cmodel, hmodel, ...)
    if (p$fixed) {
        p$foundse <- FALSE
        p$covmat <- NULL
    } else {
        p$params <- msm.rep.constraints(p$params, p, hmodel)
        hess <- if (hessian) p$opt$hessian else p$information
        if (!is.null(hess) &&
            all(!is.na(hess)) && all(!is.nan(hess)) && all(is.finite(hess)) && all(eigen(hess)$values > 0))
        {
            p$foundse <- TRUE
            p$covmat <- matrix(0, nrow=p$npars, ncol=p$npars)
            p$covmat[p$optpars,p$optpars] <- solve(0.5 * hess)
            p$covmat <- p$covmat[!duplicated(abs(p$constr)),!duplicated(abs(p$constr)), drop=FALSE][abs(p$constr),abs(p$constr), drop=FALSE]
            p$ci <- cbind(p$params - qnorm(1 - 0.5*(1-cl))*sqrt(diag(p$covmat)),
                          p$params + qnorm(1 - 0.5*(1-cl))*sqrt(diag(p$covmat)))
            p$ci[p$fixedpars,] <- NA
            for (i in 1:2)
                p$ci[,i] <- gexpit(p$ci[,i], p$ranges[,"lower",drop=FALSE], p$ranges[,"upper",drop=FALSE])
        }
        else {
            p$foundse <- FALSE
            p$covmat <- p$ci <- NULL
            if (!is.null(hess))
                warning("Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite.")
        }
    }
    p$estimates.t <- p$params  # Calculate estimates and CIs on natural scale
    p$estimates.t <- msm.inv.transform(p$params, hmodel, p$ranges)
    ## calculate CIs for misclassification probabilities (needs multivariate transform and delta method)
    if (any(p$plabs=="p") && p$foundse){
        p.se <- p.se.msm(x=list(data=msmdata,qmodel=qmodel,emodel=emodel,hmodel=hmodel,
                         qcmodel=qcmodel,ecmodel=ecmodel,paramdata=p,center=center),
                         covariates = if(center) "mean" else 0)
        p$ci[p$plabs %in% c("p","pbase"),] <- as.numeric(unlist(p.se[,c("LCL","UCL")]))
    }
    ## calculate CIs for initial state probabilities in HMMs (using normal simulation method)
    if (p$foundse && any(p$plabs=="initp"))  p <- initp.ci.msm(p, cl)

### FORM A MSM OBJECT FROM THE RESULTS
    msmobject <- list (
                       call = match.call(),
                       minus2loglik = p$lik,
                       deriv = p$deriv,
                       estimates = p$params,
                       estimates.t = p$estimates.t,
                       fixedpars = p$fixedpars,
                       center = center,
                       covmat = p$covmat,
                       ci = p$ci,
                       opt = p$opt,
                       foundse = p$foundse,
                       data = msmdata,
                       qmodel = qmodel,
                       emodel = emodel,
                       qcmodel = qcmodel,
                       ecmodel = ecmodel,
                       hmodel = hmodel,
                       cmodel = cmodel,
                       pci = pci,
                       paramdata=p,
                       cl=cl,
        covariates=covariates,
        misccovariates=misccovariates,
        hcovariates=hcovariates,
        initcovariates=initcovariates
                       )
    attr(msmobject, "fixed") <- p$fixed
    class(msmobject) <- "msm"

### Form lists of matrices from parameter estimates 
    msmobject <- msm.form.output(msmobject, "intens")
   
### Include intensity and misclassification matrices on natural scales
    q <- qmatrix.msm(msmobject, covariates=(if(center) "mean" else 0))
    msmobject$Qmatrices$baseline <- q$estimates
    msmobject$QmatricesSE$baseline <- q$SE
    msmobject$QmatricesL$baseline <- q$L
    msmobject$QmatricesU$baseline <- q$U

    if (hmodel$hidden) {
        msmobject$hmodel <- msm.form.houtput(hmodel, p, msmdata, cmodel)
    }
    if (emodel$misc) {
        msmobject <- msm.form.output(msmobject, "misc")
        e <- ematrix.msm(msmobject, covariates=(if(center) "mean" else 0))
        msmobject$Ematrices$baseline <- e$estimates
        msmobject$EmatricesSE$baseline <- e$SE
        msmobject$EmatricesL$baseline <- e$L
        msmobject$EmatricesU$baseline <- e$U
    }
    msmobject$msmdata[!names(msmobject$msmdata)=="mf"] <- NULL # only keep model frame in returned data.  drop at last minute, as might be needed in msm.form.houtput.
### Include mean sojourn times
    msmobject$sojourn <- sojourn.msm(msmobject, covariates=(if(center) "mean" else 0))
    msmobject
}

msm.check.qmatrix <- function(qmatrix)
{
    if (!is.numeric(qmatrix) || ! is.matrix(qmatrix))
        stop("qmatrix should be a numeric matrix")
    if (nrow(qmatrix) != ncol(qmatrix))
        stop("Number of rows and columns of qmatrix should be equal")
    q2 <- qmatrix; diag(q2) <- 0
    if (any(q2 < 0))
        stop("off-diagonal entries of qmatrix should not be negative")
    invisible()
}

msm.fixdiag.qmatrix <- function(qmatrix)
{
    diag(qmatrix) <- 0
    diag(qmatrix) <- - rowSums(qmatrix)
    qmatrix
}

msm.fixdiag.ematrix <- function(ematrix)
{
    diag(ematrix) <- 0
    diag(ematrix) <- 1 - rowSums(ematrix)
    ematrix
}

msm.form.qmodel <- function(qmatrix, qconstraint=NULL, analyticp=TRUE, use.expm=FALSE, phase.states=NULL)
{
    msm.check.qmatrix(qmatrix)
    nstates <- dim(qmatrix)[1]
    qmatrix <- msm.fixdiag.qmatrix(qmatrix)
    if (is.null(rownames(qmatrix)))
        rownames(qmatrix) <- colnames(qmatrix) <- paste("State", seq(nstates))
    else if (is.null(colnames(qmatrix))) colnames(qmatrix) <- rownames(qmatrix)
    imatrix <- ifelse(qmatrix > 0, 1, 0)
    inits <- t(qmatrix)[t(imatrix)==1]
    npars <- sum(imatrix)
    if (npars==0) stop("qmatrix contains all zeroes off the diagonal, so no transitions are permitted in the model")
    ## for phase-type models, leave processing qconstraint until after phased Q matrix has been formed in msm.phase2qmodel
    if (!is.null(qconstraint) && is.null(phase.states)) { 
        if (!is.numeric(qconstraint)) stop("qconstraint should be numeric")
        if (length(qconstraint) != npars)
            stop("baseline intensity constraint of length " ,length(qconstraint), ", should be ", npars)
        constr <- match(qconstraint, unique(qconstraint))
    }
    else
        constr <- 1:npars
    ndpars <- max(constr)
    ipars <- t(imatrix)[t(lower.tri(imatrix) | upper.tri(imatrix))]
    graphid <- paste(which(ipars==1), collapse="-")
    if (analyticp && graphid %in% names(.msm.graphs[[paste(nstates)]])) {
        ## analytic P matrix is implemented for this particular intensity matrix
        iso <- .msm.graphs[[paste(nstates)]][[graphid]]$iso
        perm <- .msm.graphs[[paste(nstates)]][[graphid]]$perm
        qperm <- order(perm) # diff def in 1.2.3, indexes q matrices not vectors
    }
    else {
        iso <- 0
        perm <- qperm <- NA
    }
    qmodel <- list(nstates=nstates, iso=iso, perm=perm, qperm=qperm,
                   npars=npars, imatrix=imatrix, qmatrix=qmatrix, inits=inits,
                   constr=constr, ndpars=ndpars, expm=as.numeric(use.expm))
    class(qmodel) <- "msmqmodel"
    qmodel
}

msm.check.ematrix <- function(ematrix, nstates)
{
    if (!is.numeric(ematrix) || ! is.matrix(ematrix))
        stop("ematrix should be a numeric matrix")
    if (nrow(ematrix) != ncol(ematrix))
        stop("Number of rows and columns of ematrix should be equal")
    if (!all(dim(ematrix) == nstates))
        stop("Dimensions of qmatrix and ematrix should be the same")
    if (!all ( ematrix >= 0 | ematrix <= 1) )
        stop("Not all elements of ematrix are between 0 and 1")
    invisible()
}

msm.form.emodel <- function(ematrix, econstraint=NULL, initprobs=NULL, est.initprobs, qmodel)
{
    diag(ematrix) <- 0
    imatrix <- ifelse(ematrix > 0 & ematrix < 1, 1, 0) # don't count as parameters if perfect misclassification (1.4.2 bug fix)
    diag(ematrix) <- 1 - rowSums(ematrix)
    if (is.null(rownames(ematrix)))
        rownames(ematrix) <- colnames(ematrix) <- paste("State", seq(qmodel$nstates))
    else if (is.null(colnames(ematrix))) colnames(ematrix) <- rownames(ematrix)
    dimnames(imatrix) <- dimnames(ematrix)
    npars <- sum(imatrix)
    nstates <- nrow(ematrix)
    inits <- t(ematrix)[t(imatrix)==1]
    if (is.null(initprobs)) {
        initprobs <- if (est.initprobs) rep(1/qmodel$nstates, qmodel$nstates) else c(1, rep(0, qmodel$nstates-1))
    }
    else {
        if (!is.numeric(initprobs)) stop("initprobs should be numeric")
        if (is.matrix(initprobs)) {
            if (ncol(initprobs) != qmodel$nstates) stop("initprobs matrix has ", ncol(initprobs), " columns, should be number of states = ", qmodel$nstates)
            if (est.initprobs) { warning("Not estimating initial state occupancy probabilities since supplied as a matrix") }
            initprobs <- initprobs / rowSums(initprobs)
            est.initprobs <- FALSE
        }
        else {
            if (length(initprobs) != qmodel$nstates) stop("initprobs vector of length ", length(initprobs), ", should be vector of length ", qmodel$nstates, " or a matrix")
            initprobs <- initprobs / sum(initprobs)
            if (est.initprobs && any(initprobs==1)) {
                est.initprobs <- FALSE
                warning("Not estimating initial state occupancy probabilities, since some are fixed to 1")
            }
        }
    }
    nipars <- if (est.initprobs) qmodel$nstates else 0
    if (!is.null(econstraint)) {
        if (!is.numeric(econstraint)) stop("econstraint should be numeric")
        if (length(econstraint) != npars)
            stop("baseline misclassification constraint of length " ,length(econstraint), ", should be ", npars)
        constr <- match(econstraint, unique(econstraint))
    }
    else
        constr <- seq(length.out=npars)
    ndpars <- if(npars>0) max(constr) else 0

    nomisc <- isTRUE(all.equal(as.vector(diag(ematrix)),rep(1,nrow(ematrix)))) # degenerate ematrix with no misclassification: all 1 on diagonal
    if (nomisc)
        emodel <- list(misc=FALSE, npars=0, ndpars=0)
    else                                     
        emodel <- list(misc=TRUE, npars=npars, nstates=nstates, imatrix=imatrix, ematrix=ematrix, inits=inits,
                       constr=constr, ndpars=ndpars, nipars=nipars, initprobs=initprobs, est.initprobs=est.initprobs)
    class(emodel) <- "msmemodel"
    emodel
}

### Check elements of state vector. For simple models and misc models specified with ematrix
### Only check performed for hidden models is that data and model dimensions match

msm.check.state <- function(nstates, state, censor, hmodel)
{
    if (hmodel$hidden){
        if (!is.null(ncol(state))) {
            pl1 <- if (ncol(state) > 1) "s" else ""
            pl2 <- if (max(hmodel$nout) > 1) "s" else ""
            if ((ncol(state) != max(hmodel$nout)) && (max(hmodel$nout) > 1))
                stop(sprintf("outcome matrix in data has %d column%s, but outcome models have a maximum of %d dimension%s", ncol(state), pl1, max(hmodel$nout), pl2))
        } else {
            if (max(hmodel$nout) > 1) {
                stop("Only one column in state outcome data, but multivariate hidden model supplied")
            }
        }

    }  else { 
        states <- c(1:nstates, censor)
        state <- na.omit(state) # NOTE added in 1.4
        if (!is.null(ncol(state)) && ncol(state) > 1) stop("Matrix outcomes only allowed for hidden Markov models")
        if (!is.null(state)) {
            statelist <- if (nstates==2) "1, 2" else if (nstates==3) "1, 2, 3" else paste("1, 2, ... ,",nstates)
            if (length(setdiff(unique(state), states)) > 0)
                stop("State vector contains elements not in ",statelist)
            miss.state <- setdiff(states, unique(state))
            ## Don't do this for misclassification models: it's OK to have particular state not observed by chance
            if (length(miss.state) > 0 && (!hmodel$hidden || !hmodel$ematrix))
                warning("State vector doesn't contain observations of ",paste(miss.state, collapse=","))
        }
    }
    invisible()
}

msm.check.times <- function(time, subject, state=NULL)
{
    final.rows <- !is.na(subject) & !is.na(time)
    if (!is.null(state)) {
        nas <- if (is.matrix(state)) apply(state, 1, function(x)all(is.na(x))) else is.na(state)
        final.rows <- final.rows & !nas
        state <- if (is.matrix(state)) state[final.rows, ,drop=FALSE] else state[final.rows]
    }
    final.rows <- which(final.rows)
    time <- time[final.rows]; subject <- subject[final.rows]
### Check if any individuals have only one observation (after excluding missing data)
### Note this shouldn't happen after 1.2
    subj.num <- match(subject,unique(subject)) # avoid problems with factor subjects with empty levels
    nobspt <- table(subj.num)
    if (any (nobspt == 1)) {
        badsubjs <- unique(subject)[ nobspt == 1 ]
        andothers <- if (length(badsubjs)>3) " and others" else ""
        if (length(badsubjs)>3) badsubjs <- badsubjs[1:3]
        badlist <- paste(badsubjs, collapse=", ")
        plural <- if (length(badsubjs)==1) "" else "s"
        has <-  if (length(badsubjs)==1) "has" else "have"
        warning ("Subject", plural, " ", badlist, andothers, " only ", has, " one complete observation, which doesn't give any information")
    }
### Check if observations within a subject are adjacent
    ind <- tapply(seq_along(subj.num), subj.num, length)
    imin <- tapply(seq_along(subj.num), subj.num, min)
    imax <- tapply(seq_along(subj.num), subj.num, max)
    adjacent <- (ind == imax-imin+1)
    if (any (!adjacent)) {
        badsubjs <- unique(subject)[ !adjacent ]
        andothers <- if (length(badsubjs)>3) " and others" else ""
        if (length(badsubjs)>3) badsubjs <- badsubjs[1:3]
        badlist <- paste(badsubjs, collapse=", ")
        plural <- if (length(badsubjs)==1) "" else "s"
        stop ("Observations within subject", plural, " ", badlist, andothers, " are not adjacent in the data")
    }
### Check if observations are ordered in time within subject
    orderedpt <- ! tapply(time, subj.num, is.unsorted)
    if (any (!orderedpt)) {
        badsubjs <- unique(subject)[ !orderedpt ]
        andothers <- if (length(badsubjs)>3) " and others" else ""
        if (length(badsubjs)>3) badsubjs <- badsubjs[1:3]
        badlist <- paste(badsubjs, collapse=", ")
        plural <- if (length(badsubjs)==1) "" else "s"
        stop ("Observations within subject", plural, " ", badlist, andothers, " are not ordered by time")
    }
### Check if any consecutive observations are made at the same time, but with different states
    if (!is.null(state)){
        if (is.matrix(state)) state <- apply(state, 1, paste, collapse=",")
        prevsubj <- c(-Inf, subj.num[seq_along(subj.num)-1])
        prevtime <- c(-Inf, time[1:length(time)-1])
        prevstate <- c(-Inf, state[1:length(state)-1])
        sametime <- final.rows[subj.num==prevsubj & prevtime==time & prevstate!=state]
        badlist <- paste(paste(sametime-1, sametime, sep=" and "), collapse=", ")
        if (any(sametime))
            warning("Different states observed at the same time on the same subject at observations ", badlist)
    }
    invisible()
}

msm.form.obstype <- function(mf, obstype, dmodel, exacttimes)
{
    if (!is.null(obstype)) {
        if (!is.numeric(obstype)) stop("obstype should be numeric")
        if (length(obstype) == 1) obstype <- rep(obstype, nrow(mf))
        else if (length(obstype)!=nrow(mf)) stop("obstype of length ", length(obstype), ", should be length 1 or ",nrow(mf))
        if (any(! obstype[duplicated(mf$"(subject)")] %in% 1:3)) stop("elements of obstype should be 1, 2, or 3") # ignore obstypes at subject's first observation
    }
    else if (!is.null(exacttimes) && exacttimes)
        obstype <- rep(2, nrow(mf))
    else {
        obstype <- rep(1, nrow(mf))
        if (dmodel$ndeath > 0){
            dobs <- if (is.matrix(mf$"(state)")) apply(mf$"(state)", 1, function(x) any(x %in% dmodel$obs)) else (mf$"(state)" %in% dmodel$obs)
            obstype[dobs] <- 3
        }
    }
    obstype
}

### On exit, obstrue will contain the true state (if known) or 0 (if unknown)
### Any NAs should be replaced by 0 - logically if you don't know whether the state is known or not, that means you don't know the state

### If "obstrue" and "censor" both specified, obstrue=1 specifies that the true state is within "censor.states", and we have no other outcome generated conditionally on that true state at that time.  In HMMs with "censor" where obstrue not specified, then it is the observed state which is assumed to be within "censor.states.

msm.form.obstrue <- function(mf, hmodel, cmodel) {
    obstrue <- mf$"(obstrue)"
    if (!is.null(obstrue)) {
        if (!hmodel$hidden) {
            warning("Specified obstrue for a non-hidden model, ignoring.")
            obstrue <- rep(1, nrow(mf))
        }
        else if (!is.numeric(obstrue) && !is.logical(obstrue)) stop("obstrue should be logical or numeric")
        else {
            if (is.logical(obstrue) || (all(na.omit(obstrue) %in% 0:1) && !any(is.na(obstrue)))){
                ## obstrue is an indicator.  Actual state should then be supplied in the outcome vector
                ## (typically misclassification models)
                ## interpret presence of NAs as indicating true state supplied here
                if (!is.null(ncol(mf$"(state)")) && ncol(mf$"(state)") > 1)
                    stop("obstrue must contain NA or the true state for a multiple outcome HMM, not an 0/1 indicator")
                obstrue <- ifelse(obstrue, mf$"(state)", 0)
            } else {
                ## obstrue contains the actual state (used when we have another outcome conditionally on this)
                if (!all(na.omit(obstrue) %in% 0:hmodel$nstates)){
                    stop("Interpreting \"obstrue\" as containing true states, but it contains values not in 0,1,...,", hmodel$nstates)
                }
                obstrue[is.na(obstrue)] <- 0 # true state assumed unknown if NA

                ## If misclassification model specified through "ematrix", interpret the state data as the true state.

                ## If specified through hmodel (including hmmCat), interpret the state as a misclassified/HMM outcome generated conditionally on the true state.
            }
        }
    }
    else if (hmodel$hidden) obstrue <- rep(0, nrow(mf))
    else obstrue <- mf$"(state)"
    if (cmodel$ncens > 0){
        ## If censoring and obstrue, put the first of the possible states into obstrue
        ## Used in Viterbi
        for (i in seq_along(cmodel$censor))
            obstrue[obstrue==cmodel$censor[i] & obstrue > 0] <- cmodel$states[cmodel$index[i]]
    }
    obstrue
}

## Replace censored states by state with highest probability that they
## could represent. Used in msm.check.model to check consistency of
## data with transition parameters

msm.impute.censored <- function(fromstate, tostate, Pmat, cmodel)
{
    ## e.g. cmodel$censor 99,999;  cmodel$states 1,2,1,2,3;  cmodel$index 1, 3, 6
    ## Both from and to are censored
    wb <- which ( fromstate %in% cmodel$censor & tostate %in% cmodel$censor)
    for (i in wb) {
        si <- which(cmodel$censor==fromstate[i])
        fc <- cmodel$states[(cmodel$index[si]) : (cmodel$index[si+1]-1)]
        ti <- which(cmodel$censor==tostate[i])
        tc <- cmodel$states[(cmodel$index[ti]) : (cmodel$index[ti+1]-1)]
        mp <- which.max(Pmat[fc, tc])
        fromstate[i] <- fc[row(Pmat[fc, tc])[mp]]
        tostate[i] <- tc[col(Pmat[fc, tc])[mp]]
    }
    ## Only from is censored
    wb <- which(fromstate %in% cmodel$censor)
    for (i in wb) {
        si <- which(cmodel$censor==fromstate[i])
        fc <- cmodel$states[(cmodel$index[si]) : (cmodel$index[si+1]-1)]
        fromstate[i] <- fc[which.max(Pmat[fc, tostate[i]])]
    }
    ## Only to is censored
    wb <- which(tostate %in% cmodel$censor)
    for (i in wb) {
        si <- which(cmodel$censor==tostate[i])
        tc <- cmodel$states[(cmodel$index[si]) : (cmodel$index[si+1]-1)]
        tostate[i] <- tc[which.max(Pmat[fromstate[i], tc])]
    }
    list(fromstate=fromstate, tostate=tostate)
}

### CHECK THAT TRANSITION PROBABILITIES FOR DATA ARE ALL NON-ZERO
### (e.g. check for backwards transitions when the model is irreversible)
### obstype 1 must have unitprob > 0
### obstype 2 must have qunit != 0, and unitprob > 0.
### obstype 3 must have unitprob > 0

msm.check.model <- function(fromstate, tostate, obs, subject, obstype=NULL, qmatrix, cmodel)
{
    n <- length(fromstate)
    qmatrix <- qmatrix / mean(qmatrix[qmatrix>0]) # rescale to avoid false warnings with small rates
    Pmat <- MatrixExp(qmatrix)
    Pmat[Pmat < .Machine$double.eps] <- 0
    imputed <- msm.impute.censored(fromstate, tostate, Pmat, cmodel)
    fs <- imputed$fromstate; ts <- imputed$tostate
    unitprob <- apply(cbind(fs, ts), 1, function(x) { Pmat[x[1], x[2]] } )
    qunit <- apply(cbind(fs, ts), 1, function(x) { qmatrix[x[1], x[2]] } )

    if (identical(all.equal(min(unitprob, na.rm=TRUE), 0),  TRUE))
    {
        badobs <- which.min(unitprob)
        warning ("Data may be inconsistent with transition matrix for model without misclassification:\n",
                 "individual ", if(is.null(subject)) "" else subject[badobs], " moves from state ", fromstate[badobs],
                 " to state ", tostate[badobs], " at observation ", obs[badobs], "\n")
    }
    if (any(qunit[obstype==2]==0)) {
        badobs <- min (obs[qunit==0 & obstype==2], na.rm = TRUE)
        warning ("Data may be inconsistent with intensity matrix for observations with exact transition times and no misclassification:\n",
                 "individual ", if(is.null(subject)) "" else subject[obs==badobs], " moves from state ", fromstate[obs==badobs],
                 " to state ", tostate[obs==badobs], " at observation ", badobs)
    }
    absorbing <- absorbing.msm(qmatrix=qmatrix)
    absabs <- (fromstate %in% absorbing) & (tostate %in% absorbing)
    if (any(absabs)) {
        badobs <- min( obs[absabs] )
        warning("Absorbing - absorbing transition at observation ", badobs)
    }
    invisible()
}

msm.check.constraint <- function(constraint, mm){
    if (is.null(constraint)) return(invisible())
    covlabels <- colnames(mm)[-1]
    if (!is.list(constraint)) stop(deparse(substitute(constraint)), " should be a list")
    if (!all(sapply(constraint, is.numeric)))
        stop(deparse(substitute(constraint)), " should be a list of numeric vectors")
    ## check and parse the list of constraints on covariates
    for (i in names(constraint))
        if (!(is.element(i, covlabels))){
            factor.warn <- if (i %in% names(attr(mm,"contrasts")))
                "\n\tFor factor covariates, specify constraints using covnameCOVVALUE = c(...)"
            else ""
            stop("Covariate \"", i, "\" in constraint statement not in model.", factor.warn)
        }
}

msm.check.covinits <- function(covinits, covlabels){
    if (!is.list(covinits)) warning(deparse(substitute(covinits)), " should be a list")
    else if (!all(sapply(covinits, is.numeric)))
        warning(deparse(substitute(covinits)), " should be a list of numeric vectors")
    else {
        notin <- setdiff(names(covinits), covlabels)
        if (length(notin) > 0){
            plural <- if (length(notin) > 1) "s " else " "
            warning("covariate", plural, paste(notin, collapse=", "), " in ", deparse(substitute(covinits)), " unknown")
        }
    }
}

### Process covariates constraints, in preparation for being passed to the likelihood optimiser
### This function is called for both sets of covariates (transition rates and the misclassification probs)

msm.form.covmodel <- function(mf,
                              mm,
                              constraint,
                              covinits,
                              covmeans,
                              nmatrix,     # number of transition intensities / misclassification probs
                              cri
                              )
{
    if (!is.null(cri))
        return(msm.form.covmodel.byrate(mf, mm, constraint, covinits, covmeans, nmatrix, cri))
    ncovs <- ncol(mm) - 1
    covlabels <- colnames(mm)[-1]
    if (is.null(constraint)) {
        constraint <- rep(list(1:nmatrix), ncovs)
        names(constraint) <- covlabels
        constr <- 1:(nmatrix*ncovs)
    }
    else {
        msm.check.constraint(constraint, mm)
        constr <- inits <- numeric()
        maxc <- 0
        for (i in seq_along(covlabels)){
            ## build complete vectorised list of constraints for covariates in covariates statement
            ## so. e.g. constraints = (x1=c(3,3,4,4,5), x2 = (0.1,0.2,0.3,0.4,0.4))
            ##     turns into constr = c(1,1,2,2,3,4,5,6,7,7) with seven distinct covariate effects
            ## Allow constraints such as: some elements are minus others. Use negative elements of constr to do this.
            ## e.g. constr = c(1,1,-1,-1,2,3,4,5)
            ## obtained by match(abs(x), unique(abs(x))) * sign(x)
            if (is.element(covlabels[i], names(constraint))) {
                if (length(constraint[[covlabels[i]]]) != nmatrix)
                    stop("\"",covlabels[i],"\" constraint of length ",
                         length(constraint[[covlabels[i]]]),", should be ",nmatrix)
            }
            else
                constraint[[covlabels[i]]] <- seq(nmatrix)
            constr <- c(constr, (maxc + match(abs(constraint[[covlabels[i]]]),
                                             unique(abs(constraint[[covlabels[i]]]))))*sign(constraint[[covlabels[i]]]) )
            maxc <- max(abs(constr))
        }
    }
    inits <- numeric()
    if (!is.null(covinits))
        msm.check.covinits(covinits, covlabels)
    for (i in seq_along(covlabels)) {
        if (!is.null(covinits) && is.element(covlabels[i], names(covinits))) {
            thisinit <- covinits[[covlabels[i]]]
            if (!is.numeric(thisinit)) {
                warning("initial values for covariates should be numeric, ignoring")
                thisinit <- rep(0, nmatrix)
            }
            if (length(thisinit) != nmatrix) {
                warning("\"", covlabels[i], "\" initial values of length ", length(thisinit), ", should be ", nmatrix, ", ignoring")
                thisinit <- rep(0, nmatrix)
            }
            inits <- c(inits, thisinit)
        }
        else {
            inits <- c(inits, rep(0, nmatrix))
        }
    }
    npars <- ncovs*nmatrix
    ndpars <- max(unique(abs(constr)))
    list(npars=npars,
         ndpars=ndpars,    # number of distinct covariate effect parameters
         ncovs=ncovs,
         constr=constr,
         covlabels=covlabels, # factors as separate contrasts
         inits = inits,
         covmeans = attr(mm, "means")
         )
}

## Process constraints and initial values for covariates supplied as a
## list of transition-specific formulae.  Convert to form needed for a
## single covariates formula common to all transitions, which can be
## processed with msm.form.covmodel.

msm.form.covmodel.byrate <- function(mf, mm,
                                     constraint, # as supplied by user
                                     covinits,  # as supplied by user
                                     covmeans,
                                     nmatrix,
                                     cri
                                     ){
    covs <- colnames(mm)[-1]
    ## Convert short form constraints to long form
    msm.check.constraint(constraint, mm)
    constr <- inits <- numeric()
    for (i in seq_along(covs)){
        if (covs[i] %in% names(constraint)){
            if (length(constraint[[covs[i]]]) != sum(cri[,i]))
                stop("\"",covs[i],"\" constraint of length ",
                     length(constraint[[covs[i]]]),", should be ",sum(cri[,i]))
            con <- match(constraint[[covs[i]]], unique(constraint[[covs[i]]])) + 1
            constraint[[covs[i]]] <- rep(1, nmatrix)
            constraint[[covs[i]]][cri[,i]==1] <- con
        }
        else constraint[[covs[i]]] <- seq(length=nmatrix)
    }
    ## convert short to long initial values in the same way
    if (!is.null(covinits)) msm.check.covinits(covinits, covs)
    for (i in seq_along(covs)) {
        if (!is.null(covinits) && (covs[i] %in% names(covinits))) {
            if (!is.numeric(covinits[[covs[i]]])) {
                warning("initial values for covariates should be numeric, ignoring")
                covinits[[covs[i]]] <- rep(0, nmatrix)
            }
            thisinit <- rep(0, nmatrix)
            if (length(covinits[[covs[i]]]) != sum(cri[,i])) {
                warning("\"", covs[i], "\" initial values of length ", length(covinits[[covs[i]]]), ", should be ", sum(cri[,i]), ", ignoring")
                covinits[[covs[i]]] <- rep(0, nmatrix)
            }
            else thisinit[cri[,i]==1] <- covinits[[covs[i]]]
            covinits[[covs[i]]] <- thisinit
        }
    }
    qcmodel <- msm.form.covmodel(mf, mm, constraint, covinits, covmeans, nmatrix, cri=NULL)
    qcmodel$cri <- cri
    qcmodel
}


msm.form.dmodel <- function(death, qmodel, hmodel)
{
    nstates <- qmodel$nstates
    statelist <- if (nstates==2) "1, 2" else if (nstates==3) "1, 2, 3" else paste("1, 2, ... ,",nstates)
    if (is.null(death)) death <- FALSE
    if (is.logical(death) && death==TRUE)
        states <- nstates
    else if (is.logical(death) && death==FALSE)
        states <- numeric(0) ## Will be changed to -1 when passing to C
    else if (!is.numeric(death)) stop("Exact death states indicator must be numeric")
    else if (length(setdiff(death, 1:nstates)) > 0)
        stop("Exact death states indicator contains states not in ",statelist)
    else states <- death
    ndeath <- length(states)
    if (hmodel$hidden) {
        ## Form death state info from hmmIdent parameters.
        ## Special observations in outcome data which denote death states
        ## are given as the parameter to hmmIdent()
        mods <- if(is.matrix(hmodel$models)) hmodel$models[1,] else hmodel$models
        if (!all(mods[states] == match("identity", .msm.HMODELS)))
            stop("States specified in \"deathexact\" should have the identity hidden distribution hmmIdent()")
        obs <- ifelse(hmodel$npars[states]>0, 
                      hmodel$pars[hmodel$parstate %in% states],
                      states)
    }
    else obs <- states
    if (any (states %in% transient.msm(qmatrix=qmodel$qmatrix)))
        stop("Not all the states specified in \"deathexact\" are absorbing")
    list(ndeath=ndeath, states=states, obs=obs)
}

msm.form.cmodel <- function(censor=NULL, censor.states=NULL, qmatrix, hmodel=NULL)
{
    if (is.null(censor)) {
        ncens <- 0
        if (!is.null(censor.states)) warning("censor.states supplied but censor not supplied")
    }
    else {
        if (isTRUE(hmodel$mv)) stop("`censor` not supported with multivariate hidden Markov models")
        if (!is.numeric(censor)) stop("censor must be numeric")
        if (any(censor %in% 1:nrow(qmatrix))) warning("some censoring indicators are the same as actual states")
        ncens <- length(censor)
        if (is.null(censor.states)) {
            if (ncens > 1) {
                warning("more than one type of censoring given, but censor.states not supplied. Assuming only one type of censoring")
                ncens <- 1; censor <- censor[1]
            }
            censor.states <- transient.msm(qmatrix=qmatrix)
            states.index <- c(1, length(censor.states)+1)
            states_list <- list(censor.states)
        }
        else {
            if (ncens == 1) {
                if (!is.vector(censor.states) ||
                    (is.list(censor.states) && (length(censor.states) > 1)) )
                    stop("if one type of censoring, censor.states should be a vector, or a list with one vector element")
                if (!is.numeric(unlist(censor.states))) stop("censor.states should be all numeric")
                states.index <- c(1, length(unlist(censor.states))+1)
            }
            else {
                if (!is.list(censor.states)) stop("censor.states should be a list")
                if (length(censor.states) != ncens) stop("expected ", ncens, " elements in censor.states list, found ", length(censor.states))
                states.index <- cumsum(c(0, lapply(censor.states, length))) + 1
            }
            states_list <- censor.states
            censor.states <- unlist(states_list)
        }
        names(states_list) <- censor
    }
    if (ncens==0) censor <- censor.states <- states_list <- states.index <- NULL
    ## Censoring information to be passed to C
    list(ncens = ncens, # number of censoring states
         censor = censor, # vector of their labels in the data
         states = censor.states, # possible true states that the censoring represents (in unlisted form, for C) 
         states_list = states_list, # (same in list form, for R)
         index = states.index # index into censor.states for the start of each true-state set, including an extra length(censor.states)+1.  For C
         )
}

### Transform set of sets of probs {prs} to {log(prs/pr1)}
msm.mnlogit.transform <- function(pars, hmodel){
    res <- pars
    plabs <- hmodel$plabs
    states <- hmodel$parstate
    outcomes <- hmodel$parout
    if (any(plabs=="p")) {
        for (i in unique(states)) {
            for (j in unique(outcomes)) {
                pfree <- which(plabs=="p" & states==i & outcomes==j)
                pbase <- which(plabs=="pbase" & states==i & outcomes==j)
                ## recalculate baseline prob if necessary, e.g. if constraints applied
                res[pbase] <- 1 - sum(res[pfree])
                res[pfree] <- log(res[pfree] / res[pbase])
            }
        }
    }
    res
}

### Transform set of sets of murs = {log(prs/pr1)} to probs {prs}
### ie psum = sum(exp(mus)),  pr1 = 1 / (1 + psum),  prs = exp(mus) / (1 + psum)
msm.mninvlogit.transform <- function(pars, hmodel) {
    res <- pars
    plabs <- hmodel$plabs
    states <- hmodel$parstate
    outcomes <- hmodel$parout 
    if (any(plabs=="p")) { # p's are unconstrained probabilities 
        ## indicator for which p's belong to which state 
        whichst <- match(states[plabs=="p"], unique(states[plabs=="p"]))
        ## indicator for which p's belong to which outcome in multivariate HMMs
        whichout <- match(outcomes[plabs=="p"], unique(outcomes[plabs=="p"]))
        for (i in unique(whichst)) {
            for (j in unique(whichout)) { 
                pfree <- which(plabs=="p")[whichst==i & whichout==j] 
                pbase <- which(plabs=="pbase" & states==i & outcomes==j)

                if (is.matrix(pars)) {# will be used when applying covariates
                    psum <- colSums(exp(pars[pfree,,drop=FALSE]))
                    res[pbase,] <- 1 / (1 + psum)
                    res[pfree,] <- exp(pars[pfree,,drop=FALSE]) /
                      rep(1 + psum, each=sum(whichst==i & whichout==j))
                } else {
                    psum <- sum(exp(pars[pfree]))
                    res[pbase] <- 1 / (1 + psum)
                    res[pfree] <- exp(pars[pfree]) / (1 + psum)
                }
            }
        }
    }
    res
}

## transform parameters from natural scale to real-line optimisation scale

msm.transform <- function(pars, hmodel, ranges){
    labs <- names(pars)
    pars <- glogit(pars, ranges[,"lower"], ranges[,"upper"])
    hpinds <- which(!(labs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov")))
    hpars <- pars[hpinds]
    hpars <- msm.mnlogit.transform(hpars, hmodel)
    pars[hpinds] <- hpars
    pars[labs=="initp"] <- log(pars[labs=="initp"] / pars[labs=="initpbase"])
    pars
}

## transform parameters from real-line optimisation scale to natural scale

msm.inv.transform <- function(pars, hmodel, ranges){
    labs <- names(pars)
    pars <- gexpit(pars, ranges[,"lower"], ranges[,"upper"])
    hpinds <- which(!(labs %in% c("qbase","qcov","hcov","initp","initp0","initpcov")))
    hpars <- pars[hpinds]
    hpars <- msm.mninvlogit.transform(hpars, hmodel)
    pars[hpinds] <- hpars
    ep <- exp(pars[labs=="initp"])
    pars[labs=="initp"] <- ep / (1 + sum(ep))
    pars[labs=="initpbase"] <- 1 / (1 + sum(ep))
    pars
}

## Collect all model parameters together ready for optimisation
## Handle parameters fixed at initial values or constrained to equal other parameters

msm.form.params <- function(qmodel, qcmodel, emodel, hmodel, fixedpars)
{
    ## Transition intensities
    ni <- qmodel$npars
    ## Covariates on transition intensities
    nc <- qcmodel$npars
    ## HMM response parameters
    nh <- sum(hmodel$npars)
    ## Covariates on HMM response distribution
    nhc <- sum(hmodel$ncoveffs)
    ## Initial state occupancy probabilities in HMM
    nip <- hmodel$nipars
    ## Covariates on initial state occupancy probabilities.
    nipc <- hmodel$nicoveffs
    npars <- ni + nc + nh + nhc + nip + nipc
    inits <- as.numeric(c(qmodel$inits, qcmodel$inits, hmodel$pars, unlist(hmodel$coveffect)))
    plabs <- c(rep("qbase",ni), rep("qcov", nc), hmodel$plabs, rep("hcov", nhc))
    if (nip > 0) {
        inits <- c(inits, hmodel$initprobs)
        initplabs <- c("initpbase", rep("initp",nip-1))
        initplabs[hmodel$initprobs==0] <- "initp0" # those initialised to zero will be fixed at zero
        plabs <- c(plabs, initplabs)
        if (nipc > 0) {
            inits <- c(inits, unlist(hmodel$icoveffect))
            plabs <- c(plabs, rep("initpcov",nipc))
        }
    }
    ## store indicator for which parameters are HMM location parameters (not HMM cov effects or initial state probs)
    hmmpars <- which(!(plabs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov")))
    hmmparscov <- which(!(plabs %in% c("qbase","qcov","initpbase","initp","initp0","initpcov")))
    names(inits) <- plabs
    ranges <- .msm.PARRANGES[plabs,,drop=FALSE]
    if (!is.null(hmodel$ranges)) ranges[hmmparscov,] <- hmodel$ranges
    inits <- msm.transform(inits, hmodel, ranges)
    ## Form constraint vector for complete set of parameters
    ## No constraints allowed on initprobs and their covs for the moment
    constr <- c(qmodel$constr,
                if(is.null(qcmodel$constr)) NULL else (ni + abs(qcmodel$constr))*sign(qcmodel$constr),
                ni + nc + hmodel$constr,
                ni + nc + nh + hmodel$covconstr,
                ni + nc + nh + nhc + seq(length.out=nip),
                ni + nc + nh + nhc + nip + seq(length.out=nipc))
    constr <- match(abs(constr), unique(abs(constr)))*sign(constr)
    ## parameters which are always fixed and not included in user-supplied fixedpars
    auxpars <- which(plabs %in% .msm.AUXPARS)
    duppars <- which(duplicated(abs(constr)))
    realpars <- setdiff(seq(npars), union(auxpars, duppars))
    nrealpars <- npars - length(auxpars) - length(duppars)
    ## if transition-specific covariates, then fixedpars indices generally smaller
    nshortpars <- nrealpars - sum(qcmodel$cri[!duplicated(qcmodel$constr)]==0)
    if (is.logical(fixedpars))
        fixedpars <- if (fixedpars == TRUE) seq(nshortpars) else numeric()
    if (any(! (fixedpars %in% seq(length.out=nshortpars))))
        stop ( "Elements of fixedpars should be in 1, ..., ", nshortpars)
    if (!is.null(qcmodel$cri)) {
        ## Convert user-supplied fixedpars indexing transition-specific covariates
        ## to fixedpars indexing transition-common covariates
        inds <- rep(1, nrealpars)
        inds[qmodel$ndpars + qcmodel$constr[!duplicated(qcmodel$constr)]] <-
            qcmodel$cri[!duplicated(qcmodel$constr)]
        inds[inds==1] <- seq(length.out=nshortpars)
        fixedpars <- match(fixedpars, inds)
        ## fix covariate effects not included in model to zero
        fixedpars <- sort(c(fixedpars, which(inds==0)))
    }
    notfixed <- realpars[setdiff(seq_along(realpars),fixedpars)]
    fixedpars <- sort(c(realpars[fixedpars], auxpars)) # change fixedpars to index unconstrained pars
    allinits <- inits
    optpars <- intersect(notfixed, which(!duplicated(abs(constr))))
    inits <- inits[optpars]
    fixed <- (length(fixedpars) + length(duppars) == npars) # TRUE if all parameters are fixed, then no optimisation needed, just evaluate likelihood
    names(allinits) <- plabs; names(fixedpars) <- plabs[fixedpars]; names(plabs) <- NULL
    paramdata <- list(inits=inits, plabs=plabs, allinits=allinits, hmmpars=hmmpars,
                      fixed=fixed, fixedpars=fixedpars,
                      optpars=optpars, auxpars=auxpars,
                      constr=constr, npars=npars, duppars=duppars,
                      nfix=length(fixedpars), nopt=length(optpars), ndup=length(duppars),
                      ranges=ranges)
    paramdata
}

## Unfix all fixed parameters in a paramdata object p (as
## returned by msm.form.params).  Used for calculating deriv /
## information over all parameters for models with fixed parameters.
## Don't unfix auxiliary pars such as binomial denominators
## Don't unconstrain constraints

msm.unfixallparams <- function(p) {
    npars <- length(p$allinits)
    p$fixed <- FALSE
    p$optpars <- setdiff(1:npars, union(p$auxpars,p$duppars))
    p$fixedpars <- p$auxpars
    p$nopt <- length(p$optpars)
    p$nfix <- npars - length(p$optpars)
    p$inits <- p$allinits[p$optpars]
    p
}

msm.rep.constraints <- function(pars, # transformed pars
                                paramdata,
                                hmodel){
    plabs <- names(pars)
    p <- paramdata
    ## Handle constraints on misclassification probs / HMM cat probs separately
    ## This is fiddly. First replicate within states, on log(pr/pbase) scale
    if (any(hmodel$plabs=="p")) {
        for (i in 1:hmodel$nstates){
            inds <- (hmodel$parstate==i & hmodel$plabs=="p")
            hpc <- hmodel$constr[inds]
            pars[p$hmmpars][inds] <- pars[p$hmmpars][inds][match(hpc, unique(hpc))]
        }
    }
    ## ...then replicate them between states, on pr scale, so that constraint applies
    ## to pr not log(pr/pbase). After, transform back
    pars[p$hmmpars] <- msm.mninvlogit.transform(pars[p$hmmpars], hmodel)
    plabs <- plabs[!duplicated(abs(p$constr))][abs(p$constr)]
    pars <- pars[!duplicated(abs(p$constr))][abs(p$constr)]*sign(p$constr)
    pars[p$hmmpars] <- msm.mnlogit.transform(pars[p$hmmpars], hmodel)
    names(pars) <- plabs
    pars
}

## Apply covariates to transition intensities. Parameters enter this
## function already log transformed and replicated, and exit on
## natural scale.

msm.add.qcovs <- function(qmodel, pars, mm){
    labs <- names(pars)
    beta <- rbind(pars[labs=="qbase"],
                  matrix(pars[labs=="qcov"], ncol=qmodel$npars, byrow=TRUE))
    qvec <- exp(mm %*% beta)
    imat <- t(qmodel$imatrix); row <- col(imat)[imat==1]; col <- row(imat)[imat==1]
    qmat <- array(0, dim=c(qmodel$nstates, qmodel$nstates, nrow(mm)))
    for (i in 1:qmodel$npars) {
        qmat[row[i],col[i],] <- qvec[,i]
    }
    for (i in 1:qmodel$nstates)
        ## qmat[i,i,] <- -apply(qmat[i,,,drop=FALSE], 3, sum)
        qmat[i,i,] <- -colSums(qmat[i,,,drop=FALSE], , 2)
    qmat
}

## Derivatives of intensity matrix Q wrt unique log q and beta, after
## applying constraints.  By observation with covariates applied.
## e.g. qmodel$constr 1 1 2, qcmodel$constr 1 -1 2 3 3 -3
## returns derivs w.r.t pars named p1 p2 p3 p4 p5
## i.e. with baseline q on the log scale
## qo = exp(p1 + p3x1 + p5x2)
## q1 = exp(p1 + -p3x1 + p5x2)
## q2 = exp(p2 + p4x1  - p5x2)

msm.form.dq <- function(qmodel, qcmodel, pars, paramdata, mm){
    labs <- names(pars)
    q0 <- pars[labs=="qbase"]
    beta <- rbind(q0, matrix(pars[labs=="qcov"], ncol=qmodel$npars, byrow=TRUE))
    qvec <- exp(mm %*% beta)
    covs <- mm[,-1,drop=FALSE]
    qrvec <- exp(covs %*% beta[-1,,drop=FALSE])
    nopt <- qmodel$ndpars + qcmodel$ndpars # after constraint but before omitting fixed
    dqvec <- array(pars[labs=="qbase"],
                   dim=c(nrow(mm), qmodel$npars, nopt))
    for (i in seq_len(qmodel$npars)) {
        ind <- rep(0, qmodel$ndpars); ind[qmodel$constr[i]] <- 1
        cind <- 1:qmodel$ndpars
#        dqvec[,i,cind] = qrvec[,i] * rep(ind, each=nrow(mm))
        dqvec[,i,cind] = qrvec[,i] * rep(ind, each=nrow(mm)) * exp(q0[i])
        if (qcmodel$npars > 0)
            for (k in 1:qcmodel$ncovs) {
                con <- qcmodel$constr[(k-1)*qmodel$npars + 1:qmodel$npars]
                ucon <- match(abs(con), unique(abs(con)))
                ind <- rep(0, max(ucon)); ind[ucon[i]] <- sign(con[i])
                cind <- max(cind) + unique(ucon)
                dqvec[,i,cind] <- covs[,k] * qvec[,i] * rep(ind, each=nrow(mm))
            }
    }
    dqmat <- array(0, dim=c(qmodel$nstates, qmodel$nstates, nopt, nrow(mm)))
    imat <- t(qmodel$imatrix); row <- col(imat)[imat==1]; col <- row(imat)[imat==1]
    for (i in 1:qmodel$npars)
        dqmat[row[i],col[i],,] <- t(dqvec[,i,])
    for (i in 1:qmodel$nstates)
        ## dqmat[i,i,,] <- -apply(dqmat[i,,,,drop=FALSE], c(3,4), sum)
        dqmat[i,i,,] <- -colSums(dqmat[i,,,,drop=FALSE], , 2)
    p <- paramdata # leave fixed parameters out
    fixed <- p$fixedpars[p$plabs[p$fixedpars] %in% c("qbase","qcov")]
    con <- abs(p$constr[p$plabs %in% c("qbase","qcov")])
    ## FIXME
    if (any(fixed)) dqmat <- dqmat[,,-con[fixed],,drop=FALSE]
    dqmat
}


## Apply covariates to HMM location parameters
## Parameters enter transformed, and exit on natural scale

msm.add.hmmcovs <- function(hmodel, pars, mml){
    labs <- names(pars)
    hpinds <- which(!(labs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov")))
    n <- nrow(mml[[1]])
    hpars <- matrix(rep(pars[hpinds], n), ncol=n,  dimnames=list(labs[hpinds], NULL))
    ito <- 0
    for (i in which(hmodel$ncovs > 0)) {
        mm <- mml[[hmodel$parstate[i]]]
        ## TODO is this cleaner with hmodel$coveffstate ?
        ## TODO is this correct for misc covs: two matrices
        ifrom <- ito + 1; ito <- ito + hmodel$ncovs[i]
        beta <- pars[names(pars)=="hcov"][ifrom:ito]
        hpars[i,] <- hpars[i,] + mm[,-1,drop=FALSE] %*% beta
    }
    for (i in seq_along(hpinds)){
        hpars[i,] <- gexpit(hpars[i,], hmodel$ranges[i,"lower"], hmodel$ranges[i,"upper"])
    }
    hpars <- msm.mninvlogit.transform(hpars, hmodel)
    hpars
}

### Derivatives of HMM pars w.r.t HMM baseline pars (on transformed,
### e.g. log, scale) and covariate effects
## Account for parameter constraints, e.g.
## hmodel$constr 1 2 3 2 4, hmodel$covconstr 1 2 1 3
## mu0 = g(p1 + p5x1 + p6x2);  g is id for mean, g'(p)=1, or exp for sigma
## sd0 = g(p2)
## mu1 = g(p3 + p5x1 + p7x2)
## sd1 = g(p2)
## id1 = g(p4)
## want deriv wrt p1,p2,...,p7

msm.form.dh <- function(hmodel,
                        pars, # before adding covariates, and on transformed scale
                        newpars, # after adding covariates, and on natural scale
                        paramdata,
                        mml){
    labs <- names(pars)
    hpinds <- which(!(labs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov"))) # baseline HMM parameters
    n <- nrow(mml[[1]])
    hpars <- matrix(rep(pars[hpinds], n), ncol=n,  dimnames=list(labs[hpinds], NULL))
    nh <- sum(hmodel$npars)
    nopt <- length(unique(hmodel$constr)) + length(unique(hmodel$covconstr)) # TODO store in hmodel
    dh <- array(0, dim=c(nh, nopt, n))
    for (i in 1:nh){
        ind <- rep(0, length(unique(hmodel$constr)))
        ind[hmodel$constr[i]] <- 1
        cind <- 1:length(unique(hmodel$constr))
        if (labs[hpinds][i] != "p") {
            a <- hmodel$ranges[i,"lower"]; b <- hmodel$ranges[i,"upper"]
            gdash <- dgexpit(glogit(newpars[i,], a, b), a, b)
            dh[i,cind,] <- rep(ind, n)*rep(gdash, each=length(cind))
            if (hmodel$ncovs[i] > 0) {
                mm <- mml[[hmodel$parstate[i]]][,-1,drop=FALSE]
                con <- hmodel$covconstr[hmodel$coveffstate == hmodel$parstate[i]]
                cind <- max(cind) + unique(con)
                dh[i,cind,] <- t(mm) * rep(gdash, each=length(cind))
            }
        }
    }
    dh <- msm.dmninvlogit(hmodel, hpars, mml, pars[names(pars)=="hcov"], dh)
    p <- paramdata # leave fixed parameters out
    dhinds <- which(!(labs %in% c("qbase","qcov","initpbase","initp","initp0","initpcov")))
    fixed <- intersect(p$fixedpars, dhinds) # index into full par vec
    fixed <- match(fixed, setdiff(dhinds, p$duppars)) # index into constrained hmm par vec
    if (any(fixed)) dh <- dh[,-fixed,,drop=FALSE]
    dh
}

### Derivatives of outcome probabilities w.r.t baseline log odds and covariate effects
### in a HMM categorical outcome distribution
### Constraints on probabilities not supported, seems too fiddly.
### Could generalize for use in est.initprobs, if ever support derivatives for that

### log(p_r/pbase) = lp_r + beta_r ' x
### p_r / pbase  =  exp(lp_r + beta_r'x)
### (p1 +...+ pR) / pbase  = 1/pbase = sum_allR(exp()), so
### pbase = 1/sum_allR(exp())
### p_r = exp(lp_r + beta_r'x) / sum_all(exp(lp_s + beta_s'x)), where pars 0 for s=pbase
### d/dlp_s pbase    = -exp(lp_s + beta_s'x)/ (sum_all())^2
### d/dbeta_sk pbase = (-x_sk * exp(lp_s + beta_s'x)/ (sum_all())^2

### for r!=s, d/dlp_s p_r  =  exp(lp_r + beta_r'x) * -exp(lp_s + beta_s'x)/ (sum_all())^2
### for r=s,  d/dlp_s p_r  =  exp(lp_r + beta_r'x) / sum_all()  +
###                                   exp(lp_r + beta_r'x) * -exp(lp_r + beta_r'x)/ (sum_all())^2
### for r!=s, d/dbeta_sk p_r  =  exp(lp_r + beta_r'x) * (-x_sk * exp(lp_s + beta_s'x)/ (sum_all())^2
### for r=s,  d/dbeta_sk p_r  =  x_sk exp(lp_r + beta_r'x) / sum_all()  +
###                                   exp(lp_r + beta_r'x) * (-x_sk * exp(lp_r + beta_r'x)/ (sum_all())^2

### if no covs:  p_r = exp(lp_r) / sum_all(exp(lp_s)).   lp_r = log(p_r/pbase)
### d/dlp_s p_r =
### for r!=s,  exp(lp_r) * (-exp(lp_s)/ (sum_all())^2
### for r=s,  d/dlp_s p_r  =  exp(lp_r) / sum_all() - exp(lp_r)^2 / (sum_all())^2

msm.dmninvlogit <- function(hmodel, pars, mml, hcov, dh){
    plabs <- hmodel$plabs
    states <- hmodel$parstate
    n <- nrow(mml[[1]])
    nh <- sum(hmodel$npars)
    if (any(plabs=="p")) {
        whichst <- match(states[plabs=="p"], unique(states[plabs=="p"]))
        for (i in unique(whichst)) {
            labsi <- plabs[states==i]
            ppars <- which(labsi %in% c("p","pbase","p0"))
            lp <- matrix(pars[states==i ,1][ppars], nrow=n, ncol=length(ppars),
                         byrow=TRUE, dimnames=list(NULL,labsi[ppars]))
            ncovs <- hmodel$ncovs[states==i][ppars]
            beta <- hcov[hmodel$coveffstate==i]
            X <- mml[[states[i]]][,-1,drop=FALSE]
            ito <- 0
            for (j in which(ncovs > 0)){
                ifrom <- ito + 1; ito <- ito + ncovs[j]
                betaj <- beta[ifrom:ito]
                lp[,j] <- lp[,j] + X %*% betaj
            }
            elp <- lp
            elp[,labsi[ppars]=="p"] <- exp(elp[,labsi[ppars]=="p"])
            psum <- 1 + rowSums(elp[,labsi[ppars]=="p",drop=FALSE])
            pil <- which(states==i & plabs=="p")
            pis <- which(colnames(elp) == "p")
            for (r in seq_along(pil)){
                ito <- 0
                for (s in seq_along(pil)){
                    dh[pil[r],pil[s],] <-
                        if (r == s)
                            elp[,pis[r]] / psum * (1 - elp[,pis[r]] / psum)
                        else
                            - elp[,pis[r]] * elp[,pis[s]] / psum^2

                    if (ncovs[pis[s]] > 0){
                        ifrom <- ito + 1; ito <- ito + ncovs[pis[s]]
                        bil <- nh + which(hmodel$coveffstate==i)[ifrom:ito]
                        for (k in seq_along(bil))
                            dh[pil[r],bil[k],] <-
                                if (r == s)
                                    X[,k] * elp[,pis[r]] / psum * (1 - elp[,pis[r]] / psum)
                                else
                                    - X[,k] *  elp[,pis[r]] * elp[,pis[s]] / psum^2
                    }
                }
            }
            pib <- which(states==i & plabs=="pbase")
            ito <- 0
            for (s in seq_along(pil)){
                dh[pib,pil[s],] <- - elp[,pis[s]] / psum^2
                if (ncovs[pis[s]] > 0){
                    ifrom <- ito + 1; ito <- ito + ncovs[pis[s]]
                    bil <- nh + which(hmodel$coveffstate==i)[ifrom:ito]
                    for (k in seq_along(bil))
                        dh[pib,bil[k],] <- - X[,k] * elp[,pis[s]] / psum^2
                }
            }
        }
    }
    dh
}

msm.initprobs2mat <- function(hmodel, pars, mm, mf, cmodel){
    npts <- attr(mf, "npts")
    ## Convert vector initial state occupancy probs to matrix by patient
    if (!hmodel$hidden) return(0)
    if (hmodel$est.initprobs) {
        initp <- pars[names(pars) %in% c("initpbase","initp","initp0")]
        initp <- matrix(rep(initp, each=npts), nrow=npts,
                        dimnames=list(NULL,
                        names(pars)[names(pars) %in% c("initpbase","initp","initp0")]))
        ## Multiply baselines (entering on mnlogit scale) by current covariate
        ## effects, giving matrix of patient-specific initprobs
        est <- which(colnames(initp)=="initp")
        ip <- initp[,est,drop=FALSE]
        if (hmodel$nicoveffs > 0) {
            ## cov effs ordered by states (excluding state 1) within covariates
            coveffs <- pars[names(pars)=="initpcov"]
            coveffs <- matrix(coveffs, nrow=max(hmodel$nicovs), byrow=TRUE)
            ip <- ip + as.matrix(mm[,-1,drop=FALSE]) %*% coveffs
        }
        initp[,est] <- exp(ip) / (1 + rowSums(exp(ip)))
        initp[,"initpbase"] <- 1 / (1 + rowSums(exp(ip)))
    }
    else if (!is.matrix(hmodel$initprobs))
        initp <- matrix(rep(hmodel$initprobs,each=npts),nrow=npts)
    else initp <- hmodel$initprobs

    ## If the initial state is fully or partially known for some people but not others
    ## can specify this by setting both "censor" and "obstrue" at those times 
    ## Adjust initprobs in those cases so it is zero for states outside the censor set
    ## and reweight other entries to sum to 1. 
    initstate <- mf$"(state)"[!duplicated(mf$"(subject)")]
    initobstrue <- mf$"(obstrue)"[!duplicated(mf$"(subject)")]
    initknown <- (initstate %in% cmodel$censor) & (initobstrue != 0) 
    if (cmodel$ncens > 0 && any(initknown)) {
        for (i in 1:npts) {
            if (initknown[i]){
                cs <- cmodel$states_list[[as.character(initstate[i])]]
                ip <- initp[i,cs]
                initp[i,] <- 0
                initp[i,cs] <- ip / sum(ip)
            }
        }
    }
    initp
}

## Entry point to C code for calculating the likelihood and related quantities

Ccall.msm <- function(params, do.what="lik", msmdata, qmodel, qcmodel, cmodel, hmodel, paramdata)
{
    p <- paramdata
    pars <- p$allinits
    pars[p$optpars] <- params
    pars <- msm.rep.constraints(pars, paramdata, hmodel)
    agg <- if (!hmodel$hidden && cmodel$ncens==0 &&
               do.what %in% c("lik","deriv","info")) TRUE else FALSE # data as aggregate transition counts, as opposed to individual observations
    ## Add covariates to hpars and q here. Inverse-transformed to natural scale on exit
    mm.cov <- if (agg) msmdata$mm.cov.agg else msmdata$mm.cov
    Q <- msm.add.qcovs(qmodel, pars, mm.cov)
    DQ <- if (do.what %in% c("deriv","info","deriv.subj","dpmat")) msm.form.dq(qmodel, qcmodel, pars, p, mm.cov) else NULL
    H <- if (hmodel$hidden) msm.add.hmmcovs(hmodel, pars, msmdata$mm.hcov) else NULL
    DH <- if (hmodel$hidden && (do.what %in% c("deriv","info","deriv.subj"))) msm.form.dh(hmodel, pars, H, paramdata, msmdata$mm.hcov) else NULL
    initprobs <- msm.initprobs2mat(hmodel, pars, msmdata$mm.icov, msmdata$mf, cmodel)

   mf <- msmdata$mf; mf.agg <- msmdata$mf.agg
   ## In R, ordinal variables indexed from 1.  In C, these are indexed from 0.
   mf.agg$"(fromstate)" <- mf.agg$"(fromstate)" - 1
   mf.agg$"(tostate)" <- mf.agg$"(tostate)" - 1
   firstobs <- c(which(!duplicated(model.extract(mf, "subject"))), nrow(mf)+1) - 1
   mf$"(subject)" <- match(mf$"(subject)", unique(mf$"(subject)"))
   ntrans <- sum(duplicated(model.extract(mf, "subject")))
   hmodel$models <- hmodel$models - 1
   nagg <- if(is.null(mf.agg)) 0 else nrow(mf.agg)
   mf$"(pcomb)" <- mf$"(pcomb)" - 1
   npcombs <- length(unique(na.omit(model.extract(mf, "pcomb"))))
   qmodel$nopt <- if (is.null(DQ)) 0 else dim(DQ)[3]
   hmodel$nopt <- if (is.null(DH)) 0 else dim(DH)[2]
   nopt <- qmodel$nopt + hmodel$nopt

   ## coerce types here to avoid PROTECT faff with doing this in C
   mfac <- list("(fromstate)" = as.integer(mf.agg$"(fromstate)"),
                  "(tostate)" = as.integer(mf.agg$"(tostate)"),
                  "(timelag)" = as.double(mf.agg$"(timelag)"),
                  "(nocc)" = as.integer(mf.agg$"(nocc)"),
                  "(noccsum)" = as.integer(mf.agg$"(noccsum)"),
                  "(whicha)" = as.integer(mf.agg$"(whicha)"),
                  "(obstype)" = as.integer(mf.agg$"(obstype)"))
   mfc <- list("(subject)" = as.integer(mf$"(subject)"),  "(time)" = as.double(mf$"(time)"),
               ## supply matrix outcomes to C by row so multivariate outcomes are together
               ## Also, state data for misclassification HMMs are indexed from 1 not 0 in C
               "(state)" = as.double(t(mf$"(state)")), 
               "(obstype)" = as.integer(mf$"(obstype)"),
              "(obstrue)" = as.integer(mf$"(obstrue)"), "(pcomb)" = as.integer(mf$"(pcomb)"))
    auxdata <- list(nagg=as.integer(nagg),n=as.integer(nrow(mf)),npts=as.integer(attr(mf,"npts")),
                   ntrans=as.integer(ntrans), npcombs=as.integer(npcombs),
                   nout = as.integer(if(is.null(ncol(mf$"(state)"))) 1 else ncol(mf$"(state)")),
                   nliks=as.integer(get("nliks",msm.globals)),firstobs=as.integer(firstobs))
   qmodel <- list(nstates=as.integer(qmodel$nstates), npars=as.integer(qmodel$npars),
                  nopt=as.integer(qmodel$nopt), iso=as.integer(qmodel$iso),
                  perm=as.integer(qmodel$perm), qperm=as.integer(qmodel$qperm),
                  expm=as.integer(qmodel$expm))
   cmodel <- list(ncens=as.integer(cmodel$ncens), censor=as.integer(cmodel$censor),
                  states=as.integer(cmodel$states), index=as.integer(cmodel$index - 1))
   hmodel <- list(hidden=as.integer(hmodel$hidden), mv=as.integer(hmodel$mv),
                  models=as.integer(hmodel$models),
                  totpars=as.integer(hmodel$totpars), firstpar=as.integer(hmodel$firstpar),
                  npars=as.integer(hmodel$npars), nopt=as.integer(hmodel$nopt),
                  ematrix=as.integer(hmodel$ematrix))
    pars <- list(Q=as.double(Q),DQ=as.double(DQ),H=as.double(H),DH=as.double(DH),
                initprobs=as.double(initprobs),nopt=as.integer(nopt))
    .Call("msmCEntry",  as.integer(match(do.what, .msm.CTASKS) - 1),
         mfac, mfc, auxdata, qmodel, cmodel, hmodel, pars, PACKAGE="msm")
}


lik.msm <- function(params, ...)
{
    ## number of likelihood evaluations so far including this one
    ## used for error message for iffy initial values in HMMs
    assign("nliks", get("nliks",msm.globals) + 1, envir=msm.globals)
    Ccall.msm(params, do.what="lik", ...)
}

deriv.msm <- function(params, ...)
{
    Ccall.msm(params, do.what="deriv", ...)
}

information.msm <- function(params, ...)
{
    Ccall.msm(params, do.what="info", ...)
}

## Convert vector of MLEs into matrices and append them to the model object

msm.form.output <- function(x, whichp)
{
    model <- if (whichp=="intens") x$qmodel else x$emodel
    cmodel <- if (whichp=="intens") x$qcmodel else x$ecmodel
    p <- x$paramdata
    Matrices <- MatricesSE <- MatricesL <- MatricesU <- MatricesFixed <- list()
    basename <- if (whichp=="intens") "logbaseline" else "logitbaseline"
    fixedpars.logical <- p$constr %in% p$constr[p$fixedpars]
    for (i in 0:cmodel$ncovs) {
        matrixname <- if (i==0) basename else cmodel$covlabels[i] # name of the current output matrix.
        mat <- t(model$imatrix) # state matrices filled by row, while R fills them by column.
        if (whichp=="intens")
            parinds <- if (i==0) which(p$plabs=="qbase") else which(p$plabs=="qcov")[(i-1)*model$npars + 1:model$npars]
        if (whichp=="misc")
            parinds <- if (i==0) which(p$plabs=="p") else which(p$plabs=="hcov")[i + cmodel$ncovs*(1:model$npars - 1)]
        if (any(parinds)) mat[t(model$imatrix)==1] <- p$params[parinds]
        else mat[mat==1] <- Inf ## if no parinds are "p", then there are off-diag 1s in ematrix
        mat <- t(mat)
        dimnames(mat) <- dimnames(model$imatrix)
        fixed <- array(FALSE, dim=dim(model$imatrix))
        if (p$foundse && !p$fixed){
            intenscov <- p$covmat[parinds, parinds]
            intensse <- sqrt(diag(as.matrix(intenscov)))
            semat <- lmat <- umat <- t(model$imatrix)
            if (any(parinds)){
                semat[t(model$imatrix)==1] <- intensse
                lmat[t(model$imatrix)==1] <- p$ci[parinds,1]
                umat[t(model$imatrix)==1] <- p$ci[parinds,2]
                fixed[t(model$imatrix)==1] <- fixedpars.logical[parinds]
            }
            else semat[semat==1] <- lmat[lmat==1] <- umat[umat==1] <- Inf
            semat <- t(semat); lmat <- t(lmat); umat <- t(umat); fixed <- t(fixed)
            diag(semat) <- diag(lmat) <- diag(umat) <- 0
            for (i in 1:nrow(fixed)){
                foff <- fixed[i,-i][model$imatrix[i,-i]==1]
                fixed[i,i] <- (length(foff)>1) && all(foff)
            }
            if (whichp=="misc")
                fixed[which(x$hmodel$model==match("identity", .msm.HMODELS)),] <- TRUE
           
            dimnames(semat)  <- dimnames(mat)
        }
        else {
            semat <- lmat <- umat <- NULL
        }
        Matrices[[matrixname]] <- mat
        MatricesSE[[matrixname]] <- semat
        MatricesL[[matrixname]] <- lmat
        MatricesU[[matrixname]] <- umat
        MatricesFixed[[matrixname]] <- fixed
    }
    nam <- if(whichp=="intens") "Qmatrices" else "Ematrices"
    x[[nam]] <- Matrices; x[[paste0(nam, "SE")]] <- MatricesSE; x[[paste0(nam, "L")]] <- MatricesL
    x[[paste0(nam, "U")]] <- MatricesU; x[[paste0(nam, "Fixed")]] <- MatricesFixed
    x
}

## Format hidden Markov model estimates and CIs

msm.form.houtput <- function(hmodel, p, msmdata, cmodel)
{
    hmodel$pars <- p$estimates.t[!(p$plabs %in% c("qbase","qcov","hcov","initp","initpbase","initp0","initpcov"))]
    hmodel$coveffect <- p$estimates.t[p$plabs == "hcov"]
    hmodel$fitted <- !p$fixed
    hmodel$foundse <- p$foundse
    if (hmodel$nip > 0) {
        iplabs <- p$plabs[p$plabs %in% c("initp","initp0")]
        whichst <- which(iplabs == "initp") + 1  # init probs for which states have covs on them (not the zero probs)
        if (hmodel$foundse) {
            hmodel$initprobs <- rbind(cbind(p$estimates.t[p$plabs %in% c("initpbase","initp","initp0")],
                                            p$ci[p$plabs %in% c("initpbase","initp","initp0"),,drop=FALSE]))
            rownames(hmodel$initprobs) <- paste("State",1:hmodel$nstates)
            colnames(hmodel$initprobs) <- c("Estimate", "LCL", "UCL")
            if (any(hmodel$nicovs > 0)) {
                covnames <- names(hmodel$icoveffect)
                hmodel$icoveffect <- cbind(p$estimates.t[p$plabs == "initpcov"],  p$ci[p$plabs == "initpcov",,drop=FALSE])
                rownames(hmodel$icoveffect) <- paste(covnames, paste("State",whichst), sep=", ")
                colnames(hmodel$icoveffect) <- c("Estimate", "LCL", "UCL")
            }
        }
        else {
            hmodel$initprobs <- c(1 - sum(p$estimates.t[p$plabs == "initp"]),
                                  p$estimates.t[p$plabs %in% c("initp","initp0")])
            names(hmodel$initprobs) <- paste("State", 1:hmodel$nstates)
            if (any(hmodel$nicovs > 0)) {
                covnames <- names(hmodel$icoveffect)
                hmodel$icoveffect <- p$estimates.t[p$plabs == "initpcov"]
#                names(hmodel$icoveffect) <- paste(covnames, paste("State",2:hmodel$nstates), sep=", ")
                names(hmodel$icoveffect) <- paste(covnames, paste("State",whichst), sep=", ")
            }
        }
    }
    hmodel$initpmat <- msm.initprobs2mat(hmodel, p$estimates, msmdata$mm.icov, msmdata$mf, cmodel)
    if (hmodel$foundse) {
        hmodel$ci <- p$ci[!(p$plabs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov")), , drop=FALSE]
        hmodel$covci <- p$ci[p$plabs %in% c("hcov"), ]
    }
    names(hmodel$pars) <- hmodel$plabs
    hmodel
}



## Table of 'transitions': previous state versus current state



#' Table of transitions
#' 
#' Calculates a frequency table counting the number of times each pair of
#' states were observed in successive observation times.  This can be a useful
#' way of summarising multi-state data.
#' 
#' If the data are intermittently observed (panel data) this table should not
#' be used to decide what transitions should be allowed in the \eqn{Q} matrix,
#' which works in continuous time.  This function counts the transitions
#' between states over a time interval, not in real time.  There can be
#' observed transitions between state \eqn{r} and \eqn{s} over an interval even
#' if \eqn{q_{rs}=0}, because the process may have passed through one or more
#' intermediate states in the middle of the interval.
#' 
#' @param state Observed states, assumed to be ordered by time within each
#' subject.
#' @param subject Subject identification numbers corresponding to \code{state}.
#' If not given, all observations are assumed to be on the same subject.
#' @param data An optional data frame in which the variables represented by
#' \code{subject} and \code{state} can be found.
#' @return A frequency table with starting states as rows and finishing states
#' as columns.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{crudeinits.msm}}
#' @keywords models
#' @examples
#' 
#' ## Heart transplant data
#' data(cav)
#' 
#' ## 148 deaths from state 1, 48 from state 2 and 55 from state 3.
#' statetable.msm(state, PTNUM, data=cav)
#' 
#' 
#' @export statetable.msm
statetable.msm <- function(state, subject, data=NULL)
{
    if(!is.null(data)) {
        data <- as.data.frame(data)
        state <- eval(substitute(state), data, parent.frame())
    }
    n <- length(state)
    if (!is.null(data))
        subject <-
            if(missing(subject)) rep(1,n) else eval(substitute(subject), data, parent.frame())
    subject <- match(subject, unique(subject))
    prevsubj <- c(NA, subject[1:(n-1)])
    previous <- c(NA, state[1:(n-1)])
    previous[prevsubj!=subject] <- NA
    ntrans <- table(previous, state)
# or simpler as
# ntrans <- table(state[duplicated(subject,fromLast=TRUE)], state[duplicated(subject)])
    names(dimnames(ntrans)) <- c("from", "to")
    ntrans
}

## Calculate crude initial values for transition intensities by assuming observations represent the exact transition times



#' Calculate crude initial values for transition intensities
#' 
#' Calculates crude initial values for transition intensities by assuming that
#' the data represent the exact transition times of the Markov process.
#' 
#' 
#' Suppose we want a crude estimate of the transition intensity
#' \eqn{q_{rs}}{q_rs} from state \eqn{r} to state \eqn{s}.  If we observe
#' \eqn{n_{rs}}{n_rs} transitions from state \eqn{r} to state \eqn{s}, and a
#' total of \eqn{n_r} transitions from state \eqn{r}, then \eqn{q_{rs} / }{q_rs
#' / q_rr}\eqn{ q_{rr}}{q_rs / q_rr} can be estimated by \eqn{n_{rs} /
#' n_r}{n_rs / n_r}. Then, given a total of \eqn{T_r} years spent in state
#' \eqn{r}, the mean sojourn time \eqn{1 / q_{rr}}{1 / q_rr} can be estimated
#' as \eqn{T_r / n_r}.  Thus, \eqn{n_{rs} / T_r}{n_rs / T_r} is a crude
#' estimate of \eqn{q_{rs}}{q_rs}.
#' 
#' If the data do represent the exact transition times of the Markov process,
#' then these are the exact maximum likelihood estimates.
#' 
#' Observed transitions which are incompatible with the given \code{qmatrix}
#' are ignored.  Censored states are ignored.
#' 
#' @param formula A formula giving the vectors containing the observed states
#' and the corresponding observation times. For example,
#' 
#' \code{state ~ time}
#' 
#' Observed states should be in the set \code{1, \dots{}, n}, where \code{n} is
#' the number of states.  Note hidden Markov models are not supported by this
#' function.
#' @param subject Vector of subject identification numbers for the data
#' specified by \code{formula}. If missing, then all observations are assumed
#' to be on the same subject. These must be sorted so that all observations on
#' the same subject are adjacent.
#' @param qmatrix Matrix of indicators for the allowed transitions.  An initial
#' value will be estimated for each value of qmatrix that is greater than zero.
#' Transitions are taken as disallowed for each entry of \code{qmatrix} that is
#' 0.
#' @param data An optional data frame in which the variables represented by
#' \code{subject} and \code{state} can be found.
#' @param censor A state, or vector of states, which indicates censoring.  See
#' \code{\link{msm}}.
#' @param censor.states Specifies the underlying states which censored
#' observations can represent.  See \code{\link{msm}}.
#' @return The estimated transition intensity matrix.  This can be used as the
#' \code{qmatrix} argument to \code{\link{msm}}.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{statetable.msm}}
#' @keywords models
#' @examples
#' 
#' data(cav)
#' twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
#' c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
#' statetable.msm(state, PTNUM, data=cav)
#' crudeinits.msm(state ~ years, PTNUM, data=cav, qmatrix=twoway4.q)
#' 
#' @export crudeinits.msm
crudeinits.msm <- function(formula, subject, qmatrix, data=NULL, censor=NULL, censor.states=NULL)
{
    cens <- msm.form.cmodel(censor, censor.states, qmatrix)
    mf <- model.frame(formula, data=data, na.action=NULL)
    state <- mf[,1]
    if (is.factor(state)) state <- as.numeric(as.character(state))
    time <- mf[,2]
    n <- length(state)
    if (missing(subject)) subject <- rep(1, n)
    else if (!is.null(data))
        subject <- eval(substitute(subject), as.list(data), parent.frame())
    if (is.null(subject)) subject <- rep(1, n)
    notna <- !is.na(subject) & !is.na(time) & !is.na(state)
    subject <- subject[notna]; time <- time[notna]; state <- state[notna]
    msm.check.qmatrix(qmatrix)
    msm.check.state(nrow(qmatrix), state, cens$censor, list(hidden=FALSE))
    msm.check.times(time, subject, state)
    nocens <- (! (state %in% cens$censor) )
    state <- state[nocens]; subject <- subject[nocens]; time <- time[nocens]
    n <- length(state)
    lastsubj <- !duplicated(subject, fromLast=TRUE)
    timecontrib <- ifelse(lastsubj, NA, c(time[2:n], 0) - time)
    tottime <- tapply(timecontrib[!lastsubj], state[!lastsubj], sum) # total time spent in each state
    ntrans <- statetable.msm(state, subject, data=NULL) # table of transitions
    nst <- nrow(qmatrix)
    estmat <- matrix(0, nst, nst)
    rownames(estmat) <- colnames(estmat) <- paste(1:nst)
    tab <- sweep(ntrans, 1, tottime, "/")
    for (i in 1:nst) 
        for (j in 1:nst)
            if ((paste(i) %in% rownames(tab)) && (paste(j) %in% colnames(tab)))
                estmat[paste(i), paste(j)] <- tab[paste(i),paste(j)]#
    ## If no observed transitions but transition is permitted,
    ## set initial value to small number relative to other initial values
    estmat[estmat == 0 & qmatrix==1] <- mean(tab[tab>0]) / 100 # 
    estmat[qmatrix == 0] <- 0 #
    estmat <- msm.fixdiag.qmatrix(estmat)
    rownames(estmat) <- rownames(qmatrix)
    colnames(estmat) <- colnames(qmatrix)
    estmat
}

### Construct a model with time-dependent transition intensities.
### Form a new dataset with censored states and extra covariate, and
### form a new censor model, given change times in tcut

msm.pci <- function(tcut, mf, qmodel, cmodel, covariates)
{
    if (!is.numeric(tcut)) stop("Expected \"tcut\" to be a numeric vector of change points")
    ## Make new dataset with censored observations at time cut points
    ntcut <- length(tcut)
    npts <- length(unique(model.extract(mf, "subject")))
    nextra <- ntcut*npts
    basenames <- c("(state)","(time)","(subject)","(obstype)","(obstrue)","(obs)","(pci.imp)")
    covnames <- setdiff(colnames(mf), basenames)
    extra <- mf[rep(1,nextra),]
    extra$"(state)" = rep(NA, nextra)
    extra$"(time)" = rep(tcut, npts)
    extra$"(subject)" <- rep(unique(model.extract(mf, "subject")), each=ntcut)
    extra$"(obstype)" = rep(1, nextra)
    extra$"(obstrue)" = rep(TRUE, nextra)
    extra$"(obs)" <- NA
    extra$"(pci.imp)" = 1
    extra[, covnames] <- NA
    mf$"(pci.imp)" <- 0

    ## Merge new and old observations
    new <- rbind(mf, extra)
    new <- new[order(new$"(subject)", new$"(time)"),]
    label <- if (cmodel$ncens > 0) max(cmodel$censor)*2 else qmodel$nstates + 1
    new$"(state)"[is.na(new$"(state)")] <- label
    ## Only keep cutpoints within range of each patient's followup
    mintime <- tapply(mf$"(time)", mf$"(subject)", min)
    maxtime <- tapply(mf$"(time)", mf$"(subject)", max)
    ptminmax <- data.frame(subject = names(mintime), mintime, maxtime) 
    new$"(mintime)" <- ptminmax$mintime[match(new$`(subject)`, ptminmax$subject)]
    new$"(maxtime)" <- ptminmax$maxtime[match(new$`(subject)`, ptminmax$subject)]
    new <- new[new$"(time)" >= new$"(mintime)" & new$"(time)" <= new$"(maxtime)", ]
    
    prevsubj <- c(NA,new$"(subject)"[1:(nrow(new)-1)]); nextsubj <- c(new$"(subject)"[2:nrow(new)], NA)
    prevtime <- c(NA,new$"(time)"[1:(nrow(new)-1)]); nexttime <- c(new$"(time)"[2:nrow(new)], NA)
    prevstate <- c(NA,new$"(state)"[1:(nrow(new)-1)]); nextstate <- c(new$"(state)"[2:nrow(new)], NA)

    ## Don't label imputed states as censored if the next observation has obstype=2,
    ## because we know the state at the imputed time is the same as the previous observation 
    nextobstype <- c(new$"(obstype)"[2:nrow(new)], NA)
    ot2 <- new$"(pci.imp)"==1 & nextobstype==2
    new$"(state)"[ot2] <- prevstate[ot2]
    
    ## Drop imputed observations at times when there was already an observation
    ## assumes there wasn't already duplicated obs times
    new <- new[!((new$"(subject)"==prevsubj & new$"(time)"==prevtime & new$"(state)"==label & prevstate!=label) |
                 (new$"(subject)"==nextsubj & new$"(time)"==nexttime & new$"(state)"==label & nextstate!=label))
               ,]

    ## Carry last value forward for other covariates
    if (length(covnames) > 0) {
        eind <- which(is.na(new[,covnames[1]]) & new$"(pci.imp)"==1)
        while(length(eind) > 0){
            new[eind,covnames] <- new[eind - 1, covnames]
            eind <- which(is.na(new[,covnames[1]]))
        }
    }

    ## Check range of cut points
    if (any(tcut <= min(mf$"(time)")))
      warning("Time cut point", if (sum(tcut <= min(mf$"(time)")) > 1) "s " else " ",
                paste(tcut[tcut<=min(mf$"(time)")],collapse=","),
                " less than or equal to minimum observed time of ",min(mf$"(time)"))
    if (any(tcut >= max(mf$"(time)")))
        warning("Time cut point", if (sum(tcut >= max(mf$"(time)")) > 1) "s " else " ",
                paste(tcut[tcut>=max(mf$"(time)")],collapse=","),
                " greater than or equal to maximum observed time of ",max(mf$"(time)"))
    tcut <- tcut[tcut > min(mf$"(time)") & tcut < max(mf$"(time)")]
    ntcut <- length(tcut)

    if (ntcut==0)
        res <- NULL # no cut points in range of data, continue with no time-dependent model
    else {
        ## Insert new covariate in data representing time period
        tcovlabel <- "timeperiod"
        if (any(covnames=="timeperiod")) stop("Cannot have a covariate called \"timeperiod\" if \"pci\" is supplied")
        tcov <- factor(cut(new$"(time)", c(-Inf,tcut,Inf), right=FALSE))
        levs <- levels(tcov)
        levels(tcov) <- gsub(" ","", levs) # get rid of spaces in e.g. [10, Inf) levels
        assign(tcovlabel, tcov)
        new[,tcovlabel] <- tcov

        ## Add "+ timeperiod" to the Q covariates formula.
        current.covs <- if(attr(mf,"ncovs")>0) attr(terms(covariates), "term.labels") else NULL
        covariates <- reformulate(c(current.covs, tcovlabel))
        attr(new, "covnames") <- c(attr(mf, "covnames"), "timeperiod")
        attr(new, "covnames.q") <- c(attr(mf, "covnames.q"), "timeperiod")
        attr(new, "ncovs") <- attr(mf, "ncovs") + 1
        ## New censoring model
        cmodel$ncens <- cmodel$ncens + 1
        cmodel$censor <- c(cmodel$censor, label)
        cmodel$states <- c(cmodel$states, transient.msm(qmatrix=qmodel$imatrix))
        cmodel$index <- if (is.null(cmodel$index)) 1 else cmodel$index
        cmodel$index <- c(cmodel$index, length(cmodel$states) + 1)
        res <- list(mf=new, covariates=covariates, cmodel=cmodel, tcut=tcut)
    }
    res
}

msm.check.covlist <- function(covlist, qemodel) {
    check.numnum <- function(str)
        length(grep("^[0-9]+-[0-9]+$", str)) == length(str)
    num <- sapply(names(covlist), check.numnum)
    if (!all(num)) {
        badnums <- which(!num)
        plural1 <- if (length(badnums)>1) "s" else "";
        plural2 <- if (length(badnums)>1) "e" else "";
        badnames <- paste(paste("\"",names(covlist)[badnums],"\"",sep=""), collapse=",")
        badnums <- paste(badnums, collapse=",")
        stop("Name", plural1, " ", badnames, " of \"covariates\" formula", plural2, " ", badnums, " not in format \"number-number\"")
    }
    for (i in seq_along(covlist))
        if (!inherits(covlist[[i]], "formula"))
            stop("\"covariates\" should be a formula or list of formulae")
    trans <- sapply(strsplit(names(covlist), "-"), as.numeric)
    tm <- if(inherits(qemodel,"msmqmodel")) "transition" else "misclassification"
    qe <- if(inherits(qemodel,"msmqmodel")) "qmatrix" else "ematrix"
    imat <- qemodel$imatrix
    for (i in seq(length.out=ncol(trans))){
        if (imat[trans[1,i],trans[2,i]] != 1)
            stop("covariates on ", names(covlist)[i], " ", tm, " requested, but this is not permitted by the ", qe, ".")
    }
}

## Form indicator matrix for effects that will be fixed to zero when
## "covariates" specified as a list of transition-specific formulae

msm.form.cri <- function(covlist, qmodel, mf, mm, tdmodel) {

    imat <- t(qmodel$imatrix) # order named transitions / misclassifications by row
    tnames <- paste(col(imat)[imat==1],row(imat)[imat==1],sep="-")
    covlabs <- colnames(mm)[-1]
    npars <- qmodel$npars
    cri <- matrix(0, nrow=npars, ncol=length(covlabs), dimnames = list(tnames, covlabs))
    ## time effects specified through "pci" will always be applied to all transitions
    if (!is.null(tdmodel)){
        grep("timeperiod\\[.+,.+\\)", covlabs)
        tdcovs <- grep("timeperiod\\[.+,.+\\)", covlabs)
        cri[,tdcovs] <- 1
    }
    sorti <- function(x) {
        ## converts, e.g. c("b:a:c","d:f","f:e") to c("a:b:c", "d:f", "e:f")
        sapply(lapply(strsplit(x, ":"), sort), paste, collapse=":")
    }
    for (i in 1:npars) {
        if (tnames[i] %in% names(covlist)) {
            covlabsi <- colnames(model.matrix(covlist[[tnames[i]]], data=mf))[-1]
            cri[i, match(sorti(covlabsi), sorti(covlabs))] <- 1
        }
    }
    cri
}

## adapted from stats:::na.omit.data.frame.  ignore handling of
## non-atomic, matrix within df

na.omit.msmdata <- function(object, hidden=FALSE, misc=FALSE, ...) {
    omit <- na.find.msmdata(object, hidden=hidden, misc=misc)
    xx <- object[!omit, , drop = FALSE]
    if (any(omit > 0L)) {
        temp <- setNames(seq(omit)[omit], attr(object, "row.names")[omit])
        attr(temp, "class") <- "omit"
        attr(xx, "na.action") <- temp
    }
    xx
}

na.fail.msmdata <- function(object, hidden=FALSE, misc=FALSE, ...) {
    omit <- na.find.msmdata(object, hidden=hidden, misc=misc)
    if (any(omit))
        stop("Missing values or subjects with only one observation in data")
    else object
}

na.find.msmdata <- function(object, hidden=FALSE, misc=FALSE, ...) {
    subj <- as.character(object[,"(subject)"])
    firstobs <- !duplicated(subj)
    lastobs <- !duplicated(subj, fromLast=TRUE)
    if (misc)
        obstrue <- object[,"(obstrue)"]
    else if (hidden)
        obstrue <- !is.na(object[,"(obstrue)"])
    nm <- names(object)
    omit <- FALSE
    for (j in seq_along(object)) {
        ## Drop all NAs in time, subject as usual
        if (nm[j] %in% c("(time)", "(subject)"))
            omit <- omit | is.na(object[[j]])
        if (nm[j] == "(state)") {
            ## Indicator for missing outcome 
            ## For matrix HMM outcomes ("states"), only drop a row if all columns are NA
            nas <- if (is.matrix(object[[j]])) apply(object[[j]], 1, function(x)all(is.na(x)))
                   else is.na(object[[j]])
            ## Don't drop missing outcomes in HMMs at first obs or if true state known
            ## since there is information then 
            if (hidden) nas[obstrue | firstobs] <- FALSE 
            omit <- omit | nas
        }
        ## Don't drop NAs in obstype at first observation for a subject
        else if (nm[j]=="(obstype)")
            omit <- omit | (is.na(object[[j]]) & !firstobs)
        ## covariates on initial state probs - only drop if NA at initial observation
        else if (j %in% attr(object, "icovi"))
            omit <- omit | (is.na(object[[j]]) & firstobs)
        ## Don't drop NAs in covariates at last observation for a subject
        ## Note NAs in obstrue should have previously been replaced by zeros in msm.form.obstrue, so could assert for this here. 
        else if (nm[j]!="(obstrue)")
            omit <- omit | (is.na(object[[j]]) & !lastobs)
    }
    ## Drop obs with only one subject remaining after NAs have been omitted
    nobspt <- table(subj[!omit])[subj]
    omit <- omit | (nobspt==1)
    omit
}

msm.form.mf.agg <- function(x){
    mf <- x$data$mf; qmodel <- x$qmodel; hmodel <- x$hmodel; cmodel <- x$cmodel
    if (!hmodel$hidden && cmodel$ncens==0){
        mf.trans <- msm.obs.to.fromto(mf)
        msm.check.model(mf.trans$"(fromstate)", mf.trans$"(tostate)", mf.trans$"(obs)", mf.trans$"(subject)", mf.trans$"(obstype)", qmodel$qmatrix, cmodel)
        ## Aggregate over unique from/to/timelag/cov/obstype
        mf.agg <- msm.aggregate.data(mf.trans)
    } else mf.agg <- NULL
    mf.agg
}

msm.obs.to.fromto <- function(mf)
{
    n <- nrow(mf)
    subj <- model.extract(mf, "subject")
    time <- model.extract(mf, "time")
    state <- model.extract(mf, "state")
    firstsubj <- !duplicated(subj)
    lastsubj <- !duplicated(subj, fromLast=TRUE)
    mf.trans <- mf[!lastsubj,,drop=FALSE]   ## retains all covs corresp to the start of the transition
    names(mf.trans)[names(mf.trans)=="(state)"] <- "(fromstate)"
    mf.trans$"(tostate)" <- if(is.matrix(state)) state[!firstsubj,,drop=FALSE] else state[!firstsubj]
    mf.trans$"(obstype)" <- model.extract(mf, "obstype")[!firstsubj] # obstype matched with end of transition
    mf.trans$"(obs)" <- model.extract(mf, "obs")[!firstsubj]
    mf.trans$"(timelag)" <- diff(time)[!firstsubj[-1]]
    ## NOTE not kept constants npts, ncovs, covlabels, covdata, hcovdata, covmeans
    ## NOTE: orig returned dat$time, dat$obstype.obs of orig length, assume not needed
    ## NOTE: firstsubj in old only used in msm.aggregate.hmmdata, which is not used
    mf.trans
}

### Aggregate the data by distinct values of time lag, covariate values, from state, to state, observation type
### Result is passed to the C likelihood function (for non-hidden multi-state models)

msm.aggregate.data <- function(mf.trans)
{
    n <- nrow(mf.trans)
    apaste <- do.call("paste", mf.trans[,c("(fromstate)","(tostate)","(timelag)","(obstype)", attr(mf.trans, "covnames"))])
    ma <- mf.trans[!duplicated(apaste),]
    ma <- ma[order(unique(apaste)),]
    ma$"(nocc)" <- as.numeric(table(apaste))
    apaste2 <- ma[,"(timelag)"]
    if (attr(ma, "ncovs") > 0)
        apaste2 <- paste(apaste2,  do.call("paste", ma[,attr(ma, "covnames"),drop=FALSE]))
    ## which unique timelag/cov combination each row of aggregated data corresponds to
    ## lik.c needs this to know when to recalculate the P matrix.
    ma$"(whicha)" <- match(apaste2, sort(unique(apaste2)))
    ## for Fisher information: number of obs over timelag/covs starting in
    ## fromstate, replicated for all tostates.
    apaste3 <- paste(ma$"(fromstate)", apaste2)
    ma$"(noccsum)" <- unname(tapply(ma$"(nocc)", apaste2, sum)[apaste3])
    ma <- ma[order(apaste2,ma$"(fromstate)",ma$"(tostate)"),]
    ## NOTE not kept covdata, hcovdata, npts, covlabels, covmeans, nobs, ntrans
    ma
}

### Form indicator for which unique timelag/obstype/cov each observation belongs to
### Used in HMMs/censoring to calculate P matrices efficiently in C

msm.form.hmm.agg <- function(mf){
    mf.trans <- msm.obs.to.fromto(mf)
    mf.trans$"(apaste)" <- do.call("paste", mf.trans[,c("(timelag)","(obstype)",attr(mf.trans,"covnames"))])
    mf.trans$"(pcomb)" <- match(mf.trans$"(apaste)", unique(mf.trans$"(apaste)"))
    mf$"(pcomb)" <- NA
    mf$"(pcomb)"[duplicated(model.extract(mf,"subject"))] <- mf.trans$"(pcomb)"
    mf$"(pcomb)"
}

### FORM DESIGN MATRICES FOR COVARIATE MODELS.
msm.form.mm.cov <- function(x){
    mm.cov <- model.matrix.wrap(x$covariates, x$data$mf)
    msm.center.covs(mm.cov, attr(x$data$mf,"covmeans"), x$center)
}

msm.form.mm.cov.agg <- function(x){
    mm.cov.agg <- if (x$hmodel$hidden || (x$cmodel$ncens > 0)) NULL else model.matrix(x$covariates, x$data$mf.agg)
    msm.center.covs(mm.cov.agg, attr(x$data$mf.agg,"covmeans"), x$center)
}

msm.form.mm.mcov <- function(x){
    mm.mcov <- if (x$emodel$misc) model.matrix.wrap(x$misccovariates, x$data$mf) else NULL
    msm.center.covs(mm.mcov, attr(x$data$mf,"covmeans"), x$center)
}

msm.form.mm.hcov <- function(x){
    hcov <- x$hcovariates
    nst <- x$qmodel$nstates
    if (x$hmodel$hidden) {
        mm.hcov <- vector(mode="list", length=nst)
        if (is.null(hcov))
            hcov <- rep(list(~1), nst)
        for (i in seq_len(nst)){
            if (is.null(hcov[[i]])) hcov[[i]] <- ~1
            mm.hcov[[i]] <- model.matrix.wrap(hcov[[i]], x$data$mf)
            mm.hcov[[i]] <- msm.center.covs(mm.hcov[[i]], attr(x$data$mf,"covmeans"), x$center)
        }
    } else mm.hcov <- NULL
    mm.hcov
}

msm.form.mm.icov <- function(x){
    if (x$hmodel$hidden) {
        if (is.null(x$initcovariates)) x$initcovariates <- ~1
        mm.icov <- model.matrix.wrap(x$initcovariates, x$data$mf[!duplicated(x$data$mf$"(subject)"),])
    } else mm.icov <- NULL
    msm.center.covs(mm.icov, attr(x$data$mf,"covmeans"), x$center)
}

model.matrix.wrap <- function(formula, data){
    mm <- model.matrix(formula, data)
    polys <- unlist(attr(mm, "contrasts") == "contr.poly")
    covlist <- paste(names(polys),collapse=",")
    if (any(polys))
        warning(sprintf("Polynomial factor contrasts (found for covariates \"%s\") not supported in msm output functions.  Use treatment contrasts for ordered factors", covlist))   
    mm
}

msm.center.covs <- function(covmat, cm, center=TRUE){
    if (is.null(covmat)||is.null(cm)) return(covmat)
    means <- cm[colnames(covmat)]
    attr(covmat, "means") <- means[-1]
    if (center)
        covmat <- sweep(covmat, 2, means)
    covmat
}

expand.data <- function(x){
    x$data$mf.agg <- msm.form.mf.agg(x)
    x$data$mm.cov <- msm.form.mm.cov(x)
    x$data$mm.cov.agg <- msm.form.mm.cov.agg(x)
    x$data$mm.mcov <- msm.form.mm.mcov(x)
    x$data$mm.hcov <- msm.form.mm.hcov(x)
    x$data$mm.icov <- msm.form.mm.icov(x)
    x$data
}



#' Extract original data from \code{msm} objects.
#' 
#' Extract the data from a multi-state model fitted with \code{msm}.
#' 
#' 
#' @aliases model.frame.msm model.matrix.msm
#' @param formula A fitted multi-state model object, as returned by
#' \code{\link{msm}}.
#' @param agg Return the model frame in the efficient aggregated form used to
#' calculate the likelihood internally for non-hidden Markov models.  This has
#' one row for each unique combination of from-state, to-state, time lag,
#' covariate value and observation type.  The variable named \code{"(nocc)"}
#' counts how many observations of that combination there are in the original
#' data.
#' @param object A fitted multi-state model object, as returned by
#' \code{\link{msm}}.
#' @param model \code{"intens"} to return the design matrix for covariates on
#' intensities, \code{"misc"} for misclassification probabilities, \code{"hmm"}
#' for a general hidden Markov model, and \code{"inits"} for initial state
#' probabilities in hidden Markov models.
#' @param state State corresponding to the required covariate design matrix in
#' a hidden Markov model.
#' @param ... Further arguments (not used).
#' @return \code{model.frame} returns a data frame with all the original
#' variables used for the model fit, with any missing data removed (see
#' \code{na.action} in \code{\link{msm}}).  The state, time, subject,
#' \code{obstype} and \code{obstrue} variables are named \code{"(state)"},
#' \code{"(time)"}, \code{"(subject)"}, \code{"(obstype)"} and
#' \code{"(obstrue)"} respectively (note the brackets).  A variable called
#' \code{"(obs)"} is the observation number from the original data before any
#' missing data were dropped.  The variable \code{"(pcomb)"} is used for
#' computing the likelihood for hidden Markov models, and identifies which
#' distinct time difference, \code{obstype} and covariate values (thus which
#' distinct interval transition probability matrix) each observation
#' corresponds to.
#' 
#' The model frame object has some other useful attributes, including
#' \code{"usernames"} giving the user's original names for these variables
#' (used for model refitting, e.g. in bootstrapping or cross validation) and
#' \code{"covnames"} identifying which ones are covariates.
#' 
#' \code{model.matrix} returns a design matrix for a part of the model that
#' includes covariates.  The required part is indicated by the \code{"model"}
#' argument.
#' 
#' For time-inhomogeneous models fitted with \code{"pci"}, these datasets will
#' have imputed observations at each time change point, indicated where the
#' variable \code{"(pci.imp)"} in the model frame is 1.  The model matrix for
#' intensities will have factor contrasts for the \code{timeperiod} covariate.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}, \code{\link{model.frame}},
#' \code{\link{model.matrix}}.
#' @keywords models
#' @export 
model.frame.msm <- function(formula, agg=FALSE, ...){
    x <- formula
    if (agg) x$data$mf.agg else x$data$mf
}

#' @rdname model.frame.msm 
#' @export 
model.matrix.msm <- function(object, model="intens", state=1, ...){
    switch(model,
           intens=msm.form.mm.cov(object),
           misc=msm.form.mm.mcov(object),
           hmm=msm.form.mm.hcov(object)[[state]],
           init=msm.form.mm.icov(object))
}

.onLoad <- function(libname, pkgname) {
    assign("msm.globals", new.env(), envir=parent.env(environment()))
    assign("nliks", 0, envir=msm.globals)
}
