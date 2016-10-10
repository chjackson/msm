\name{hmodel.object}
\alias{hmodel.object}
\title{Developer documentation: hidden Markov model structure object}
\description{
A list giving information about the models for the outcome data
conditionally on the states of a hidden Markov model.  Used in internal
computations, and returned in a fitted \code{\link{msm}} model object.  
}
\value{
\item{hidden}{\code{TRUE} for hidden Markov models, \code{FALSE} otherwise.}
\item{nstates}{Number of states, the same as \code{qmodel$nstates}.}
\item{fitted}{\code{TRUE} if the parameter values in \code{pars} are the maximum
  likelihood estimates, \code{FALSE} if they are the initial values.}
\item{models}{The outcome distribution for each hidden state.  A vector
  of length \code{nstates} whose \eqn{r}th entry is the index of the
  state \eqn{r} outcome distributions in the vector of supported distributions.
  The vector of
  supported distributions is given in full by \code{msm:::.msm.HMODELS}: the
  first few are 1 for categorical outcome, 2 for identity, 3 for uniform
  and 4 for normal.
}
\item{labels}{String identifying each distribution in \code{models}.}
\item{npars}{Vector of length \code{nstates} giving the number of
  parameters in each outcome distribution, excluding covariate effects.}
\item{nipars}{Number of initial state occupancy probabilities being
    estimated.  This is zero if \code{est.initprobs=FALSE}, otherwise equal to
    the number of states.}
\item{totpars}{Total number of parameters, equal to 
  \code{sum(npars)}. }
\item{pars}{A vector of length \code{totpars}, made from concatenating a
  list of length \code{nstates} whose \eqn{r}th component is
  vector of the parameters for the state \eqn{r} outcome distribution. 
}
\item{plabs}{List with the names of the parameters in \code{pars}.}
\item{parstate}{A vector of length \code{totpars}, whose \eqn{i}th
  element is the state corresponding to the \eqn{i}th parameter.}
\item{firstpar}{A vector of length \code{nstates} giving the index in
  \code{pars} of the first parameter for each state.}
\item{locpars}{Index in \code{pars} of parameters which can
  have covariates on them. }
\item{initprobs}{Initial state occupancy probabilities, as supplied to
  \code{\link{msm}} (initial values before estimation, if \code{est.initprobs=TRUE}.)}
\item{est.initprobs}{Are initial state occupancy probabilities
  estimated (\code{TRUE} or \code{FALSE}), as supplied in the
  \code{est.initprobs} argument of \code{\link{msm}}.}
\item{ncovs}{Number of covariate effects per parameter in \code{pars},
  with, e.g. factor contrasts expanded.}
\item{coveffect}{Vector of covariate effects, of length \code{sum(ncovs)}.}
\item{covlabels}{Labels of these effects.}
\item{coveffstate}{Vector indicating state corresponding to each element of \code{coveffect}.}
\item{ncoveffs}{Number of covariate effects on HMM outcomes, equal to \code{sum(ncovs)}.}
\item{nicovs}{Vector of length \code{nstates-1} giving the number of
  covariate effects on each initial state occupancy probability
  (log relative to the baseline probability).}
\item{icoveffect}{Vector of length \code{sum(nicovs)} giving covariate effects on initial state occupancy probabilities.}
\item{nicoveffs}{Number of covariate effects on initial state occupancy
  probabilities, equal to \code{sum(nicovs)}.}
\item{constr}{Constraints on (baseline) hidden Markov model outcome parameters,
  as supplied in the \code{hconstraint} argument of \code{\link{msm}},
  excluding covariate effects, converted to a vector 
  and mapped to the set 1,2,3,\ldots if necessary.}
\item{covconstr}{Vector of constraints on covariate effects in hidden Markov outcome models,
  as supplied in the \code{hconstraint} argument of \code{\link{msm}},
  excluding baseline parameters, converted to a vector 
  and mapped to the set 1,2,3,\ldots if necessary.}
\item{ranges}{Matrix of range restrictions for HMM parameters, including
  those given to the \code{hranges} argument to \code{\link{msm}}.}
\item{foundse}{\code{TRUE} if standard errors are available for the estimates.}
\item{initpmat}{Matrix of initial state occupancy probabilities with one
  row for each subject (estimated if \code{est.initprobs=TRUE}).}
\item{ci}{Confidence intervals for baseline HMM outcome parameters.}
\item{covci}{Confidence intervals for covariate effects in HMM outcome models.}
}
\seealso{
  \code{\link{msm.object}},\code{\link{qmodel.object}}, \code{\link{emodel.object}}.
}
