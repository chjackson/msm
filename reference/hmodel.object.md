# Developer documentation: hidden Markov model structure object

A list giving information about the models for the outcome data
conditionally on the states of a hidden Markov model. Used in internal
computations, and returned in a fitted
[`msm`](https://chjackson.github.io/msm/reference/msm.md) model object.

## Value

- hidden:

  `TRUE` for hidden Markov models, `FALSE` otherwise.

- nstates:

  Number of states, the same as `qmodel$nstates`.

- fitted:

  `TRUE` if the parameter values in `pars` are the maximum likelihood
  estimates, `FALSE` if they are the initial values.

- models:

  The outcome distribution for each hidden state. A vector of length
  `nstates` whose \\r\\th entry is the index of the state \\r\\ outcome
  distributions in the vector of supported distributions. The vector of
  supported distributions is given in full by `msm:::.msm.HMODELS`: the
  first few are 1 for categorical outcome, 2 for identity, 3 for uniform
  and 4 for normal.

- labels:

  String identifying each distribution in `models`.

- npars:

  Vector of length `nstates` giving the number of parameters in each
  outcome distribution, excluding covariate effects.

- nipars:

  Number of initial state occupancy probabilities being estimated. This
  is zero if `est.initprobs=FALSE`, otherwise equal to the number of
  states.

- totpars:

  Total number of parameters, equal to `sum(npars)`.

- pars:

  A vector of length `totpars`, made from concatenating a list of length
  `nstates` whose \\r\\th component is vector of the parameters for the
  state \\r\\ outcome distribution.

- plabs:

  List with the names of the parameters in `pars`.

- parstate:

  A vector of length `totpars`, whose \\i\\th element is the state
  corresponding to the \\i\\th parameter.

- firstpar:

  A vector of length `nstates` giving the index in `pars` of the first
  parameter for each state.

- locpars:

  Index in `pars` of parameters which can have covariates on them.

- initprobs:

  Initial state occupancy probabilities, as supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md) (initial
  values before estimation, if `est.initprobs=TRUE`.)

- est.initprobs:

  Are initial state occupancy probabilities estimated (`TRUE` or
  `FALSE`), as supplied in the `est.initprobs` argument of
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- ncovs:

  Number of covariate effects per parameter in `pars`, with, e.g. factor
  contrasts expanded.

- coveffect:

  Vector of covariate effects, of length `sum(ncovs)`.

- covlabels:

  Labels of these effects.

- coveffstate:

  Vector indicating state corresponding to each element of `coveffect`.

- ncoveffs:

  Number of covariate effects on HMM outcomes, equal to `sum(ncovs)`.

- nicovs:

  Vector of length `nstates-1` giving the number of covariate effects on
  each initial state occupancy probability (log relative to the baseline
  probability).

- icoveffect:

  Vector of length `sum(nicovs)` giving covariate effects on initial
  state occupancy probabilities.

- nicoveffs:

  Number of covariate effects on initial state occupancy probabilities,
  equal to `sum(nicovs)`.

- constr:

  Constraints on (baseline) hidden Markov model outcome parameters, as
  supplied in the `hconstraint` argument of
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), excluding
  covariate effects, converted to a vector and mapped to the set
  1,2,3,... if necessary.

- covconstr:

  Vector of constraints on covariate effects in hidden Markov outcome
  models, as supplied in the `hconstraint` argument of
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), excluding
  baseline parameters, converted to a vector and mapped to the set
  1,2,3,... if necessary.

- ranges:

  Matrix of range restrictions for HMM parameters, including those given
  to the `hranges` argument to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- foundse:

  `TRUE` if standard errors are available for the estimates.

- initpmat:

  Matrix of initial state occupancy probabilities with one row for each
  subject (estimated if `est.initprobs=TRUE`).

- ci:

  Confidence intervals for baseline HMM outcome parameters.

- covci:

  Confidence intervals for covariate effects in HMM outcome models.

## See also

[`msm.object`](https://chjackson.github.io/msm/reference/msm.object.md),[`qmodel.object`](https://chjackson.github.io/msm/reference/qmodel.object.md),
[`emodel.object`](https://chjackson.github.io/msm/reference/emodel.object.md).
