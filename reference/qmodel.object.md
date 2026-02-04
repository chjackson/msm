# Developer documentation: transition model structure object

A list giving information about the structure of states and allowed
transitions in a multi-state model, and options for likelihood
calculation. Used in internal computations, and returned in a fitted
[`msm`](https://chjackson.github.io/msm/reference/msm.md) model object.

## Value

- nstates:

  Number of states

- iso:

  Label for which basic structure the model is isomorphic to in the list
  of structures for which analytic formulae for the transition
  probabilities are implemented in the source file `src/analyticp.c`.
  This list is given by the internal object `msm:::.msm.graphs` which is
  defined and documented in the source file `R/constants.R`.

  `iso` is 0 if the analytic P matrix is not implemented for this
  structure, or if analytic P matrix calculations are disabled using
  `use.analyticp=FALSE` in the call to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- perm:

  Permutation required to convert the base isomorphism into the
  structure of this model. A vector of integers whose \\r\\th element is
  the state number in the base structure representing state \\r\\ in the
  current structure.

- qperm:

  Inverse permutation: vector whose \\r\\th element is the state number
  in the current structure representing the \\r\\th state in the base
  structure.

- npars:

  Number of allowed instantaneous transitions, equal to `sum(imatrix)`.

- imatrix:

  Indicator matrix for allowed instantaneous transitions. This has
  \\(r,s)\\ entry 1 if the transition from \\r\\ to \\s\\ is permitted
  in continuous time, and 0 otherwise. The diagonal entries are
  arbitrarily set to 0.

- qmatrix:

  Matrix of initial values for the transition intensities, supplied as
  the `qmatrix` argument of
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- inits:

  Vector of these initial values, reading across rows of `qmatrix` and
  excluding the diagonal and disallowed transitions.

- constr:

  Indicators for equality constraints on baseline intensities, taken
  from the `qconstraint` argument to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), and mapped
  if necessary to the set (1,2,3,...).

- ndpars:

  Number of distinct allowed instantaneous transitions, after applying
  equality constraints.

- expm:

  Use expm package to calculate matrix exponentials for likelihoods, as
  supplied to the `use.expm` argument of
  [`msm`](https://chjackson.github.io/msm/reference/msm.md). `TRUE` or
  `FALSE`.

## See also

[`msm.object`](https://chjackson.github.io/msm/reference/msm.object.md),[`emodel.object`](https://chjackson.github.io/msm/reference/emodel.object.md),
[`hmodel.object`](https://chjackson.github.io/msm/reference/hmodel.object.md).
