# Developer documentation: misclassification model structure object

A list giving information about the misclassifications assumed in a
multi-state model fitted with the `ematrix` argument of
[`msm`](https://chjackson.github.io/msm/reference/msm.md). Returned in a
fitted [`msm`](https://chjackson.github.io/msm/reference/msm.md) model
object. This information is converted internally to a `hmodel` object
(see
[`hmodel.object`](https://chjackson.github.io/msm/reference/hmodel.object.md))
for use in likelihood computations.

## Value

- nstates:

  Number of states (same as `qmodel$nstates`).

- npars:

  Number of allowed misclassifications, equal to `sum(imatrix)`.

- imatrix:

  Indicator matrix for allowed misclassifications. This has \\(r,s)\\
  entry 1 if misclassification of true state \\r\\ as observed state
  \\s\\ is possible. diagonal entries are arbitrarily set to 0.

- ematrix:

  Matrix of initial values for the misclassification probabilities,
  supplied as the `ematrix` argument of
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- inits:

  Vector of these initial values, reading across rows of `qmatrix` and
  excluding the diagonal and disallowed transitions.

- constr:

  Indicators for equality constraints on baseline misclassification
  probabilities, taken from the `econstraint` argument to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), and mapped
  if necessary to the set (1,2,3,...)

- ndpars:

  Number of distinct misclassification probabilities, after applying
  equality constraints.

- nipars:

  Number of initial state occupancy probabilities being estimated. This
  is zero if `est.initprobs=FALSE`, otherwise equal to the number of
  states.

- initprobs:

  Initial state occupancy probabilities, as supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md) (initial
  values before estimation, if `est.initprobs=TRUE`.)

- est.initprobs:

  Are initial state occupancy probabilities estimated (`TRUE` or
  `FALSE`), as supplied in the `est.initprobs` argument of
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

## See also

[`msm.object`](https://chjackson.github.io/msm/reference/msm.object.md),[`qmodel.object`](https://chjackson.github.io/msm/reference/qmodel.object.md),
[`hmodel.object`](https://chjackson.github.io/msm/reference/hmodel.object.md).
