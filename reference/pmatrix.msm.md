# Transition probability matrix

Extract the estimated transition probability matrix from a fitted
continuous-time multi-state model for a given time interval, at a given
set of covariate values.

## Usage

``` r
pmatrix.msm(
  x = NULL,
  t = 1,
  t1 = 0,
  covariates = "mean",
  ci = c("none", "normal", "bootstrap"),
  cl = 0.95,
  B = 1000,
  cores = NULL,
  qmatrix = NULL,
  ...
)
```

## Arguments

- x:

  A fitted multi-state model, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- t:

  The time interval to estimate the transition probabilities for, by
  default one unit.

- t1:

  The starting time of the interval. Used for models `x` with
  piecewise-constant intensities fitted using the `pci` option to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md). The
  probabilities will be computed on the interval \[t1, t1+t\].

- covariates:

  The covariate values at which to estimate the transition
  probabilities. This can either be:  

  the string `"mean"`, denoting the means of the covariates in the data
  (this is the default),  

  the number `0`, indicating that all the covariates should be set to
  zero,  

  or a list of values, with optional names. For example

  `list (60, 1)`

  where the order of the list follows the order of the covariates
  originally given in the model formula, or a named list,

  `list (age = 60, sex = 1)`

  If some covariates are specified but not others, the missing ones
  default to zero.

  For time-inhomogeneous models fitted using the `pci` option to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md),
  "covariates" here include only those specified using the `covariates`
  argument to [`msm`](https://chjackson.github.io/msm/reference/msm.md),
  and exclude the artificial covariates representing the time period.

  For time-inhomogeneous models fitted "by hand" by using a
  time-dependent covariate in the `covariates` argument to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), the
  function
  [`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md)
  should be used to to calculate transition probabilities.

- ci:

  If `"normal"`, then calculate a confidence interval for the transition
  probabilities by simulating `B` random vectors from the asymptotic
  multivariate normal distribution implied by the maximum likelihood
  estimates (and covariance matrix) of the log transition intensities
  and covariate effects, then calculating the resulting transition
  probability matrix for each replicate. See, e.g. Mandel (2013) for a
  discussion of this approach.

  If `"bootstrap"` then calculate a confidence interval by
  non-parametric bootstrap refitting. This is 1-2 orders of magnitude
  slower than the `"normal"` method, but is expected to be more
  accurate. See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details of bootstrapping in msm.

  If `"none"` (the default) then no confidence interval is calculated.

- cl:

  Width of the symmetric confidence interval, relative to 1.

- B:

  Number of bootstrap replicates, or number of normal simulations from
  the distribution of the MLEs

- cores:

  Number of cores to use for bootstrapping using parallel processing.
  See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details.

- qmatrix:

  A transition intensity matrix. Either this or a fitted model `x` must
  be supplied. No confidence intervals are available if `qmatrix` is
  supplied.

- ...:

  Optional arguments to be passed to
  [`MatrixExp`](https://chjackson.github.io/msm/reference/MatrixExp.md)
  to control the method of computing the matrix exponential.

## Value

The matrix of estimated transition probabilities \\P(t)\\ in the given
time. Rows correspond to "from-state" and columns to "to-state".

Or if `ci="normal"` or `ci="bootstrap"`, `pmatrix.msm` returns a list
with components `estimates` and `ci`, where `estimates` is the matrix of
estimated transition probabilities, and `ci` is a list of two matrices
containing the upper and lower confidence limits.

## Details

For a continuous-time homogeneous Markov process with transition
intensity matrix \\Q\\, the probability of occupying state \\s\\ at time
\\u + t\\ conditionally on occupying state \\r\\ at time \\u\\ is given
by the \\(r,s)\\ entry of the matrix \\P(t) = \exp(tQ)\\, where
\\\exp()\\ is the matrix exponential.

For non-homogeneous processes, where covariates and hence the transition
intensity matrix \\Q\\ are piecewise-constant in time, the transition
probability matrix is calculated as a product of matrices over a series
of intervals, as explained in
[`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md).

The
[`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md)
function is only necessary for models fitted using a time-dependent
covariate in the `covariates` argument to
[`msm`](https://chjackson.github.io/msm/reference/msm.md). For
time-inhomogeneous models fitted using "pci", `pmatrix.msm` can be used,
with arguments `t` and `t1`, to calculate transition probabilities over
any time period.

## References

Mandel, M. (2013). "Simulation based confidence intervals for functions
with complicated derivatives." The American Statistician 67(2):76-81

## See also

[`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md),
[`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md),
[`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>.
