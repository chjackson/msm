# Mean sojourn times from a multi-state model

Estimate the mean sojourn times in the transient states of a multi-state
model and their confidence limits.

## Usage

``` r
sojourn.msm(
  x,
  covariates = "mean",
  ci = c("delta", "normal", "bootstrap", "none"),
  cl = 0.95,
  B = 1000
)
```

## Arguments

- x:

  A fitted multi-state model, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covariates:

  The covariate values at which to estimate the mean sojourn times. This
  can either be:  

  the string `"mean"`, denoting the means of the covariates in the data
  (this is the default),  

  the number `0`, indicating that all the covariates should be set to
  zero,  

  a list of values, with optional names. For example,

  `list(60, 1)`, where the order of the list follows the order of the
  covariates originally given in the model formula, or a named list,
  e.g.

  `list (age = 60, sex = 1)`

- ci:

  If `"delta"` (the default) then confidence intervals are calculated by
  the delta method, or by simple transformation of the Hessian in the
  very simplest cases.

  If `"normal"`, then calculate a confidence interval by simulating `B`
  random vectors from the asymptotic multivariate normal distribution
  implied by the maximum likelihood estimates (and covariance matrix) of
  the log transition intensities and covariate effects, then
  transforming.

  If `"bootstrap"` then calculate a confidence interval by
  non-parametric bootstrap refitting. This is 1-2 orders of magnitude
  slower than the `"normal"` method, but is expected to be more
  accurate. See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details of bootstrapping in msm.

- cl:

  Width of the symmetric confidence interval to present. Defaults to
  0.95.

- B:

  Number of bootstrap replicates, or number of normal simulations from
  the distribution of the MLEs

## Value

A data frame with components:

- estimates:

  Estimated mean sojourn times in the transient states.

- SE:

  Corresponding standard errors.

- L:

  Lower confidence limits.

- U:

  Upper confidence limits.

## Details

The mean sojourn time in a transient state \\r\\ is estimated by \\- 1 /
q\_{rr}\\, where \\q\_{rr}\\ is the \\r\\th entry on the diagonal of the
estimated transition intensity matrix.

A continuous-time Markov model is fully specified by the mean sojourn
times and the probability that each state is next
([`pnext.msm`](https://chjackson.github.io/msm/reference/pnext.msm.md)).
This is a more intuitively meaningful description of a model than the
transition intensity matrix
([`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md)).

Time dependent covariates, or time-inhomogeneous models, are not
supported. This would require the mean of a piecewise exponential
distribution, and the package author is not aware of any general
analytic form for that.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),
[`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md),
[`deltamethod`](https://chjackson.github.io/msm/reference/deltamethod.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
