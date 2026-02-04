# Transition intensity matrix

Extract the estimated transition intensity matrix, and the corresponding
standard errors, from a fitted multi-state model at a given set of
covariate values.

## Usage

``` r
qmatrix.msm(
  x,
  covariates = "mean",
  sojourn = FALSE,
  ci = c("delta", "normal", "bootstrap", "none"),
  cl = 0.95,
  B = 1000,
  cores = NULL
)
```

## Arguments

- x:

  A fitted multi-state model, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covariates:

  The covariate values at which to estimate the intensity matrix. This
  can either be:  

  the string `"mean"`, denoting the means of the covariates in the data
  (this is the default),  

  the number `0`, indicating that all the covariates should be set to
  zero,  

  or a list of values, with optional names. For example

  `list (60, 1)`

  where the order of the list follows the order of the covariates
  originally given in the model formula. Or more clearly, a named list,

  `list (age = 60, sex = 1)`

  If some covariates are specified but not others, the missing ones
  default to zero.

  With `covariates="mean"`, for factor / categorical variables, the mean
  of the 0/1 dummy variable for each factor level is used, representing
  an average over all values in the data, rather than a specific factor
  level.

- sojourn:

  Set to TRUE if the estimated sojourn times and their standard errors
  should also be returned.

- ci:

  If `"delta"` (the default) then confidence intervals are calculated by
  the delta method, or by simple transformation of the Hessian in the
  very simplest cases. Normality on the log scale is assumed.

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
  the distribution of the MLEs.

- cores:

  Number of cores to use for bootstrapping using parallel processing.
  See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details.

## Value

A list with components:

- estimate:

  Estimated transition intensity matrix.

- SE:

  Corresponding approximate standard errors.

- L:

  Lower confidence limits

- U:

  Upper confidence limits

Or if `ci="none"`, then `qmatrix.msm` just returns the estimated
transition intensity matrix.

If `sojourn` is `TRUE`, extra components called `sojourn`, `sojournSE`,
`sojournL` and `sojournU` are included, containing the estimates,
standard errors and confidence limits, respectively, of the mean sojourn
times in each transient state.

The default print method for objects returned by `qmatrix.msm` presents
estimates and confidence limits. To present estimates and standard
errors, do something like

`qmatrix.msm(x)[c("estimates","SE")]`

## Details

Transition intensities and covariate effects are estimated on the log
scale by [`msm`](https://chjackson.github.io/msm/reference/msm.md). A
covariance matrix is estimated from the Hessian of the maximised
log-likelihood.

A more practically meaningful parameterisation of a continuous-time
Markov model with transition intensities \\q\_{rs}\\ is in terms of the
mean sojourn times \\-1 / q\_{rr}\\ in each state \\r\\ and the
probabilities that the next move of the process when in state \\r\\ is
to state \\s\\, \\-q\_{rs} / q\_{rr}\\.

## See also

[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md),
[`sojourn.msm`](https://chjackson.github.io/msm/reference/sojourn.msm.md),
[`deltamethod`](https://chjackson.github.io/msm/reference/deltamethod.md),
[`ematrix.msm`](https://chjackson.github.io/msm/reference/ematrix.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
