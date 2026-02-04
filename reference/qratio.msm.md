# Estimated ratio of transition intensities

Compute the estimate and approximate standard error of the ratio of two
estimated transition intensities from a fitted multi-state model at a
given set of covariate values.

## Usage

``` r
qratio.msm(
  x,
  ind1,
  ind2,
  covariates = "mean",
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

- ind1:

  Pair of numbers giving the indices in the intensity matrix of the
  numerator of the ratio, for example, `c(1,2)`.

- ind2:

  Pair of numbers giving the indices in the intensity matrix of the
  denominator of the ratio, for example, `c(2,1)`.

- covariates:

  The covariate values at which to estimate the intensities. This can
  either be:  

  the string `"mean"`, denoting the means of the covariates in the data
  (this is the default),  

  the number `0`, indicating that all the covariates should be set to
  zero,  

  or a list of values, with optional names. For example

  `list (60, 1)`

  where the order of the list follows the order of the covariates
  originally given in the model formula, or a named list,

  `list (age = 60, sex = 1)`

- ci:

  If `"delta"` (the default) then confidence intervals are calculated by
  the delta method.

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

- cores:

  Number of cores to use for bootstrapping using parallel processing.
  See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details.

## Value

A named vector with elements `estimate`, `se`, `L` and `U` containing
the estimate, standard error, lower and upper confidence limits,
respectively, of the ratio of intensities.

## Details

For example, we might want to compute the ratio of the progression rate
and recovery rate for a fitted model `disease.msm` with a health state
(state 1) and a disease state (state 2). In this case, the progression
rate is the (1,2) entry of the intensity matrix, and the recovery rate
is the (2,1) entry. Thus to compute this ratio with covariates set to
their means, we call

`qratio.msm(disease.msm, c(1,2), c(2,1))` .

Standard errors are estimated by the delta method. Confidence limits are
estimated by assuming normality on the log scale.

## See also

[`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
