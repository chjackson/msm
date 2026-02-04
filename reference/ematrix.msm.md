# Misclassification probability matrix

Extract the estimated misclassification probability matrix, and
corresponding confidence intervals, from a fitted multi-state model at a
given set of covariate values.

## Usage

``` r
ematrix.msm(
  x,
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
  [`msm`](https://chjackson.github.io/msm/reference/msm.md)

- covariates:

  The covariate values for which to estimate the misclassification
  probability matrix. This can either be:  

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
  the delta method, or by simple transformation of the Hessian in the
  very simplest cases.

  If `"normal"`, then calculate a confidence interval by simulating `B`
  random vectors from the asymptotic multivariate normal distribution
  implied by the maximum likelihood estimates (and covariance matrix) of
  the multinomial-logit-transformed misclassification probabilities and
  covariate effects, then transforming back.

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

A list with components:

- estimate:

  Estimated misclassification probability matrix. The rows correspond to
  true states, and columns observed states.

- SE:

  Corresponding approximate standard errors.

- L:

  Lower confidence limits.

- U:

  Upper confidence limits.

Or if `ci="none"`, then `ematrix.msm` just returns the estimated
misclassification probability matrix.

The default print method for objects returned by `ematrix.msm` presents
estimates and confidence limits. To present estimates and standard
errors, do something like

`ematrix.msm(x)[c("estimates","SE")]`

## Details

Misclassification probabilities and covariate effects are estimated on
the multinomial-logit scale by
[`msm`](https://chjackson.github.io/msm/reference/msm.md). A covariance
matrix is estimated from the Hessian of the maximised log-likelihood.
From these, the delta method can be used to obtain standard errors of
the probabilities on the natural scale at arbitrary covariate values.
Confidence intervals are estimated by assuming normality on the
multinomial-logit scale.

## See also

[`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
