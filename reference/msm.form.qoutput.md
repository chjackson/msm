# Extract msm model parameter estimates in compact format

Extract estimates and confidence intervals for transition intensities
(or misclassification probabilities), and their covariate effects, in a
tidy matrix format with one row per transition. This is used by the
print method
([`print.msm`](https://chjackson.github.io/msm/reference/print.msm.md))
for `msm` objects. Covariate effects are returned as hazard or odds
ratios, not on the log scale.

## Usage

``` r
msm.form.qoutput(x, covariates = "mean", cl = 0.95, digits = 4, ...)

msm.form.eoutput(x, covariates = "mean", cl = 0.95, digits = 4, ...)
```

## Arguments

- x:

  A fitted multi-state model object, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covariates:

  Covariate values defining the "baseline" parameters (see
  [`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md)).

- cl:

  Width of the symmetric confidence interval to present. Defaults to
  0.95.

- digits:

  Minimum number of significant digits for the formatted character
  matrix returned as an attribute. This is passed to
  [`format`](https://rdrr.io/r/base/format.html). Defaults to 4.

- ...:

  Other arguments to be passed to
  [`format`](https://rdrr.io/r/base/format.html).

## Value

A numeric matrix with one row per transition, and one column for each
estimate or confidence limit. The `"formatted"` attribute contains the
same results formatted for pretty printing. `msm.form.qoutput` returns
the transition intensities and their covariates, and `msm.form.eoutput`
returns the misclassification probabilities and their covariates.

## See also

[`print.msm`](https://chjackson.github.io/msm/reference/print.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
