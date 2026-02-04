# Summarise a fitted multi-state model

Summary method for fitted
[`msm`](https://chjackson.github.io/msm/reference/msm.md) models. This
is simply a wrapper around
[`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md)
which produces a table of observed and expected state prevalences for
each time, and for models with covariates,
[`hazard.msm`](https://chjackson.github.io/msm/reference/hazard.msm.md)
to print hazard ratios with 95% confidence intervals for covariate
effects.

## Usage

``` r
# S3 method for class 'msm'
summary(object, hazard.scale = 1, ...)
```

## Arguments

- object:

  A fitted multi-state model object, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- hazard.scale:

  Vector with same elements as number of covariates on transition rates.
  Corresponds to the increase in each covariate used to calculate its
  hazard ratio. Defaults to all 1.

- ...:

  Further arguments passed to
  [`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md).

## Value

A list of class `summary.msm`, with components:

- prevalences:

  Output from
  [`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md).

- hazard:

  Output from
  [`hazard.msm`](https://chjackson.github.io/msm/reference/hazard.msm.md).

- hazard.scale:

  Value of the `hazard.scale` argument.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),[`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md),
[`hazard.msm`](https://chjackson.github.io/msm/reference/hazard.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
