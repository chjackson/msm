# Tidy the output of prevalence.msm

Note this should be called as `tidy()` not `tidy.msm.prevalence()` or
anything else, as this is a method for the generic `tidy()` function.

## Usage

``` r
# S3 method for class 'msm.prevalence'
tidy(x, ...)
```

## Arguments

- x:

  Output of
  [`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md).

- ...:

  Further arguments (unused).

## Value

A tibble with one row per combination of output type (count or
percentage) and state, and columns for observed value, expected value
and confidence limits for the expected value (if available).
