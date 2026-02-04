# Tidy the output of totlos.msm and similar functions

Note this should be called as `tidy()` not `tidy.msm.totlos()` or
anything else, as this is a method for the generic `tidy()` function.

## Usage

``` r
# S3 method for class 'msm.estbystate'
tidy(x, ...)
```

## Arguments

- x:

  Output of
  [`totlos.msm`](https://chjackson.github.io/msm/reference/totlos.msm.md),
  [`envisits.msm`](https://chjackson.github.io/msm/reference/totlos.msm.md)
  or
  [`efpt.msm`](https://chjackson.github.io/msm/reference/efpt.msm.md),
  which return objects of class `"msm.estbystate"`.

- ...:

  Further arguments (unused).

## Value

A tibble with one row per state, and columns for the estimate, and
confidence intervals if available.
