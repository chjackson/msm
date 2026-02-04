# Tidy the output of pmatrix.msm and similar functions

This is the method for the generic \`tidy\` function that is used for
tidying the output of
[`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md),
[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md),
[`ematrix.msm`](https://chjackson.github.io/msm/reference/ematrix.msm.md),
[`pnext.msm`](https://chjackson.github.io/msm/reference/pnext.msm.md) or
[`ppass.msm`](https://chjackson.github.io/msm/reference/ppass.msm.md).
This should be called as `tidy()`, not `tidy.msm.est()` or
`tidy.qmatrix()` or anything else.

## Usage

``` r
# S3 method for class 'msm.est'
tidy(x, ...)
```

## Arguments

- x:

  Output of
  [`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md),
  [`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md),
  [`ematrix.msm`](https://chjackson.github.io/msm/reference/ematrix.msm.md),
  [`pnext.msm`](https://chjackson.github.io/msm/reference/pnext.msm.md)
  or
  [`ppass.msm`](https://chjackson.github.io/msm/reference/ppass.msm.md),
  which all return objects of class `"msm.est"`.

- ...:

  Further arguments (unused).
