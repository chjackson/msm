# Transient and absorbing states

Returns the transient and absorbing states of either a fitted model or a
transition intensity matrix.

## Usage

``` r
transient.msm(x = NULL, qmatrix = NULL)

absorbing.msm(x = NULL, qmatrix = NULL)
```

## Arguments

- x:

  A fitted multi-state model as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- qmatrix:

  A transition intensity matrix. The diagonal is ignored and taken to be
  minus the sum of the rest of the row.

## Value

A vector of the ordinal indices of the transient or absorbing states.

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
