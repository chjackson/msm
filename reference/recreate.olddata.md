# Convert data stored in msm object to old format

Converts the `data` element of msm objects to the old format.

## Usage

``` r
recreate.olddata(x)
```

## Arguments

- x:

  Object returned by the
  [`msm`](https://chjackson.github.io/msm/reference/msm.md) function,
  representing a fitted multi-state model.

## Value

A list of vectors and matrices in the undocumented ad-hoc format used
for the `data` component of `msm` objects in msm versions 1.3.1 and
earlier.

## Details

This is just provided for convenience and to illustrate the changes. It
is not guaranteed to be complete, and is liable to be withdrawn. Users
who were relying on the previous undocumented format are advised to
upgrade their code to use the new format, which uses model frames and
model design matrices in the standard format used in version 1.4, based
on [`model.frame`](https://rdrr.io/r/stats/model.frame.html) and
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html).
