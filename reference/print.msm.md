# Print a fitted msm model object

Print a fitted msm model object

## Usage

``` r
# S3 method for class 'msm'
print(x, covariates = NULL, digits = 4, ...)

printnew.msm(x, covariates = NULL, digits = 4, ...)
```

## Arguments

- x:

  Output from [`msm`](https://chjackson.github.io/msm/reference/msm.md),
  representing a fitted multi-state model object.

- covariates:

  Covariates for which to print “baseline” transition intensities or
  misclassification probabilities. See
  [`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md)
  for more details.

- digits:

  Minimum number of significant digits, passed to
  [`format`](https://rdrr.io/r/base/format.html). Defaults to 4.

- ...:

  Other arguments to be passed to
  [`format`](https://rdrr.io/r/base/format.html).

## Value

The object returned by `print.msm` is a numeric matrix with one column
for each estimate or confidence limit for intensities and their
covariates, in the same arrangement as printed, but with the underlying
numbers in full precision. The results formatted for printing are stored
in the `"formatted"` attribute of the object, as a character matrix.
These can alternatively be produced by
[`msm.form.qoutput`](https://chjackson.github.io/msm/reference/msm.form.qoutput.md),
which has no printing side-effect.
[`msm.form.eoutput`](https://chjackson.github.io/msm/reference/msm.form.qoutput.md)
produces the same arrangement for misclassification probabilities
instead of intensities.

## Details

This is the new method of formatting msm objects for printing. The old
method was based on printing lists of matrices. That produced a lot of
wasted space for parameters which were zero, and it was difficult to
match corresponding numbers between matrices. The new method presents
all the transition intensities and covariate effects as a single compact
table, and likewise for misclassification matrices.

Also in the old method, covariate effects were presented as log hazard
ratios or log odds ratios. The log scale is more convenient
mathematically, but unnatural to interpret. The new method presents
hazard ratios for covariates on transition intensities and odds ratios
for misclassification probabilities.

`printnew.msm` is an alias for `print.msm`.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),
[`printold.msm`](https://chjackson.github.io/msm/reference/printold.msm.md),
[`msm.form.qoutput`](https://chjackson.github.io/msm/reference/msm.form.qoutput.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
