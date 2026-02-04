# Extract model coefficients

Extract the estimated log transition intensities and the corresponding
linear effects of each covariate.

## Usage

``` r
# S3 method for class 'msm'
coef(object, ...)
```

## Arguments

- object:

  A fitted multi-state model object, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- ...:

  (unused) further arguments passed to or from other methods.

## Value

If there is no misclassification, `coef.msm` returns a list of matrices.
The first component, labelled `logbaseline`, is a matrix containing the
estimated transition intensities on the log scale with any covariates
fixed at their means in the data. Each remaining component is a matrix
giving the linear effects of the labelled covariate on the matrix of log
intensities.  

For misclassification models, `coef.msm` returns a list of lists. The
first component, `Qmatrices`, is a list of matrices as described in the
previous paragraph. The additional component `Ematrices` is a list of
similar format containing the logit-misclassification probabilities and
any estimated covariate effects.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
