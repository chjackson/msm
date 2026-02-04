# Extract model log-likelihood

Extract the log-likelihood and the number of parameters of a model
fitted with [`msm`](https://chjackson.github.io/msm/reference/msm.md).

## Usage

``` r
# S3 method for class 'msm'
logLik(object, by.subject = FALSE, ...)
```

## Arguments

- object:

  A fitted multi-state model object, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- by.subject:

  Return vector of subject-specific log-likelihoods, which should sum to
  the total log-likelihood.

- ...:

  (unused) further arguments passed to or from other methods.

## Value

The log-likelihood of the model represented by 'object' evaluated at the
maximum likelihood estimates.

Akaike's information criterion can also be computed using
[`AIC`](https://rdrr.io/r/stats/AIC.html)`(object)`.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),[`lrtest.msm`](https://chjackson.github.io/msm/reference/lrtest.msm.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
