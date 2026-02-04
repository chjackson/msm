# Convert a hmodel object to HMM constructor function calls

Convert a hmodel object to HMM constructor function calls

## Usage

``` r
hmodel2list(hmodel, hmmdist = TRUE)
```

## Arguments

- hmodel:

  A list of class `hmodel`, as returned in the `hmodel` component of the
  fitted model object from
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- hmmdist:

  `TRUE` or `FALSE` (see "Value" section).

## Value

If `hmmdist=TRUE`, returns a list of objects of class `hmmdist`. These
are the kind of objects returned by HMM constructor functions such as
[`hmmNorm`](https://chjackson.github.io/msm/reference/hmm-dists.md),
[`hmmPois`](https://chjackson.github.io/msm/reference/hmm-dists.md) etc.
Therefore the list can be passed as the `hmodel` argument to
[`msm`](https://chjackson.github.io/msm/reference/msm.md).

If `hmmdist=FALSE`, returns a list comprised of the corresponding input
arguments for the constructor functions, i.e. parameter values of HMM
emission distributions. The list has one element per state. Each of
these elements has one element per parameter (for univariate HMMs), or
one element per outcome distribution, which in turn has one element per
parameter (for multivariate HMMs).

## Author

Will Hulme `https://github.com/wjchulme` and Chris Jackson.
