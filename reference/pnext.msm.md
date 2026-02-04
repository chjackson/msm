# Probability of each state being next

Compute a matrix of the probability of each state \\s\\ being the next
state of the process after each state \\r\\. Together with the mean
sojourn times in each state
([`sojourn.msm`](https://chjackson.github.io/msm/reference/sojourn.msm.md)),
these fully define a continuous-time Markov model.

## Usage

``` r
pnext.msm(
  x,
  covariates = "mean",
  ci = c("normal", "bootstrap", "delta", "none"),
  cl = 0.95,
  B = 1000,
  cores = NULL
)
```

## Arguments

- x:

  A fitted multi-state model, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covariates:

  The covariate values at which to estimate the intensities. This can
  either be:  

  the string `"mean"`, denoting the means of the covariates in the data
  (this is the default),  

  the number `0`, indicating that all the covariates should be set to
  zero,  

  or a list of values, with optional names. For example

  `list (60, 1)`

  where the order of the list follows the order of the covariates
  originally given in the model formula, or a named list,

  `list (age = 60, sex = 1)`

- ci:

  If `"normal"` (the default) then calculate a confidence interval by
  simulating `B` random vectors from the asymptotic multivariate normal
  distribution implied by the maximum likelihood estimates (and
  covariance matrix) of the log transition intensities and covariate
  effects, then transforming.

  If `"bootstrap"` then calculate a confidence interval by
  non-parametric bootstrap refitting. This is 1-2 orders of magnitude
  slower than the `"normal"` method, but is expected to be more
  accurate. See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details of bootstrapping in msm.

  If `"delta"` then confidence intervals are calculated based on the
  delta method SEs of the log rates, but this is not recommended since
  it may not respect the constraint that probabilities are less than
  one.

- cl:

  Width of the symmetric confidence interval to present. Defaults to
  0.95.

- B:

  Number of bootstrap replicates, or number of normal simulations from
  the distribution of the MLEs.

- cores:

  Number of cores to use for bootstrapping using parallel processing.
  See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details.

## Value

The matrix of probabilities that the next move of a process in state
\\r\\ (rows) is to state \\s\\ (columns).

## Details

For a continuous-time Markov process in state \\r\\, the probability
that the next state is \\s\\ is \\-q\_{rs} / q\_{rr}\\, where
\\q\_{rs}\\ is the transition intensity
([`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md)).

A continuous-time Markov model is fully specified by these probabilities
together with the mean sojourn times \\-1/q\_{rr}\\ in each state \\r\\.
This gives a more intuitively meaningful description of a model than the
intensity matrix.

Remember that msm deals with continuous-time, not discrete-time models,
so these are *not* the same as the probability of observing state \\s\\
at a fixed time in the future. Those probabilities are given by
[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md).

## See also

[`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md),[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md),[`qratio.msm`](https://chjackson.github.io/msm/reference/qratio.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
