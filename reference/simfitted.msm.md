# Simulate from a Markov model fitted using msm

Simulate a dataset from a Markov model fitted using
[`msm`](https://chjackson.github.io/msm/reference/msm.md), using the
maximum likelihood estimates as parameters, and the same observation
times as in the original data.

## Usage

``` r
simfitted.msm(x, drop.absorb = TRUE, drop.pci.imp = TRUE)
```

## Arguments

- x:

  A fitted multi-state model object as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- drop.absorb:

  Should repeated observations in an absorbing state be omitted. Use the
  default of `TRUE` to avoid warnings when using the simulated dataset
  for further [`msm`](https://chjackson.github.io/msm/reference/msm.md)
  fits. Or set to `FALSE` if exactly the same number of observations as
  the original data are needed.

- drop.pci.imp:

  In time-inhomogeneous models fitted using the `pci` option to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), censored
  observations are inserted into the data by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md) at the times
  where the intensity changes, but dropped by default when simulating
  from the fitted model using this function. Set this argument to
  `FALSE` to keep these observations and the corresponding indicator
  variable.

## Value

A dataset with variables as described in
[`simmulti.msm`](https://chjackson.github.io/msm/reference/simmulti.msm.md).

## Details

This function is a wrapper around
[`simmulti.msm`](https://chjackson.github.io/msm/reference/simmulti.msm.md),
and only simulates panel-observed data. To generate datasets with the
exact times of transition, use the lower-level
[`sim.msm`](https://chjackson.github.io/msm/reference/sim.msm.md).

Markov models with misclassified states fitted through the `ematrix`
option to [`msm`](https://chjackson.github.io/msm/reference/msm.md) are
supported, but not general hidden Markov models with `hmodel`. For
misclassification models, this function includes misclassification in
the simulated states.

This function is used for parametric bootstrapping to estimate the null
distribution of the test statistic in
[`pearson.msm`](https://chjackson.github.io/msm/reference/pearson.msm.md).

## See also

[`simmulti.msm`](https://chjackson.github.io/msm/reference/simmulti.msm.md),
[`sim.msm`](https://chjackson.github.io/msm/reference/sim.msm.md),
[`pearson.msm`](https://chjackson.github.io/msm/reference/pearson.msm.md),
[`msm`](https://chjackson.github.io/msm/reference/msm.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
