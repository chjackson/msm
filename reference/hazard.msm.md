# Calculate tables of hazard ratios for covariates on transition intensities

Hazard ratios are computed by exponentiating the estimated covariate
effects on the log-transition intensities. This function is called by
[`summary.msm`](https://chjackson.github.io/msm/reference/summary.msm.md).

## Usage

``` r
hazard.msm(x, hazard.scale = 1, cl = 0.95)
```

## Arguments

- x:

  Output from [`msm`](https://chjackson.github.io/msm/reference/msm.md)
  representing a fitted multi-state model.

- hazard.scale:

  Vector with same elements as number of covariates on transition rates.
  Corresponds to the increase in each covariate used to calculate its
  hazard ratio. Defaults to all 1.

- cl:

  Width of the symmetric confidence interval to present. Defaults to
  0.95.

## Value

A list of tables containing hazard ratio estimates, one table for each
covariate. Each table has three columns, containing the hazard ratio,
and an approximate upper and lower confidence limit respectively
(assuming normality on the log scale), for each Markov chain transition
intensity.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),
[`summary.msm`](https://chjackson.github.io/msm/reference/summary.msm.md),
[`odds.msm`](https://chjackson.github.io/msm/reference/odds.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
