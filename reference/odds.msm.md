# Calculate tables of odds ratios for covariates on misclassification probabilities

Odds ratios are computed by exponentiating the estimated covariate
effects on the logit-misclassification probabilities.

## Usage

``` r
odds.msm(x, odds.scale = 1, cl = 0.95)
```

## Arguments

- x:

  Output from [`msm`](https://chjackson.github.io/msm/reference/msm.md)
  representing a fitted multi-state model.

- odds.scale:

  Vector with same elements as number of covariates on misclassification
  probabilities. Corresponds to the increase in each covariate used to
  calculate its odds ratio. Defaults to all 1.

- cl:

  Width of the symmetric confidence interval to present. Defaults to
  0.95.

## Value

A list of tables containing odds ratio estimates, one table for each
covariate. Each table has three columns, containing the odds ratio, and
an approximate upper 95% and lower 95% confidence limit respectively
(assuming normality on the log scale), for each misclassification
probability.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),
[`hazard.msm`](https://chjackson.github.io/msm/reference/hazard.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
