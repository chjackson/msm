# Likelihood ratio test

Likelihood ratio test between two or more fitted multi-state models

## Usage

``` r
lrtest.msm(...)
```

## Arguments

- ...:

  Two or more fitted multi-state models, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), ordered by
  increasing numbers of parameters.

## Value

A matrix with three columns, giving the likelihood ratio statistic,
difference in degrees of freedom and the chi-squared p-value for a
comparison of the first model supplied with each subsequent model.

## Warning

The comparison between models will only be valid if they are fitted to
the same dataset. This may be a problem if there are missing values and
R's default of 'na.action = na.omit' is used.

The likelihood ratio statistic only has the indicated chi-squared
distribution if the models are nested. An alternative for comparing
non-nested models is Akaike's information criterion. This can be
computed for one or more fitted `msm` models `x,y,...` using
[`AIC`](https://rdrr.io/r/stats/AIC.html)`(x,y,...)`.

## See also

[`logLik.msm`](https://chjackson.github.io/msm/reference/logLik.msm.md),[`msm`](https://chjackson.github.io/msm/reference/msm.md)
