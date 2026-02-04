# Plot of observed and expected prevalences

Provides a rough indication of goodness of fit of a multi-state model,
by estimating the observed numbers of individuals occupying a state at a
series of times, and plotting these against forecasts from the fitted
model, for each state. Observed prevalences are indicated as solid
lines, expected prevalences as dashed lines.

## Usage

``` r
# S3 method for class 'prevalence.msm'
plot(
  x,
  mintime = NULL,
  maxtime = NULL,
  timezero = NULL,
  initstates = NULL,
  interp = c("start", "midpoint"),
  censtime = Inf,
  subset = NULL,
  covariates = "population",
  misccovariates = "mean",
  piecewise.times = NULL,
  piecewise.covariates = NULL,
  xlab = "Times",
  ylab = "Prevalence (%)",
  lwd.obs = 1,
  lwd.exp = 1,
  lty.obs = 1,
  lty.exp = 2,
  col.obs = "blue",
  col.exp = "red",
  legend.pos = NULL,
  ...
)
```

## Arguments

- x:

  A fitted multi-state model produced by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- mintime:

  Minimum time at which to compute the observed and expected prevalences
  of states.

- maxtime:

  Maximum time at which to compute the observed and expected prevalences
  of states.

- timezero:

  Initial time of the Markov process. Expected values are forecasted
  from here. Defaults to the minimum of the observation times given in
  the data.

- initstates:

  Optional vector of the same length as the number of states. Gives the
  numbers of individuals occupying each state at the initial time, to be
  used for forecasting expected prevalences. The default is those
  observed in the data. These should add up to the actual number of
  people in the study at the start.

- interp:

  Interpolation method for observed states, see
  [`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md).

- censtime:

  Subject-specific maximum follow-up times, see
  [`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md).

- subset:

  Vector of the subject identifiers to calculated observed prevalences
  for.

- covariates:

  Covariate values for which to forecast expected state occupancy. See
  [`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md)
  — if this function runs too slowly, as it may if there are continuous
  covariates, replace `covariates="population"` with
  `covariates="mean"`.

- misccovariates:

  (Misclassification models only) Values of covariates on the
  misclassification probability matrix. See
  [`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md).

- piecewise.times:

  Times at which piecewise-constant intensities change. See
  [`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md).

- piecewise.covariates:

  Covariates on which the piecewise-constant intensities depend. See
  [`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md).

- xlab:

  x axis label.

- ylab:

  y axis label.

- lwd.obs:

  Line width for observed prevalences. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- lwd.exp:

  Line width for expected prevalences. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- lty.obs:

  Line type for observed prevalences. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- lty.exp:

  Line type for expected prevalences. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- col.obs:

  Line colour for observed prevalences. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- col.exp:

  Line colour for expected prevalences. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- legend.pos:

  Vector of the \\x\\ and \\y\\ position, respectively, of the legend.
  If `NULL` then a default position will be used.

- ...:

  Further arguments to be passed to the generic
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) function.

## Details

See
[`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md)
for details of the assumptions underlying this method.

Observed prevalences are plotted with a solid line, and expected
prevalences with a dotted line.

## References

Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P. Multi-state
Markov models for analysing incomplete disease history data with
illustrations for HIV disease. *Statistics in Medicine* (1994) 13(3):
805–821.

## See also

[`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md)
