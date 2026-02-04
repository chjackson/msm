# Plot empirical and fitted survival curves

Plot a Kaplan-Meier estimate of the survival probability and compare it
with the fitted survival probability from a `msm` model.

## Usage

``` r
# S3 method for class 'survfit.msm'
plot(
  x,
  from = 1,
  to = NULL,
  range = NULL,
  covariates = "mean",
  interp = c("start", "midpoint"),
  ci = c("none", "normal", "bootstrap"),
  B = 100,
  legend.pos = NULL,
  xlab = "Time",
  ylab = "Survival probability",
  lty = 1,
  lwd = 1,
  col = "red",
  lty.ci = 2,
  lwd.ci = 1,
  col.ci = "red",
  mark.time = TRUE,
  col.surv = "blue",
  lty.surv = 2,
  lwd.surv = 1,
  survdata = FALSE,
  ...
)
```

## Arguments

- x:

  Output from [`msm`](https://chjackson.github.io/msm/reference/msm.md),
  representing a fitted multi-state model object.

- from:

  Non-absorbing state from which to consider survival. Defaults to
  state 1. The fitted probabilities will then be calculated as the
  transition probabilities from this state to `to`. The empirical
  survival curve plots survival from the first observation of `from`
  (where this exists) to the first entry time into `to`.

- to:

  Absorbing state to consider. Defaults to the highest-labelled
  absorbing state.

- range:

  Vector of two elements, giving the range of times to plot for.

- covariates:

  Covariate values for which to evaluate the expected probabilities.
  This can either be:  

  the string `"mean"`, denoting the means of the covariates in the data
  (this is the default),  

  the number `0`, indicating that all the covariates should be set to
  zero,  

  or a list of values, with optional names. For example

  `list (60, 1)`

  where the order of the list follows the order of the covariates
  originally given in the model formula, or a named list,

  `list (age = 60, sex = 1)`

  but note the empirical curve is plotted for the full population. To
  consider subsets for the empirical curve, set `survdata=TRUE` to
  extract the survival data and build a survival plot by hand using
  [`plot.survfit`](https://rdrr.io/pkg/survival/man/plot.survfit.html).

- interp:

  If `interp="start"` (the default) then the entry time into the
  absorbing state is assumed to be the time it is first observed in the
  data.

  If `interp="midpoint"` then the entry time into the absorbing state is
  assumed to be halfway between the time it is first observed and the
  previous observation time. This is generally more reasonable for
  "progressive" models with observations at arbitrary times.

- ci:

  If `"none"` (the default) no confidence intervals are plotted. If
  `"normal"` or `"bootstrap"`, confidence intervals are plotted based on
  the respective method in
  [`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md).
  This is very computationally-intensive, since intervals must be
  computed at a series of times.

- B:

  Number of bootstrap or normal replicates for the confidence interval.
  The default is 100 rather than the usual 1000, since these plots are
  for rough diagnostic purposes.

- legend.pos:

  Vector of the \\x\\ and \\y\\ position, respectively, of the legend.

- xlab:

  x axis label.

- ylab:

  y axis label.

- lty:

  Line type for the fitted curve. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- lwd:

  Line width for the fitted curve. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- col:

  Colour for the fitted curve. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- lty.ci:

  Line type for the fitted curve confidence limits. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- lwd.ci:

  Line width for the fitted curve confidence limits. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- col.ci:

  Colour for the fitted curve confidence limits. See
  [`par`](https://rdrr.io/r/graphics/par.html).

- mark.time:

  Mark the empirical survival curve at each censoring point, see
  [`lines.survfit`](https://rdrr.io/pkg/survival/man/lines.survfit.html).

- col.surv:

  Colour for the empirical survival curve, passed to
  [`lines.survfit`](https://rdrr.io/pkg/survival/man/lines.survfit.html).
  See [`par`](https://rdrr.io/r/graphics/par.html).

- lty.surv:

  Line type for the empirical survival curve, passed to
  [`lines.survfit`](https://rdrr.io/pkg/survival/man/lines.survfit.html).
  See [`par`](https://rdrr.io/r/graphics/par.html).

- lwd.surv:

  Line width for the empirical survival curve, passed to
  [`lines.survfit`](https://rdrr.io/pkg/survival/man/lines.survfit.html).
  See [`par`](https://rdrr.io/r/graphics/par.html).

- survdata:

  Set to `TRUE` to return the survival data frame constructed when
  plotting the empirical curve. This can be used for constructing
  survival plots by hand using
  [`plot.survfit`](https://rdrr.io/pkg/survival/man/plot.survfit.html).

- ...:

  Other arguments to be passed to the
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) function which
  draws the fitted curve, or the
  [`lines.survfit`](https://rdrr.io/pkg/survival/man/lines.survfit.html)
  function which draws the empirical curve.

## Details

If the data represent observations of the process at arbitrary times,
then the first occurrence of the absorbing state in the data will
usually be greater than the actual first transition time to that state.
Therefore the Kaplan-Meier estimate of the survival probability will be
an overestimate.

The method of Turnbull (1976) could be used to give a non-parametric
estimate of the time to an interval-censored event, and compared to the
equivalent estimate from a multi-state model. This is implemented in the
CRAN package interval (Fay and Shaw 2010).

This currently only handles time-homogeneous models.

## References

Turnbull, B. W. (1976) The empirical distribution function with
arbitrarily grouped, censored and truncated data. J. R. Statist. Soc. B
38, 290-295.

Fay, MP and Shaw, PA (2010). Exact and Asymptotic Weighted Logrank Tests
for Interval Censored Data: The interval R package. Journal of
Statistical Software. http://www.jstatsoft.org/v36/ i02/. 36 (2):1-34.

## See also

[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html),
[`plot.survfit`](https://rdrr.io/pkg/survival/man/plot.survfit.html),
[`plot.prevalence.msm`](https://chjackson.github.io/msm/reference/plot.prevalence.msm.md)
