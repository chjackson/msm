# Kaplan Meier estimates of incidence

Compute and plot Kaplan-Meier estimates of the probability that each
successive state has not occurred yet.

## Usage

``` r
plotprog.msm(
  formula,
  subject,
  data,
  legend.pos = NULL,
  xlab = "Time",
  ylab = "1 - incidence probability",
  lwd = 1,
  xlim = NULL,
  mark.time = TRUE,
  ...
)
```

## Arguments

- formula:

  A formula giving the vectors containing the observed states and the
  corresponding observation times. For example,

  `state ~ time`

  Observed states should be in the set `1, ...{}, n`, where `n` is the
  number of states.

- subject:

  Vector of subject identification numbers for the data specified by
  `formula`. If missing, then all observations are assumed to be on the
  same subject. These must be sorted so that all observations on the
  same subject are adjacent.

- data:

  An optional data frame in which the variables represented by `state`,
  `time` and `subject` can be found.

- legend.pos:

  Vector of the \\x\\ and \\y\\ position, respectively, of the legend.

- xlab:

  x axis label.

- ylab:

  y axis label.

- lwd:

  Line width. See [`par`](https://rdrr.io/r/graphics/par.html).

- xlim:

  x axis limits, e.g. c(0,10) for an axis ranging from 0 to 10. Default
  is the range of observation times.

- mark.time:

  Mark the empirical survival curve at each censoring point, see
  [`lines.survfit`](https://rdrr.io/pkg/survival/man/lines.survfit.html).

- ...:

  Other arguments to be passed to the
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) and
  [`lines.survfit`](https://rdrr.io/pkg/survival/man/lines.survfit.html)
  functions.

## Details

If the data represent observations of the process at arbitrary times,
then the first occurrence of the state in the data will usually be
greater than the actual first transition time to that state. Therefore
the probabilities plotted by `plotprog.msm` will be overestimates.

## See also

[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html),
[`plot.survfit`](https://rdrr.io/pkg/survival/man/plot.survfit.html)
