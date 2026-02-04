# Plots of multi-state models

This produces a plot of the expected probability of survival against
time, from each transient state. Survival is defined as not entering an
absorbing state.

## Usage

``` r
# S3 method for class 'msm'
plot(
  x,
  from = NULL,
  to = NULL,
  range = NULL,
  covariates = "mean",
  legend.pos = NULL,
  xlab = "Time",
  ylab = "Fitted survival probability",
  lwd = 1,
  ...
)
```

## Arguments

- x:

  Output from [`msm`](https://chjackson.github.io/msm/reference/msm.md),
  representing a fitted multi-state model object.

- from:

  States from which to consider survival. Defaults to the complete set
  of transient states.

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

- legend.pos:

  Vector of the \\x\\ and \\y\\ position, respectively, of the legend.

- xlab:

  x axis label.

- ylab:

  y axis label.

- lwd:

  Line width. See [`par`](https://rdrr.io/r/graphics/par.html).

- ...:

  Other arguments to be passed to the generic
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) and
  [`lines`](https://rdrr.io/r/graphics/lines.html) functions.

## Details

Note that while this function is only relevant to models with absorbing
states, models in msm can have any transition structure and do not
necessarily have to have an absorbing state.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
