# Tables of observed and expected prevalences

This provides a rough indication of the goodness of fit of a multi-state
model, by estimating the observed numbers of individuals occupying each
state at a series of times, and comparing these with forecasts from the
fitted model.

## Usage

``` r
prevalence.msm(
  x,
  times = NULL,
  timezero = NULL,
  initstates = NULL,
  covariates = "population",
  misccovariates = "mean",
  piecewise.times = NULL,
  piecewise.covariates = NULL,
  ci = c("none", "normal", "bootstrap"),
  cl = 0.95,
  B = 1000,
  cores = NULL,
  interp = c("start", "midpoint"),
  censtime = Inf,
  subset = NULL,
  plot = FALSE,
  ...
)
```

## Arguments

- x:

  A fitted multi-state model produced by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- times:

  Series of times at which to compute the observed and expected
  prevalences of states.

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

- covariates:

  Covariate values for which to forecast expected state occupancy. With
  the default `covariates="population"`, expected prevalences are
  produced by summing model predictions over the covariates observed in
  the original data, for a fair comparison with the observed
  prevalences. This may be slow, particularly with continuous
  covariates.

  Predictions for fixed covariates can be obtained by supplying
  covariate values in the standard way, as in
  [`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md).
  Therefore if `covariates="population"` is too slow, using the mean
  observed values through `covariates="mean"` may give a reasonable
  approximation.

  This argument is ignored if `piecewise.times` is specified. If there
  are a mixture of time-constant and time-dependent covariates, then the
  values for all covariates should be supplied in
  `piecewise.covariates`.

- misccovariates:

  (Misclassification models only) Values of covariates on the
  misclassification probability matrix for converting expected true to
  expected misclassified states. Ignored if `covariates="population"`,
  otherwise defaults to the mean values of the covariates in the data
  set.

- piecewise.times:

  Times at which piecewise-constant intensities change. See
  [`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md)
  for how to specify this. Ignored if `covariates="population"`. This is
  only required for time-inhomogeneous models specified using explicit
  time-dependent covariates, and should not be used for models specified
  using "pci".

- piecewise.covariates:

  Covariates on which the piecewise-constant intensities depend. See
  [`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md)
  for how to specify this. Ignored if `covariates="population"`.

- ci:

  If `"normal"`, then calculate a confidence interval for the expected
  prevalences by simulating `B` random vectors from the asymptotic
  multivariate normal distribution implied by the maximum likelihood
  estimates (and covariance matrix) of the log transition intensities
  and covariate effects, then calculating the expected prevalences for
  each replicate.

  If `"bootstrap"` then calculate a confidence interval by
  non-parametric bootstrap refitting. This is 1-2 orders of magnitude
  slower than the `"normal"` method, but is expected to be more
  accurate. See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details of bootstrapping in msm.

  If `"none"` (the default) then no confidence interval is calculated.

- cl:

  Width of the symmetric confidence interval, relative to 1

- B:

  Number of bootstrap replicates

- cores:

  Number of cores to use for bootstrapping using parallel processing.
  See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details.

- interp:

  Suppose an individual was observed in states \\S\_{r-1}\\ and \\S_r\\
  at two consecutive times \\t\_{r-1}\\ and \\t_r\\, and we want to
  estimate 'observed' prevalences at a time \\t\\ between \\t\_{r-1}\\
  and \\t_r\\.

  If `interp="start"`, then individuals are assumed to be in state
  \\S\_{r-1}\\ at time \\t\\, the same state as they were at
  \\t\_{r-1}\\.

  If `interp="midpoint"` then if \\t \<= (t\_{r-1} + t_r) / 2\\, the
  midpoint of \\t\_{r-1}\\ and \\t_r\\, the state at \\t\\ is assumed to
  be \\S\_{r-1}\\, otherwise \\S\_{r}\\. This is generally more
  reasonable for "progressive" models.

- censtime:

  Adjustment to the observed prevalences to account for limited
  follow-up in the data.

  If the time is greater than `censtime` and the patient has reached an
  absorbing state, then that subject will be removed from the risk set.
  For example, if patients have died but would only have been observed
  up to this time, then this avoids overestimating the proportion of
  people who are dead at later times.

  This can be supplied as a single value, or as a vector with one
  element per subject (after any `subset` has been taken), in the same
  order as the original data. This vector also only includes subjects
  with complete data, thus it excludes for example subjects with only
  one observation (thus no observed transitions), and subjects for whom
  every observation has missing values. (Note, to help construct this,
  the complete data used for the model fit can be accessed with
  `model.frame(x)`, where `x` is the fitted model object)

  This is ignored if it is less than the subject's maximum observation
  time.

- subset:

  Subset of subjects to calculate observed prevalences for.

- plot:

  Generate a plot of observed against expected prevalences. See
  [`plot.prevalence.msm`](https://chjackson.github.io/msm/reference/plot.prevalence.msm.md)

- ...:

  Further arguments to pass to
  [`plot.prevalence.msm`](https://chjackson.github.io/msm/reference/plot.prevalence.msm.md).

## Value

A list of matrices, with components:

- Observed:

  Table of observed numbers of individuals in each state at each time

- Observed percentages:

  Corresponding percentage of the individuals at risk at each time.

- Expected:

  Table of corresponding expected numbers.

- Expected percentages:

  Corresponding percentage of the individuals at risk at each time.

Or if `ci.boot = TRUE`, the component `Expected` is a list with
components `estimates` and `ci`.  
`estimates` is a matrix of the expected prevalences, and `ci` is a list
of two matrices, containing the confidence limits. The component
`Expected percentages` has a similar format.

## Details

The fitted transition probability matrix is used to forecast expected
prevalences from the state occupancy at the initial time. To produce the
expected number in state \\j\\ at time \\t\\ after the start, the number
of individuals under observation at time \\t\\ (including those who have
died, but not those lost to follow-up) is multiplied by the product of
the proportion of individuals in each state at the initial time and the
transition probability matrix in the time interval \\t\\. The proportion
of individuals in each state at the "initial" time is estimated, if
necessary, in the same way as the observed prevalences.

For misclassification models (fitted using an `ematrix`), this aims to
assess the fit of the full model for the *observed* states. That is, the
combined Markov progression model for the true states and the
misclassification model. Thus, expected prevalences of *true* states are
estimated from the assumed proportion occupying each state at the
initial time using the fitted transition probabiliy matrix. The vector
of expected prevalences of true states is then multiplied by the fitted
misclassification probability matrix to obtain the expected prevalences
of *observed* states.

For general hidden Markov models, the observed state is taken to be the
predicted underlying state from the Viterbi algorithm
([`viterbi.msm`](https://chjackson.github.io/msm/reference/viterbi.msm.md)).
The goodness of fit of these states to the underlying Markov model is
tested.

In any model, if there are censored states, then these are replaced by
imputed values of highest probability from the Viterbi algorithm in
order to calculate the observed state prevalences.

For an example of this approach, see Gentleman *et al.* (1994).

## References

Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P. Multi-state
Markov models for analysing incomplete disease history data with
illustrations for HIV disease. *Statistics in Medicine* (1994) 13(3):
805–821.

Titman, A.C., Sharples, L. D. Model diagnostics for multi-state models.
*Statistical Methods in Medical Research* (2010) 19(6):621-651.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),
[`summary.msm`](https://chjackson.github.io/msm/reference/summary.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
