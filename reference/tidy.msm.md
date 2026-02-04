# Tidy the parameter estimates from an msm model

Tidy the parameter estimates from an msm model

## Usage

``` r
# S3 method for class 'msm'
tidy(x, ...)
```

## Arguments

- x:

  Object returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md),
  representing a fitted multi-state model.

- ...:

  Other arguments (currently unused).

## Value

A "tibble", with one row for each parameter and the following columns
describing the parameter.

- `parclass`: Class of parameters: `intens` (transition intensities),
  `hr` (hazard ratios representing effects of covariates on
  intensities), and their transformed versions `logintens` (log
  intensities) and `loghr` (log hazard ratios).

  For "misclassification" models fitted with the `ematrix` argument to
  `msm`, other classes of parameters include `misc` (misclassification
  probabilities), `logitmisc` (misclassification log odds), `or_misc`
  and `logor_misc` (effects of covariates on misclassification
  probabilities, as odds ratios or log odds ratios, with the first state
  as the reference category).

  For hidden Markov models fitted with the `hmodel` argument to `msm`,
  the parameter class called `hmm` comprises the parameters of the
  distributions of the outcome conditionally on the hidden state.
  Covariates on the location parameter of these distributions are
  included in class `hmmcov`. If initial state occupancy probabilities
  are estimated, these are included in class `initp` (or `initlogodds`
  for the log odds transforms of these), and any covariates on these
  probabilities are included in class `initpcov`.

- `state`: Starting state of the transition for transition intensities,
  and true state for misclassification probabilities or hidden Markov
  model parameters.

- `tostate`: Ending state of the transition for transition intensities,
  and observed state for misclassification probabilities

- `term`: Name of the covariate for covariate effects, or "baseline" for
  the baseline intensity or analogous parameter value. Note that the
  "baseline" parameters are the parameters with covariates set to their
  mean values in the data (stored in e.g. `x$qcmodel$covmeans`), unless
  `msm` was called with `center=FALSE`.

- `estimate`, `std.error`, `conf.low`, `conf.high`: Parameter estimate,
  standard error, and lower and upper confidence limits.

- `statistic`, `p.value`: For covariate effects, the Z-test statistic
  and p-value for a test of the null hypothesis that the covariate
  effect is zero, based on the estimate and standard error.
