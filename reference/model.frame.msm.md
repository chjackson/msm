# Extract original data from `msm` objects.

Extract the data from a multi-state model fitted with `msm`.

## Usage

``` r
# S3 method for class 'msm'
model.frame(formula, agg = FALSE, ...)

# S3 method for class 'msm'
model.matrix(object, model = "intens", state = 1, ...)
```

## Arguments

- formula:

  A fitted multi-state model object, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- agg:

  Return the model frame in the efficient aggregated form used to
  calculate the likelihood internally for non-hidden Markov models. This
  has one row for each unique combination of from-state, to-state, time
  lag, covariate value and observation type. The variable named
  `"(nocc)"` counts how many observations of that combination there are
  in the original data.

- ...:

  Further arguments (not used).

- object:

  A fitted multi-state model object, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- model:

  `"intens"` to return the design matrix for covariates on intensities,
  `"misc"` for misclassification probabilities, `"hmm"` for a general
  hidden Markov model, and `"inits"` for initial state probabilities in
  hidden Markov models.

- state:

  State corresponding to the required covariate design matrix in a
  hidden Markov model.

## Value

`model.frame` returns a data frame with all the original variables used
for the model fit, with any missing data removed (see `na.action` in
[`msm`](https://chjackson.github.io/msm/reference/msm.md)). The state,
time, subject, `obstype` and `obstrue` variables are named `"(state)"`,
`"(time)"`, `"(subject)"`, `"(obstype)"` and `"(obstrue)"` respectively
(note the brackets). A variable called `"(obs)"` is the observation
number from the original data before any missing data were dropped. The
variable `"(pcomb)"` is used for computing the likelihood for hidden
Markov models, and identifies which distinct time difference, `obstype`
and covariate values (thus which distinct interval transition
probability matrix) each observation corresponds to.

The model frame object has some other useful attributes, including
`"usernames"` giving the user's original names for these variables (used
for model refitting, e.g. in bootstrapping or cross validation) and
`"covnames"` identifying which ones are covariates.

`model.matrix` returns a design matrix for a part of the model that
includes covariates. The required part is indicated by the `"model"`
argument.

For time-inhomogeneous models fitted with `"pci"`, these datasets will
have imputed observations at each time change point, indicated where the
variable `"(pci.imp)"` in the model frame is 1. The model matrix for
intensities will have factor contrasts for the `timeperiod` covariate.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),
[`model.frame`](https://rdrr.io/r/stats/model.frame.html),
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
