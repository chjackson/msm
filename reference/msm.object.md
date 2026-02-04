# Fitted msm model objects

The [`msm`](https://chjackson.github.io/msm/reference/msm.md) function
returns a list with the following components. These are intended for
developers and confident users. To extract results from fitted model
objects, functions such as
[`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md)
or [`print.msm`](https://chjackson.github.io/msm/reference/print.msm.md)
should be used instead.

## Value

- call:

  The original call to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), as returned
  by [`match.call`](https://rdrr.io/r/base/match.call.html).

- Qmatrices:

  A list of matrices. The first component, labelled `logbaseline`, is a
  matrix containing the estimated transition intensities on the log
  scale with any covariates fixed at their means in the data (or at
  zero, if `center=FALSE`). The component labelled `baseline` is the
  equivalent on the untransformed scale. Each remaining component is a
  matrix giving the linear effects of the labelled covariate on the
  matrix of log intensities. To extract an estimated intensity matrix on
  the natural scale, at an arbitrary combination of covariate values,
  use the function
  [`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md).

- QmatricesSE:

  The standard error matrices corresponding to `Qmatrices`.

- QmatricesL,QmatricesU:

  Corresponding lower and upper symmetric confidence limits, of width
  0.95 unless specified otherwise by the `cl` argument.

- Ematrices:

  A list of matrices. The first component, labelled `logitbaseline`, is
  the estimated misclassification probability matrix (expressed as as
  log odds relative to the probability of the true state) with any
  covariates fixed at their means in the data (or at zero, if
  `center=FALSE`). The component labelled `baseline` is the equivalent
  on the untransformed scale. Each remaining component is a matrix
  giving the linear effects of the labelled covariate on the matrix of
  logit misclassification probabilities. To extract an estimated
  misclassification probability matrix on the natural scale, at an
  arbitrary combination of covariate values, use the function
  [`ematrix.msm`](https://chjackson.github.io/msm/reference/ematrix.msm.md).

- EmatricesSE:

  The standard error matrices corresponding to `Ematrices`.

- EmatricesL,EmatricesU:

  Corresponding lower and upper symmetric confidence limits, of width
  0.95 unless specified otherwise by the `cl` argument.

- minus2loglik:

  Minus twice the maximised log-likelihood.

- deriv:

  Derivatives of the minus twice log-likelihood at its maximum.

- estimates:

  Vector of untransformed maximum likelihood estimates. This includes
  parameters being estimated (as returned from
  [`optim`](https://rdrr.io/r/stats/optim.html)) and those that are
  fixed or constrained by the model specification. Transition
  intensities are on the log scale and misclassification probabilities
  are given as log odds relative to the probability of the true state.

- estimates.t:

  Vector of transformed maximum likelihood estimates with intensities
  and probabilities on their natural scales.

- fixedpars:

  Indices of `estimates` which were fixed during the maximum likelihood
  estimation.

- center:

  Indicator for whether the estimation was performed with covariates
  centered on their means in the data.

- covmat:

  Covariance matrix corresponding to `estimates`.

- ci:

  Matrix of confidence intervals corresponding to `estimates.t`

- opt:

  Return value from the optimisation routine (such as
  [`optim`](https://rdrr.io/r/stats/optim.html) or
  [`nlm`](https://rdrr.io/r/stats/nlm.html)), giving information about
  the results of the optimisation.

- foundse:

  Logical value indicating whether the Hessian was positive-definite at
  the supposed maximum of the likelihood. If not, the covariance matrix
  of the parameters is unavailable. In these cases the optimisation has
  probably not converged to a maximum.

- data:

  A list giving the data used for the model fit, for use in
  post-processing. To extract it, use the methods
  [`model.frame.msm`](https://chjackson.github.io/msm/reference/model.frame.msm.md)
  or
  [`model.matrix.msm`](https://chjackson.github.io/msm/reference/model.frame.msm.md).

  The format of this element changed in version 1.4 of msm, so that it
  now contains a
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html) object `mf`
  with all the variables used in the model. The previous format (an
  ad-hoc list of vectors and matrices) can be obtained with the function
  `recreate.olddata(msmobject)`, where `msmobject` is the object
  returned by `msm`.

- qmodel:

  A list of objects representing the transition matrix structure and
  options for likelihood calculation. See
  [`qmodel.object`](https://chjackson.github.io/msm/reference/qmodel.object.md)
  for documentation of the components.

- emodel:

  A list of objects representing the misclassification model structure,
  for models specified using the `ematrix` argument to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md). See
  [`emodel.object`](https://chjackson.github.io/msm/reference/emodel.object.md).

- qcmodel:

  A list of objects representing the model for covariates on transition
  intensities. See
  [`qcmodel.object`](https://chjackson.github.io/msm/reference/qcmodel.object.md).

- ecmodel:

  A list of objects representing the model for covariates on
  misclassification probabilities. See
  [`ecmodel.object`](https://chjackson.github.io/msm/reference/ecmodel.object.md).

- hmodel:

  A list of objects representing the hidden Markov model structure. See
  [`hmodel.object`](https://chjackson.github.io/msm/reference/hmodel.object.md).

- cmodel:

  A list giving information about censored states. See
  [`cmodel.object`](https://chjackson.github.io/msm/reference/cmodel.object.md).

- pci:

  Cut points for time-varying intensities, as supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), but
  excluding any that are outside the times observed in the data.

- paramdata:

  A list giving information about the parameters of the multi-state
  model. See
  [`paramdata.object`](https://chjackson.github.io/msm/reference/paramdata.object.md).

- cl:

  Confidence interval width, as supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covariates:

  Formula for covariates on intensities, as supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- misccovariates:

  Formula for covariates on misclassification probabilities, as supplied
  to [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- hcovariates:

  Formula for covariates on hidden Markov model outcomes, as supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- initcovariates:

  Formula for covariates on initial state occupancy probabilities in
  hidden Markov models, as supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- sojourn:

  A list as returned by
  [`sojourn.msm`](https://chjackson.github.io/msm/reference/sojourn.msm.md),
  with components:

  `mean` = estimated mean sojourn times in the transient states, with
  covariates fixed at their means (if center=TRUE) or at zero (if
  center=FALSE).

  `se` = corresponding standard errors.
