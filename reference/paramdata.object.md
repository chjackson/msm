# Developer documentation: internal msm parameters object

An object giving information about the parameters of the multi-state
model. Used internally during maximum likelihood estimation and
arranging results. Returned as the `paramdata` component of a fitted
[`msm`](https://chjackson.github.io/msm/reference/msm.md) model object.

## Value

- inits:

  Vector of initial values for distinct parameters which are being
  estimated. These have been transformed to the real line (e.g. by log),
  and exclude parameters being fixed at their initial values, parameters
  defined to be always fixed (e.g. binomial denominators) and parameters
  constrained to equal previous ones.

- plabs:

  Names of parameters in `allinits`.

- allinits:

  Vector of parameter values before estimation, including those which
  are fixed or constrained to equal other parameters, and transformed to
  the real line.

- hmmpars:

  Indices of `allinits` which represent baseline parameters of hidden
  Markov outcome models (thus excluding covariate effects in HMMs and
  initial state occupancy probabilities).

- fixed:

  `TRUE` if all parameters are fixed, `FALSE` otherwise.

- fixedpars:

  Indices of parameters in `allinits` which are fixed, either by
  definition or as requested by the user in the `fixedpars` argument to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md). Excludes
  parameters fixed by constraining to equal other parameters.

- notfixed:

  Indices of parameters which are not fixed by the definition of
  `fixedpars`.

- optpars:

  Indices of parameters in `allinits` being estimated, thus those
  included in `inits`.

- auxpars:

  Indices of "auxiliary" parameters which are always fixed, for example,
  binomial denominators
  ([`hmmBinom`](https://chjackson.github.io/msm/reference/hmm-dists.md))
  and the `which` parameter in
  [`hmmIdent`](https://chjackson.github.io/msm/reference/hmm-dists.md).

- constr:

  Vector of integers, of length `npars`, indicating which sets of
  parameters are constrained to be equal to each other. If two of these
  integers are equal the corresponding parameters are equal. A negative
  element indicates that parameter is defined to be minus some other
  parameter (this is used for covariate effects on transition
  intensities).

- npars:

  Total number of parameters, equal to `length(allinits)`.

- nfix:

  Number of fixed parameters, equal to `length(fixedpars)`.

- nopt:

  Number of parameters being estimated, equal to `length(inits)` and
  `length(optpars)`.

- ndup:

  Number of parameters defined as duplicates of previous parameters by
  equality constraints (currently unused).

- ranges:

  Matrix of defined ranges for each parameter on the natural scale (e.g.
  0 to infinity for rate parameters).

- opt:

  Object returned by the optimisation routine (such as
  [`optim`](https://rdrr.io/r/stats/optim.html)).

- foundse:

  `TRUE` if standard errors are available after optimisation. If `FALSE`
  the optimisation probably hasn't converged.

- lik:

  Minus twice the log likelihood at the parameter estimates.

- deriv:

  Derivatives of the minus twice log likelihood at the parameter
  estimates, if available.

- information:

  Corresponding expected information matrix at the parameter estimates,
  if available.

- params:

  Vector of parameter values after maximum likelihood estimation,
  corresponding to `allinits`, still on the real-line transformed scale.

- covmat:

  Covariance matrix corresponding to `params`.

- ci:

  Matrix of confidence intervals corresponding to `params`, with nominal
  coverage (default 0.95) defined by the `cl` argument of
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- estimates.t:

  Vector of parameter estimates, as `params` but with parameters on
  their natural scales.

## See also

[`msm.object`](https://chjackson.github.io/msm/reference/msm.object.md)
