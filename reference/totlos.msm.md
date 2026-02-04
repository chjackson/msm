# Total length of stay, or expected number of visits

Estimate the expected total length of stay, or the expected number of
visits, in each state, for an individual in a given period of evolution
of a multi-state model.

## Usage

``` r
totlos.msm(
  x,
  start = 1,
  end = NULL,
  fromt = 0,
  tot = Inf,
  covariates = "mean",
  piecewise.times = NULL,
  piecewise.covariates = NULL,
  num.integ = FALSE,
  discount = 0,
  env = FALSE,
  ci = c("none", "normal", "bootstrap"),
  cl = 0.95,
  B = 1000,
  cores = NULL,
  ...
)

envisits.msm(
  x = NULL,
  start = 1,
  end = NULL,
  fromt = 0,
  tot = Inf,
  covariates = "mean",
  piecewise.times = NULL,
  piecewise.covariates = NULL,
  num.integ = FALSE,
  discount = 0,
  ci = c("none", "normal", "bootstrap"),
  cl = 0.95,
  B = 1000,
  cores = NULL,
  ...
)
```

## Arguments

- x:

  A fitted multi-state model, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- start:

  Either a single number giving the state at the beginning of the
  period, or a vector of probabilities of being in each state at this
  time.

- end:

  States to estimate the total length of stay (or number of visits) in.
  Defaults to all states. This is deprecated, since with the analytic
  solution (see "Details") it doesn't save any computation to only
  estimate for a subset of states.

- fromt:

  Time from which to estimate. Defaults to 0, the beginning of the
  process.

- tot:

  Time up to which the estimate is made. Defaults to infinity, giving
  the expected time spent in or number of visits to the state until
  absorption. However, the calculation will be much more efficient if a
  finite (potentially large) time is specified: see the "Details"
  section. For models without an absorbing state, `t` must be specified.

- covariates:

  The covariate values to estimate for. This can either be:  

  the string `"mean"`, denoting the means of the covariates in the data
  (this is the default),  

  the number `0`, indicating that all the covariates should be set to
  zero,  

  or a list of values, with optional names. For example

  `list (60, 1)`

  where the order of the list follows the order of the covariates
  originally given in the model formula, or a named list,

  `list (age = 60, sex = 1)`

- piecewise.times:

  Times at which piecewise-constant intensities change. See
  [`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md)
  for how to specify this. This is only required for time-inhomogeneous
  models specified using explicit time-dependent covariates, and should
  not be used for models specified using "pci".

- piecewise.covariates:

  Covariates on which the piecewise-constant intensities depend. See
  [`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md)
  for how to specify this.

- num.integ:

  Use numerical integration instead of analytic solution (see below).

- discount:

  Discount rate in continuous time.

- env:

  Supplied to `totlos.msm`. If `TRUE`, return the expected number of
  visits to each state. If `FALSE`, return the total length of stay in
  each state. `envisits.msm` simply calls `totlos.msm` with `env=TRUE`.

- ci:

  If `"normal"`, then calculate a confidence interval by simulating `B`
  random vectors from the asymptotic multivariate normal distribution
  implied by the maximum likelihood estimates (and covariance matrix) of
  the log transition intensities and covariate effects, then calculating
  the total length of stay for each replicate.

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

- ...:

  Further arguments to be passed to the
  [`integrate`](https://rdrr.io/r/stats/integrate.html) function to
  control the numerical integration.

## Value

A vector of expected total lengths of stay (`totlos.msm`), or expected
number of visits (`envisits.msm`), for each transient state.

## Details

The expected total length of stay in state \\j\\ between times \\t_1\\
and \\t_2\\, from the point of view of an individual in state \\i\\ at
time 0, is defined by the integral from \\t_1\\ to \\t_2\\ of the
\\i,j\\ entry of the transition probability matrix \\P(t) = Exp(tQ)\\,
where \\Q\\ is the transition intensity matrix.

The corresponding expected number of visits to state \\j\\ (excluding
the stay in the current state at time 0) is \\\sum\_{i!=j} T_i
Q\_{i,j}\\, where \\T_i\\ is the expected amount of time spent in state
\\i\\.

More generally, suppose that \\\pi_0\\ is the vector of probabilities of
being in each state at time 0, supplied in `start`, and we want the
vector \\\mathbf{x}\\ giving the expected lengths of stay in each state.
The corresponding integral has the following solution (van Loan 1978;
van Rosmalen et al. 2013)

\$\$\mathbf{x} = \left\[ \begin{array}{ll} 1 & \mathbf{0}\_K \end{array}
\right\] Exp(t Q') \left\[ \begin{array}{l} \mathbf{0}\_K\\I_K
\end{array} \right\] \$\$

where \$\$Q' = \left\[ \begin{array}{ll} 0 & \mathbf{\pi}\_0\\
\mathbf{0}\_K & Q - rI_K \end{array} \right\] \$\$

\\\pi_0\\ is the row vector of initial state probabilities supplied in
`start`, \\\mathbf{0}\_K\\ is the row vector of K zeros, \\r\\ is the
discount rate, \\I_K\\ is the K x K identity matrix, and \\Exp\\ is the
matrix exponential.

Alternatively, the integrals can be calculated numerically, using the
[`integrate`](https://rdrr.io/r/stats/integrate.html) function. This may
take a long time for models with many states where \\P(t)\\ is expensive
to calculate. This is required where `tot = Inf`, since the package
author is not aware of any analytic expression for the limit of the
above formula as \\t\\ goes to infinity.

With the argument `num.integ=TRUE`, numerical integration is used even
where the analytic solution is available. This facility is just provided
for checking results against versions 1.2.4 and earlier, and will be
removed eventually. Please let the package maintainer know if any
results are different.

For a model where the individual has only one place to go from each
state, and each state is visited only once, for example a progressive
disease model with no recovery or death, these are equal to the mean
sojourn time in each state. However, consider a three-state
health-disease-death model with transitions from health to disease,
health to death, and disease to death, where everybody starts healthy.
In this case the mean sojourn time in the disease state will be greater
than the expected length of stay in the disease state. This is because
the mean sojourn time in a state is conditional on entering the state,
whereas the expected total time diseased is a forecast for a healthy
individual, who may die before getting the disease.

In the above formulae, \\Q\\ is assumed to be constant over time, but
the results generalise easily to piecewise-constant intensities. This
function automatically handles models fitted using the `pci` option to
[`msm`](https://chjackson.github.io/msm/reference/msm.md). For any other
inhomogeneous models, the user must specify `piecewise.times` and
`piecewise.covariates` arguments to `totlos.msm`.

## References

C. van Loan (1978). Computing integrals involving the matrix
exponential. IEEE Transactions on Automatic Control 23(3)395-404.

J. van Rosmalen, M. Toy and J.F. O'Mahony (2013). A mathematical
approach for evaluating Markov models in continuous time without
discrete-event simulation. Medical Decision Making 33:767-779.

## See also

[`sojourn.msm`](https://chjackson.github.io/msm/reference/sojourn.msm.md),
[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md),
[`integrate`](https://rdrr.io/r/stats/integrate.html),
[`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
