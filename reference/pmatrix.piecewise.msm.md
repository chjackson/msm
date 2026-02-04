# Transition probability matrix for processes with piecewise-constant intensities

Extract the estimated transition probability matrix from a fitted
non-time-homogeneous multi-state model for a given time interval. This
is a generalisation of
[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md)
to models with time-dependent covariates. Note that
[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md)
is sufficient to calculate transition probabilities for
time-inhomogeneous models fitted using the `pci` argument to
[`msm`](https://chjackson.github.io/msm/reference/msm.md).

## Usage

``` r
pmatrix.piecewise.msm(
  x = NULL,
  t1,
  t2,
  times,
  covariates = NULL,
  ci = c("none", "normal", "bootstrap"),
  cl = 0.95,
  B = 1000,
  cores = NULL,
  qlist = NULL,
  ...
)
```

## Arguments

- x:

  A fitted multi-state model, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md). This should
  be a non-homogeneous model, whose transition intensity matrix depends
  on a time-dependent covariate.

- t1:

  The start of the time interval to estimate the transition
  probabilities for.

- t2:

  The end of the time interval to estimate the transition probabilities
  for.

- times:

  Cut points at which the transition intensity matrix changes.

- covariates:

  A list with number of components one greater than the length of
  `times`. Each component of the list is specified in the same way as
  the `covariates` argument to
  [`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md).
  The components correspond to the covariate values in the intervals

  `(t1, times[1]], (times[1], times[2]], ..., (times[length(times)], t2]`

  (assuming that all elements of `times` are in the interval
  `(t1, t2)`).

- ci:

  If `"normal"`, then calculate a confidence interval for the transition
  probabilities by simulating `B` random vectors from the asymptotic
  multivariate normal distribution implied by the maximum likelihood
  estimates (and covariance matrix) of the log transition intensities
  and covariate effects, then calculating the resulting transition
  probability matrix for each replicate.

  If `"bootstrap"` then calculate a confidence interval by
  non-parametric bootstrap refitting. This is 1-2 orders of magnitude
  slower than the `"normal"` method, but is expected to be more
  accurate. See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details of bootstrapping in msm.

  If `"none"` (the default) then no confidence interval is calculated.

- cl:

  Width of the symmetric confidence interval, relative to 1.

- B:

  Number of bootstrap replicates, or number of normal simulations from
  the distribution of the MLEs

- cores:

  Number of cores to use for bootstrapping using parallel processing.
  See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details.

- qlist:

  A list of transition intensity matrices, of length one greater than
  the length of `times`. Either this or a fitted model `x` must be
  supplied. No confidence intervals are available if (just) `qlist` is
  supplied.

- ...:

  Optional arguments to be passed to
  [`MatrixExp`](https://chjackson.github.io/msm/reference/MatrixExp.md)
  to control the method of computing the matrix exponential.

## Value

The matrix of estimated transition probabilities \\P(t)\\ for the time
interval `[t1, tn]`. That is, the probabilities of occupying state \\s\\
at time \\t_n\\ conditionally on occupying state \\r\\ at time \\t_1\\.
Rows correspond to "from-state" and columns to "to-state".

## Details

Suppose a multi-state model has been fitted, in which the transition
intensity matrix \\Q(x(t))\\ is modelled in terms of time-dependent
covariates \\x(t)\\. The transition probability matrix \\P(t_1, t_n)\\
for the time interval \\(t_1, \\\\ t_n)\\ cannot be calculated from the
estimated intensity matrix as \\\exp((t_n - t_1) Q)\\, because \\Q\\
varies within the interval \\t_1, t_n\\. However, if the covariates are
piecewise-constant, or can be approximated as piecewise-constant, then
we can calculate \\P(t_1, t_n)\\ by multiplying together individual
matrices \\P(t_i, \\\\ t\_{i+1}) = \exp((t\_{i+1} - t_i) Q)\\,
calculated over intervals where Q is constant:

\$\$P(t_1, t_n) = P(t_1, t_2) P(t_2, t_3)\ldots P(t\_{n-1}, \$\$\$\$
t_n)\$\$

## See also

[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r
if (FALSE) { # \dontrun{
## In a clinical study, suppose patients are given a placebo in the
## first 5 weeks, then they begin treatment 1 at 5 weeks, and
## a combination of treatments 1 and 2 from 10 weeks.
## Suppose a multi-state model x has been fitted for the patients'
## progress, with treat1 and treat2 as time dependent covariates.

## Cut points for when treatment covariate changes
times <- c(0, 5, 10)

## Indicators for which treatments are active in the four intervals
## defined by the three cut points
covariates <- list( list (treat1=0, treat2=0), list (treat1=0, treat2=0), list(treat1=1, treat2=0),
list(treat1=1, treat2=1) )

## Calculate transition probabilities from the start of the study to 15 weeks
pmatrix.piecewise.msm(x, 0, 15, times, covariates)
} # }
```
