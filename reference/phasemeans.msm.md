# Parameters of phase-type models in mixture form

Parameters of fitted two-phase models, in mixture model
parameterisation.

## Usage

``` r
phasemeans.msm(
  x,
  covariates = "mean",
  ci = c("none", "normal", "bootstrap"),
  cl = 0.95,
  B = 1000,
  cores = NULL
)
```

## Arguments

- x:

  A fitted multi-state model, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covariates:

  Covariate values, see
  [`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md).

- ci:

  If `"none"` (the default) no confidence intervals are calculated.
  Otherwise `"normal"`, or `"boot"` as described by
  [`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md).

- cl:

  Width of the symmetric confidence interval, relative to 1.

- B:

  Number of bootstrap replicates, or number of normal simulations from
  the distribution of the MLEs.

- cores:

  Number of cores to use for bootstrapping using parallel processing.
  See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details.

## Value

Matrix with one row for each state that has a two-phase distribution,
and three columns: the short-stay mean, long-stay mean and long-stay
probability. These are functions of the transition intensities of the
expanded hidden Markov model, defined in
[`d2phase`](https://chjackson.github.io/msm/reference/twophase.md).

## See also

[`d2phase`](https://chjackson.github.io/msm/reference/twophase.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
