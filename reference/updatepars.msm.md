# Update the maximum likelihood estimates in a fitted model object.

Update the maximum likelihood estimates in a fitted model object.
Intended for developer use only.

## Usage

``` r
updatepars.msm(x, pars)
```

## Arguments

- x:

  A fitted multi-state model object, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- pars:

  Vector of new parameters, in their untransformed real-line
  parameterisations, to substitute for the maximum likelihood estimates
  corresponding to those in the `estimates` component of the fitted
  model object
  ([`msm.object`](https://chjackson.github.io/msm/reference/msm.object.md)).
  The order of the parameters is documented in
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), argument
  `fixedpars`.

## Value

An updated [`msm`](https://chjackson.github.io/msm/reference/msm.md)
model object with the updated maximum likelihood estimates, but with the
covariances / standard errors unchanged.

Point estimates from output functions such as
[`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md),
[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md),
or any related function, can then be evaluated with the new parameters,
and at arbitrary covariate values.

This function is used, for example, when computing confidence intervals
from
[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md),
and related functions, using the `ci="normal"` method.

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
