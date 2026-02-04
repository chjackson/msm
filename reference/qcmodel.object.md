# Developer documentation: model for covariates on transition intensities

A list representing the model for covariates on transition intensities

## Value

- npars:

  Number of covariate effect parameters. This is defined as the number
  of covariates on intensities (with factors expanded as contrasts)
  multiplied by the number of allowed transitions in the model.

  Note if [`msm`](https://chjackson.github.io/msm/reference/msm.md) was
  called with `covariates` set to a list of different covariates for
  different intensities, then this will include covariate effects that
  are implicitly defined as zero by this list. The information in
  [`paramdata`](https://chjackson.github.io/msm/reference/paramdata.object.md)
  objects can be used to identify wich ones are fixed at zero.

  This also includes any `timeperiod` covariates in a time-inhomogeneous
  model defined by the `pci` option to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- ndpars:

  Number of distinct covariate effect parameters, as `npars`, but after
  any equality constraints have been applied.

- ncovs:

  Number of covariates on intensities, with factors expanded as
  contrasts.

- constr:

  List of equality constraints on these covariate effects, as supplied
  in the `constraint` argument to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covlabels:

  Names / labels of these covariates in the model matrix (see
  [`model.matrix.msm`](https://chjackson.github.io/msm/reference/model.frame.msm.md)).

- inits:

  Initial values for these covariate effects, as a vector formed from
  the `covinits` list supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covmeans:

  Means of these covariates in the data (excluding data not required to
  fit the model, such as observations with missing data in other
  elements or subjects' last observations). This includes means of 0/1
  factor contrasts as well as continuous covariates (for historic
  reasons, which may not be sensible).

## See also

[`msm.object`](https://chjackson.github.io/msm/reference/msm.object.md).
