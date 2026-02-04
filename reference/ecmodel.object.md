# Developer documentation: model for covariates on misclassification probabilities

A list representing the model for covariates on misclassification
probabilities.

## Value

- npars:

  Number of covariate effect parameters. This is defined as the number
  of covariates on misclassification (with factors expanded as
  contrasts) multiplied by the number of allowed misclassifications in
  the model.

- ndpars:

  Number of distinct covariate effect parameters, as `npars`, but after
  any equality constraints have been applied.

- ncovs:

  Number of covariates on misclassification, with factors expanded as
  contrasts.

- constr:

  List of equality constraints on these covariate effects, as supplied
  in the `miscconstraint` argument to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covlabels:

  Names / labels of these covariates in the model matrix (see
  [`model.matrix.msm`](https://chjackson.github.io/msm/reference/model.frame.msm.md)).

- inits:

  Initial values for these covariate effects, as a vector formed from
  the `misccovinits` list supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- covmeans:

  Means of these covariates in the data (excluding data not required to
  fit the model, such as observations with missing data in other
  elements or subjects' last observations). This includes means of 0/1
  factor contrasts as well as continuous covariates (for historic
  reasons, which may not be sensible).

## See also

[`msm.object`](https://chjackson.github.io/msm/reference/msm.object.md).
