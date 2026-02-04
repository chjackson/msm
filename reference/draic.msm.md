# Criteria for comparing two multi-state models with nested state spaces

A modification of Akaike's information criterion, and a leave-one-out
likelihood cross-validation criterion, for comparing the predictive
ability of two Markov multi-state models with nested state spaces. This
is evaluated based on the restricted or aggregated data which the models
have in common.

## Usage

``` r
draic.msm(
  msm.full,
  msm.coarse,
  likelihood.only = FALSE,
  information = c("expected", "observed"),
  tl = 0.95
)

drlcv.msm(
  msm.full,
  msm.coarse,
  tl = 0.95,
  cores = NULL,
  verbose = TRUE,
  outfile = NULL
)
```

## Arguments

- msm.full:

  Model on the bigger state space.

- msm.coarse:

  Model on the smaller state space.

  The two models must both be non-hidden Markov models without censored
  states.

  The two models must be fitted to the same datasets, except that the
  state space of the coarse model must be an aggregated version of the
  state space of the full model. That is, every state in the full
  dataset must correspond to a unique state in the coarse dataset. For
  example, for the full state variable `c(1,1,2,2,3,4)`, the
  corresponding coarse states could be `c(1,1,2,2,2,3)`, but not
  `c(1,2,3,4,4,4)`.

  The structure of allowed transitions in the coarse model must also be
  a collapsed version of the big model structure, but no check is
  currently made for this in the code.

  To use these functions, all objects which were used in the calls to
  fit `msm.full` and `msm.coarse` must be in the working environment,
  for example, datasets and definitions of transition matrices.

- likelihood.only:

  Don't calculate Hessians and trace term (DRAIC).

- information:

  Use observed or expected information in the DRAIC trace term. Expected
  is the default, and much faster, though is only available for models
  fitted to pure panel data (all `obstype=1` in the call to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), thus not
  exact transition times or exact death times)

- tl:

  Width of symmetric tracking interval, by default 0.95 for a 95%
  interval.

- cores:

  Number of processor cores to use in `drlcv` for cross-validation by
  parallel processing. Requires the doParallel package to be installed.
  If not specified, parallel processing is not used. If `cores` is set
  to the string `"default"`, the default methods of
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html) (on
  Windows) or
  [`registerDoParallel`](https://rdrr.io/pkg/doParallel/man/registerDoParallel.html)
  (on Unix-like) are used.

- verbose:

  Print intermediate results of each iteration of cross-validation to
  the console while running. May not work with parallel processing.

- outfile:

  Output file to print intermediate results of cross-validation. Useful
  to track execution speed when using parallel processing, where output
  to the console may not work.

## Value

A list containing \\D\_{RAIC}\\ (`draic.msm`) or \\D\_{RLCV}\\
(`drlcv.msm`), its component terms, and tracking intervals.

## Details

Note that standard AIC can be computed for one or more fitted `msm`
models `x,y,...` using
[`AIC`](https://rdrr.io/r/stats/AIC.html)`(x,y,...)`, and this can be
used to compare models fitted to the same data. `draic.msm` and
`drlcv.msm` are designed for models fitted to data with
differently-aggregated state spaces.

The difference in restricted AIC (Liquet and Commenges, 2011), as
computed by this function, is defined as

\$\$D\_{RAIC} = l(\gamma_n \|\mathbf{x}'' ) - l(\theta_n \|\mathbf{x}''
) + trace ( J(\theta_n \|\mathbf{x}'')J(\theta_n \|\mathbf{x})^{-1} -
J(\gamma_n \|\mathbf{x}'' )J(\gamma_n \|\mathbf{x}' )^{-1})\$\$

where \\\gamma\\ and \\\theta\\ are the maximum likelihood estimates of
the smaller and bigger models, fitted to the smaller and bigger data,
respectively.

\\l(\gamma_n \|x'')\\ represents the likelihood of the simpler model
evaluated on the restricted data.

\\l(\theta_n \|x'')\\ represents the likelihood of the complex model
evaluated on the restricted data. This is a hidden Markov model, with a
misclassification matrix and initial state occupancy probabilities as
described by Thom et al (2014).

\\J()\\ are the corresponding (expected or observed, as specified by the
user) information matrices.

\\\mathbf{x}\\ is the expanded data, to which the bigger model was
originally fitted, and \\\mathbf{x}'\\ is the data to which the smaller
model was originally fitted. \\\mathbf{x}''\\ is the restricted data
which the two models have in common. \\\mathbf{x}'' = \mathbf{x}'\\ in
this implementation, so the models are nested.

The difference in likelihood cross-validatory criteria (Liquet and
Commenges, 2011) is defined as

\$\$D\_{RLCV} = 1/n \sum\_{i=1}^n \log( h\_{X''}(x_i'' \| \gamma\_{-i})
/ g\_{X''}(x_i''\| \theta\_{-i}))\$\$

where \\\gamma\_{-i}\\ and \\\theta\_{-i}\\ are the maximum likelihood
estimates from the smaller and bigger models fitted to datasets with
subject \\i\\ left out, \\g()\\ and \\h()\\ are the densities of the
corresponding models, and \\x_i''\\ is the restricted data from subject
\\i\\.

Tracking intervals are analogous to confidence intervals, but not
strictly the same, since the quantity which D_RAIC aims to estimate, the
difference in expected Kullback-Leibler discrepancy for predicting a
replicate dataset, depends on the sample size. See the references.

Positive values for these criteria indicate the coarse model is
preferred, while negative values indicate the full model is preferred.

## References

Thom, H. and Jackson, C. and Commenges, D. and Sharples, L. (2015) State
selection in multistate models with application to quality of life in
psoriatic arthritis. Statistics In Medicine 34(16) 2381 - 2480.

Liquet, B. and Commenges D. (2011) Choice of estimators based on
different observations: Modified AIC and LCV criteria. Scandinavian
Journal of Statistics; 38:268-287.

## See also

[`logLik.msm`](https://chjackson.github.io/msm/reference/logLik.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>, H. H. Z. Thom
<howard.thom@bristol.ac.uk>
