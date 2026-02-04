# Score residuals

Score residuals for detecting outlying subjects.

## Usage

``` r
scoreresid.msm(x, plot = FALSE)
```

## Arguments

- x:

  A fitted multi-state model, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- plot:

  If `TRUE`, display a simple plot of the residuals in subject order,
  labelled by subject identifiers

## Value

Vector of the residuals, named by subject identifiers.

## Details

The score residual for a single subject is

\$\$U(\theta)^T I(\theta)^{-1} U(\theta)\$\$

where \\U(\theta)\\ is the vector of first derivatives of the
log-likelihood for that subject at maximum likelihood estimates
\\\theta\\, and \\I(\theta)\\ is the observed Fisher information matrix,
that is, the matrix of second derivatives of minus the log-likelihood
for that subject at theta.

Subjects with a higher influence on the maximum likelihood estimates
will have higher score residuals.

These are only available for models with analytic derivatives (which
includes all non-hidden and most hidden Markov models).

## Author

Andrew Titman <a.titman@lancaster.ac.uk> (theory), Chris Jackson
<chris.jackson@mrc-bsu.cam.ac.uk> (code)
