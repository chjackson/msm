# Passage probabilities

Probabilities of having visited each state by a particular time in a
continuous time Markov model.

## Usage

``` r
ppass.msm(
  x = NULL,
  qmatrix = NULL,
  tot,
  start = "all",
  covariates = "mean",
  piecewise.times = NULL,
  piecewise.covariates = NULL,
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

- qmatrix:

  Instead of `x`, you can simply supply a transition intensity matrix in
  `qmatrix`.

- tot:

  Finite time to forecast the passage probabilites for.

- start:

  Starting state (integer). By default (`start="all"`), this will return
  a matrix one row for each starting state.

  Alternatively, this can be used to obtain passage probabilities from a
  *set* of states, rather than single states. To achieve this, `state`
  is set to a vector of weights, with length equal to the number of
  states in the model. These weights should be proportional to the
  probability of starting in each of the states in the desired set, so
  that weights of zero are supplied for other states. The function will
  calculate the weighted average of the passage probabilities from each
  of the corresponding states.

- covariates:

  Covariate values defining the intensity matrix for the fitted model
  `x`, as supplied to
  [`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md).

- piecewise.times:

  For models with time-dependent covariates, this defines the cut points
  in time at which the transition intensity matrix changes. This is not
  required for models fitted with the `pci` option to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), which are
  handled automatically.

- piecewise.covariates:

  For models with time-dependent covariates, this is the list of
  covariates for each time period defined by `piecewise.times`, in the
  format documented for the `covariates` argument to
  [`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md).
  This is not required for models fitted with the `pci` option to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), which are
  handled automatically.

- ci:

  If `"normal"`, then calculate a confidence interval by simulating `B`
  random vectors from the asymptotic multivariate normal distribution
  implied by the maximum likelihood estimates (and covariance matrix) of
  the log transition intensities and covariate effects.

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

  Number of bootstrap replicates.

- cores:

  Number of cores to use for bootstrapping using parallel processing.
  See
  [`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md)
  for more details.

- ...:

  Arguments to pass to
  [`MatrixExp`](https://chjackson.github.io/msm/reference/MatrixExp.md).

## Value

A matrix whose \\r, s\\ entry is the probability of having visited state
\\s\\ at least once before time \\t\\, given the state at time \\0\\ is
\\r\\. The diagonal entries should all be 1.

## Details

The passage probabilities to state \\s\\ are computed by setting the
\\s\\th row of the transition intensity matrix \\Q\\ to zero, giving an
intensity matrix \\Q^\*\\ for a simplified model structure where state
\\s\\ is absorbing. The probabilities of passage are then equivalent to
row \\s\\ of the transition probability matrix \\Exp(tQ^\*)\\
([`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md))
under this simplified model for \\t=\\`tot`.

For time-inhomogenous models, this process is generalised by calculating
an intensity matrix for each time period, zeroing the appropriate row of
each, and calculating and multiplying transition probability matrices as
in
[`pmatrix.piecewise.msm`](https://chjackson.github.io/msm/reference/pmatrix.piecewise.msm.md).

Note this is different from the probability of occupying each state at
exactly time \\t\\, given by
[`pmatrix.msm`](https://chjackson.github.io/msm/reference/pmatrix.msm.md).
The passage probability allows for the possibility of having visited the
state before \\t\\, but then occupying a different state at \\t\\.

The mean of the passage distribution is the expected first passage time,
[`efpt.msm`](https://chjackson.github.io/msm/reference/efpt.msm.md).

## References

Norris, J. R. (1997) Markov Chains. Cambridge University Press.

## See also

[`efpt.msm`](https://chjackson.github.io/msm/reference/efpt.msm.md),
[`totlos.msm`](https://chjackson.github.io/msm/reference/totlos.msm.md),
[`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk> with contributions from
Jon Fintzi.

## Examples

``` r
Q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
           c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))

## ppass[1,2](t) converges to 0.5 with t, since given in state 1, the
## probability of going to the absorbing state 4 before visiting state
## 2 is 0.5, and the chance of still being in state 1 at t decreases.

ppass.msm(qmatrix=Q, tot=2)
#>            [,1]      [,2]       [,3]      [,4]
#> [1,] 1.00000000 0.3160603 0.04444417 0.3745829
#> [2,] 0.21468299 1.0000000 0.21468299 0.3091714
#> [3,] 0.04444417 0.3160603 1.00000000 0.3745829
#> [4,] 0.00000000 0.0000000 0.00000000 1.0000000
ppass.msm(qmatrix=Q, tot=20)
#>           [,1]      [,2]      [,3]      [,4]
#> [1,] 1.0000000 0.4999773 0.1990605 0.9862716
#> [2,] 0.3992305 1.0000000 0.3992305 0.9841246
#> [3,] 0.1990605 0.4999773 1.0000000 0.9862716
#> [4,] 0.0000000 0.0000000 0.0000000 1.0000000
ppass.msm(qmatrix=Q, tot=100)
#>      [,1] [,2] [,3] [,4]
#> [1,]  1.0  0.5  0.2    1
#> [2,]  0.4  1.0  0.4    1
#> [3,]  0.2  0.5  1.0    1
#> [4,]  0.0  0.0  0.0    1

Q <- Q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)

## Probability of about 1/2 of visiting state 3 by time 10.5, the
## median first passage time

ppass.msm(qmatrix=Q, tot=10.5)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.9275602 0.5000467
#> [2,] 0.6646619 1.0000000 0.6646619
#> [3,] 0.5000467 0.9275602 1.0000000

## Mean first passage time from state 2 to state 3 is 10.02: similar
## to the median

efpt.msm(qmatrix=Q, tostate=3)
#> [1] 14.0241 10.0241  0.0000
```
