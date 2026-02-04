# Expected first passage time

Expected time until first reaching a particular state or set of states
in a Markov model.

## Usage

``` r
efpt.msm(
  x = NULL,
  qmatrix = NULL,
  tostate,
  start = "all",
  covariates = "mean",
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

- tostate:

  State, or set of states supplied as a vector, for which to estimate
  the first passage time into. Can be integer, or character matched to
  the row names of the Q matrix.

- start:

  Starting state (integer). By default (`start="all"`), this will return
  a vector of expected passage times from each state in turn.

  Alternatively, this can be used to obtain the expected first passage
  time from a *set* of states, rather than single states. To achieve
  this, `state` is set to a vector of weights, with length equal to the
  number of states in the model. These weights should be proportional to
  the probability of starting in each of the states in the desired set,
  so that weights of zero are supplied for other states. The function
  will calculate the weighted average of the expected passage times from
  each of the corresponding states.

- covariates:

  Covariate values defining the intensity matrix for the fitted model
  `x`, as supplied to
  [`qmatrix.msm`](https://chjackson.github.io/msm/reference/qmatrix.msm.md).

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

A vector of expected first passage times, or "hitting times", from each
state to the desired state.

## Details

The expected first passage times from each of a set of states
\\\mathbf{i}\\ to to the remaining set of states
\\\overline{\mathbf{i}}\\ in the state space, for a model with
transition intensity matrix \\Q\\, are

\$\$-Q\_{\mathbf{i},\mathbf{i}}^{-1} \mathbf{1}\$\$

where \\\mathbf{1}\\ is a vector of ones, and
\\Q\_{\mathbf{i},\mathbf{i}}\\ is the square subset of \\Q\\ pertaining
to states \\\mathbf{i}\\.

It is equal to the sum of mean sojourn times for all states between the
"from" and "to" states in a unidirectional model. If there is non-zero
chance of reaching an absorbing state before reaching `tostate`, then it
is infinite. It is trivially zero if the "from" state equals `tostate`.

This function currently only handles time-homogeneous Markov models. For
time-inhomogeneous models it will assume that \\Q\\ equals the average
intensity matrix over all times and observed covariates. Simulation
might be used to handle time dependence.

Note this is the *expectation* of first passage time, and the confidence
intervals are CIs for this mean, not predictive intervals for the first
passage time. The full distribution of the first passage time to a set
of states can be obtained by setting the rows of the intensity matrix
\\Q\\ corresponding to that set of states to zero to make a model where
those states are absorbing. The corresponding transition probability
matrix \\Exp(Qt)\\ then gives the probabilities of having hit or passed
that state by a time \\t\\ (see the example below). This is implemented
in
[`ppass.msm`](https://chjackson.github.io/msm/reference/ppass.msm.md).

## References

Norris, J. R. (1997) Markov Chains. Cambridge University Press.

## See also

[`sojourn.msm`](https://chjackson.github.io/msm/reference/sojourn.msm.md),
[`totlos.msm`](https://chjackson.github.io/msm/reference/totlos.msm.md),
[`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
             c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
efpt.msm(qmatrix=twoway4.q, tostate=3)
#> [1] Inf Inf   0 Inf
# given in state 1, expected time to reaching state 3 is infinite
# since may die (state 4) before entering state 3

# If we remove the death state from the model, EFPTs become finite
Q <- twoway4.q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
efpt.msm(qmatrix=Q, tostate=3)
#> [1] 14.0241 10.0241  0.0000

# Suppose we cannot die or regress while in state 2, can only go to state 3
Q <- twoway4.q; Q[2,4] <- Q[2,1] <- 0; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
efpt.msm(qmatrix=Q, tostate=3)
#> [1]      Inf 6.024096 0.000000      Inf
# The expected time from 2 to 3 now equals the mean sojourn time in 2.
-1/Q[2,2]
#> [1] 6.024096

# Calculate cumulative distribution of the first passage time
# into state 3 for the following three-state model
Q <- twoway4.q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
# Firstly form a model where the desired hitting state is absorbing
Q[3,] <- 0
MatrixExp(Q, t=10)[,3]
#> [1] 0.4790663 0.6501628 1.0000000
ppass.msm(qmatrix=Q, tot=10)
#>           [,1]     [,2]      [,3]
#> [1,] 1.0000000 0.917915 0.4790663
#> [2,] 0.4819236 1.000000 0.6501628
#> [3,] 0.0000000 0.000000 1.0000000
# Given in state 1 at time 0, P(hit 3 by time 10) = 0.479
MatrixExp(Q, t=50)[,3]  # P(hit 3 by time 50) = 0.98
#> [1] 0.9812676 0.9875017 1.0000000
ppass.msm(qmatrix=Q, tot=50)
#>      [,1]      [,2]      [,3]
#> [1,]  1.0 0.9999963 0.9812676
#> [2,]  0.5 1.0000000 0.9875017
#> [3,]  0.0 0.0000000 1.0000000

```
