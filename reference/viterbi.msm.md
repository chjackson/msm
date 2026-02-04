# Calculate the probabilities of underlying states and the most likely path through them

For a fitted hidden Markov model, or a model with censored state
observations, the Viterbi algorithm recursively constructs the path with
the highest probability through the underlying states. The probability
of each hidden state is also computed for hidden Markov models, using
the forward-backward algorithm.

## Usage

``` r
viterbi.msm(x, normboot = FALSE, newdata = NULL)
```

## Arguments

- x:

  A fitted hidden Markov multi-state model, or a model with censored
  state observations, as produced by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md)

- normboot:

  If `TRUE`, then before running the algorithm, the maximum likelihood
  estimates of the model parameters are replaced by an alternative set
  of parameters drawn randomly from the asymptotic multivariate normal
  distribution of the MLEs.

- newdata:

  An optional data frame containing observations on which to construct
  the Viterbi path and forward-backward probabilities. It must be in the
  same format as the data frame used to fit `x`. If `NULL`, the data
  frame used to fit `x` is used.

## Value

A data frame with columns:

`subject` = subject identification numbers

`time` = times of observations

`observed` = corresponding observed states

`fitted` = corresponding fitted states found by Viterbi recursion. If
the model is not a hidden Markov model, and there are no censored state
observations, this is just the observed states.

For hidden Markov models, an additional matrix `pstate` is also returned
inside the data frame, giving the probability of each hidden state at
each point, conditionally on all the data. This is computed by the
forward/backward algorithm.

## References

Durbin, R., Eddy, S., Krogh, A. and Mitchison, G. *Biological sequence
analysis*, Cambridge University Press, 1998.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
