# Hidden Markov model constructors

These functions are used to specify the distribution of the response
conditionally on the underlying state in a hidden Markov model. A list
of these function calls, with one component for each state, should be
used for the `hmodel` argument to `msm`. The initial values for the
parameters of the distribution should be given as arguments. Note the
initial values should be supplied as literal values - supplying them as
variables is currently not supported.

## Usage

``` r
hmmCat(prob, basecat)

hmmIdent(x)

hmmUnif(lower, upper)

hmmNorm(mean, sd)

hmmLNorm(meanlog, sdlog)

hmmExp(rate)

hmmGamma(shape, rate)

hmmWeibull(shape, scale)

hmmPois(rate)

hmmBinom(size, prob)

hmmBetaBinom(size, meanp, sdp)

hmmNBinom(disp, prob)

hmmBeta(shape1, shape2)

hmmTNorm(mean, sd, lower = -Inf, upper = Inf)

hmmMETNorm(mean, sd, lower, upper, sderr, meanerr = 0)

hmmMEUnif(lower, upper, sderr, meanerr = 0)

hmmT(mean, scale, df)
```

## Arguments

- prob:

  (`hmmCat`) Vector of probabilities of observing category
  `1, 2, ...{}, length(prob)` respectively. Or the probability governing
  a binomial or negative binomial distribution.

- basecat:

  (`hmmCat`) Category which is considered to be the "baseline", so that
  during estimation, the probabilities are parameterised as
  probabilities relative to this baseline category. By default, the
  category with the greatest probability is used as the baseline.

- x:

  (`hmmIdent`) Code in the data which denotes the exactly-observed
  state.

- lower:

  (`hmmUnif,hmmTNorm,hmmMEUnif`) Lower limit for an Uniform or truncated
  Normal distribution.

- upper:

  (`hmmUnif,hmmTNorm,hmmMEUnif`) Upper limit for an Uniform or truncated
  Normal distribution.

- mean:

  (`hmmNorm,hmmLNorm,hmmTNorm`) Mean defining a Normal, or truncated
  Normal distribution.

- sd:

  (`hmmNorm,hmmLNorm,hmmTNorm`) Standard deviation defining a Normal, or
  truncated Normal distribution.

- meanlog:

  (`hmmNorm,hmmLNorm,hmmTNorm`) Mean on the log scale, for a log Normal
  distribution.

- sdlog:

  (`hmmNorm,hmmLNorm,hmmTNorm`) Standard deviation on the log scale, for
  a log Normal distribution.

- rate:

  (`hmmPois,hmmExp,hmmGamma`) Rate of a Poisson, Exponential or Gamma
  distribution (see [`dpois`](https://rdrr.io/r/stats/Poisson.html),
  [`dexp`](https://rdrr.io/r/stats/Exponential.html),
  [`dgamma`](https://rdrr.io/r/stats/GammaDist.html)).

- shape:

  (`hmmPois,hmmExp,hmmGamma`) Shape parameter of a Gamma or Weibull
  distribution (see [`dgamma`](https://rdrr.io/r/stats/GammaDist.html),
  [`dweibull`](https://rdrr.io/r/stats/Weibull.html)).

- scale:

  (`hmmGamma`) Scale parameter of a Gamma distribution (see
  [`dgamma`](https://rdrr.io/r/stats/GammaDist.html)), or unstandardised
  Student t distribution.

- size:

  Order of a Binomial distribution (see
  [`dbinom`](https://rdrr.io/r/stats/Binomial.html)).

- meanp:

  Mean outcome probability in a beta-binomial distribution

- sdp:

  Standard deviation describing the overdispersion of the outcome
  probability in a beta-binomial distribution

- disp:

  Dispersion parameter of a negative binomial distribution, also called
  `size` or `order`. (see
  [`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html)).

- shape1, shape2:

  First and second parameters of a beta distribution (see
  [`dbeta`](https://rdrr.io/r/stats/Beta.html)).

- sderr:

  (`hmmMETNorm,hmmUnif`) Standard deviation of the Normal measurement
  error distribution.

- meanerr:

  (`hmmMETNorm,hmmUnif`) Additional shift in the measurement error,
  fixed to 0 by default. This may be modelled in terms of covariates.

- df:

  Degrees of freedom of the Student t distribution.

## Value

Each function returns an object of class `hmmdist`, which is a list
containing information about the model. The only component which may be
useful to end users is `r`, a function of one argument `n` which returns
a random sample of size `n` from the given distribution.

## Details

`hmmCat` represents a categorical response distribution on the set
`1, 2, ...{}, length(prob)`. The Markov model with misclassification is
an example of this type of model. The categories in this case are (some
subset of) the underlying states.

The `hmmIdent` distribution is used for underlying states which are
observed exactly without error. For hidden Markov models with multiple
outcomes, (see
[`hmmMV`](https://chjackson.github.io/msm/reference/hmmMV.md)), the
outcome in the data which takes the special `hmmIdent` value must be the
first of the multiple outcomes.

`hmmUnif`, `hmmNorm`, `hmmLNorm`, `hmmExp`, `hmmGamma`, `hmmWeibull`,
`hmmPois`, `hmmBinom`, `hmmTNorm`, `hmmNBinom` and `hmmBeta` represent
Uniform, Normal, log-Normal, exponential, Gamma, Weibull, Poisson,
Binomial, truncated Normal, negative binomial and beta distributions,
respectively, with parameterisations the same as the default
parameterisations in the corresponding base R distribution functions.

`hmmT` is the Student t distribution with general mean \\\mu\\, scale
\\\sigma\\ and degrees of freedom `df`. The variance is \\\sigma^2
df/(df + 2)\\. Note the t distribution in base R
[`dt`](https://rdrr.io/r/stats/TDist.html) is a standardised one with
mean 0 and scale 1. These allow any positive (integer or non-integer)
`df`. By default, all three parameters, including `df`, are estimated
when fitting a hidden Markov model, but in practice, `df` might need to
be fixed for identifiability - this can be done using the `fixedpars`
argument to [`msm`](https://chjackson.github.io/msm/reference/msm.md).

The `hmmMETNorm` and `hmmMEUnif` distributions are truncated Normal and
Uniform distributions, but with additional Normal measurement error on
the response. These are generalisations of the distributions proposed by
Satten and Longini (1996) for modelling the progression of CD4 cell
counts in monitoring HIV disease. See
[`medists`](https://chjackson.github.io/msm/reference/medists.md) for
density, distribution, quantile and random generation functions for
these distributions. See also
[`tnorm`](https://chjackson.github.io/msm/reference/tnorm.md) for
density, distribution, quantile and random generation functions for the
truncated Normal distribution.

See the PDF manual `msm-manual.pdf` in the `doc` subdirectory for
algebraic definitions of all these distributions. New hidden Markov
model response distributions can be added to msm by following the
instructions in Section 2.17.1.

Parameters which can be modelled in terms of covariates, on the scale of
a link function, are as follows.

|                |                            |
|----------------|----------------------------|
| PARAMETER NAME | LINK FUNCTION              |
| `mean`         | identity                   |
| `meanlog`      | identity                   |
| `rate`         | log                        |
| `scale`        | log                        |
| `meanerr`      | identity                   |
| `meanp`        | logit                      |
| `prob`         | logit or multinomial logit |

Parameters `basecat, lower, upper, size, meanerr` are fixed at their
initial values. All other parameters are estimated while fitting the
hidden Markov model, unless the appropriate `fixedpars` argument is
supplied to `msm`.

For categorical response distributions `(hmmCat)` the outcome
probabilities initialized to zero are fixed at zero, and the probability
corresponding to `basecat` is fixed to one minus the sum of the
remaining probabilities. These remaining probabilities are estimated,
and can be modelled in terms of covariates via multinomial logistic
regression (relative to `basecat`).

## References

Satten, G.A. and Longini, I.M. Markov chains with measurement error:
estimating the 'true' course of a marker of the progression of human
immunodeficiency virus disease (with discussion) *Applied Statistics*
45(3): 275-309 (1996).

Jackson, C.H. and Sharples, L.D. Hidden Markov models for the onset and
progresison of bronchiolitis obliterans syndrome in lung transplant
recipients *Statistics in Medicine*, 21(1): 113–128 (2002).

Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W. and Couto,
E. Multi-state Markov models for disease progression with classification
error. *The Statistician*, 52(2): 193–209 (2003).

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
