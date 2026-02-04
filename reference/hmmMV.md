# Multivariate hidden Markov models

Constructor for a a multivariate hidden Markov model (HMM) where each of
the `n` variables observed at the same time has a (potentially
different) standard univariate distribution conditionally on the
underlying state. The `n` outcomes are independent conditionally on the
hidden state.

## Usage

``` r
hmmMV(...)
```

## Arguments

- ...:

  The number of arguments supplied should equal the maximum number of
  observations made at one time. Each argument represents the univariate
  distribution of that outcome conditionally on the hidden state, and
  should be the result of calling a univariate hidden Markov model
  constructor (see
  [`hmm-dists`](https://chjackson.github.io/msm/reference/hmm-dists.md)).

## Value

A list of objects, each of class `hmmdist` as returned by the univariate
HMM constructors documented in
[`hmm-dists`](https://chjackson.github.io/msm/reference/hmm-dists.md).
The whole list has class `hmmMVdist`, which inherits from `hmmdist`.

## Details

If a particular state in a HMM has such an outcome distribution, then a
call to `hmmMV` is supplied as the corresponding element of the `hmodel`
argument to [`msm`](https://chjackson.github.io/msm/reference/msm.md).
See Example 2 below.

A multivariate HMM where multiple outcomes at the same time are
generated from the *same* distribution is specified in the same way as
the corresponding univariate model, so that `hmmMV` is not required. The
outcome data are simply supplied as a matrix instead of a vector. See
Example 1 below.

The outcome data for such models are supplied as a matrix, with number
of columns equal to the maximum number of arguments supplied to the
`hmmMV` calls for each state. If some but not all of the variables are
missing (`NA`) at a particular time, then the observed data at that time
still contribute to the likelihood. The missing data are assumed to be
missing at random. The Viterbi algorithm may be used to predict the
missing values given the fitted model and the observed data.

Typically the outcome model for each state will be from the same family
or set of families, but with different parameters. Theoretically,
different numbers of distributions may be supplied for different states.
If a particular state has fewer outcomes than the maximum, then the data
for that state are taken from the first columns of the response data
matrix. However this is not likely to be a useful model, since the
number of observations will probably give information about the
underlying state, violating the missing at random assumption.

Models with outcomes that are dependent conditionally on the hidden
state (e.g. correlated multivariate normal observations) are not
currently supported.

## References

Jackson, C. H., Su, L., Gladman, D. D. and Farewell, V. T. (2015) On
modelling minimal disease activity. Arthritis Care and Research (early
view).

## See also

[`hmm-dists`](https://chjackson.github.io/msm/reference/hmm-dists.md),[`msm`](https://chjackson.github.io/msm/reference/msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r
## Simulate data from a Markov model 
nsubj <- 30; nobspt <- 5
sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt),
                     time = seq(0, 20, length=nobspt))
set.seed(1)
two.q <- rbind(c(-0.1, 0.1), c(0, 0))
dat <- simmulti.msm(sim.df[,1:2], qmatrix=two.q, drop.absorb=FALSE)

### EXAMPLE 1
## Generate two observations at each time from the same outcome
## distribution:
## Bin(40, 0.1) for state 1, Bin(40, 0.5) for state 2
dat$obs1[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.1)
dat$obs2[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.1)
dat$obs1[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.5)
dat$obs2[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.5)
dat$obs <- cbind(obs1 = dat$obs1, obs2 = dat$obs2)

## Fitted model should approximately recover true parameters 
msm(obs ~ time, subject=subject, data=dat, qmatrix=two.q,
    hmodel = list(hmmBinom(size=40, prob=0.2),
                  hmmBinom(size=40, prob=0.2)))
#> Warning: Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite.
#> 
#> Call:
#> msm(formula = obs ~ time, subject = subject, data = dat, qmatrix = two.q,     hmodel = list(hmmBinom(size = 40, prob = 0.2), hmmBinom(size = 40,         prob = 0.2)))
#> 
#> Optimisation probably not converged to the maximum likelihood.
#> optim() reported convergence but estimated Hessian not positive-definite.
#> 
#> -2 * log-likelihood:  3706.02 

### EXAMPLE 2
## Generate two observations at each time from different
## outcome distributions:
## Bin(40, 0.1) and Bin(40, 0.2) for state 1, 
dat$obs1 <- dat$obs2 <- NA
dat$obs1[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.1)
dat$obs2[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.2)

## Bin(40, 0.5) and Bin(40, 0.6) for state 2
dat$obs1[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.6)
dat$obs2[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.5)
dat$obs <- cbind(obs1 = dat$obs1, obs2 = dat$obs2)

## Fitted model should approximately recover true parameters 
msm(obs ~ time, subject=subject, data=dat, qmatrix=two.q,   
    hmodel = list(hmmMV(hmmBinom(size=40, prob=0.3),
                        hmmBinom(size=40, prob=0.3)),                 
                 hmmMV(hmmBinom(size=40, prob=0.3),
                       hmmBinom(size=40, prob=0.3))),
    control=list(maxit=10000))
#> 
#> Call:
#> msm(formula = obs ~ time, subject = subject, data = dat, qmatrix = two.q,     hmodel = list(hmmMV(hmmBinom(size = 40, prob = 0.3), hmmBinom(size = 40,         prob = 0.3)), hmmMV(hmmBinom(size = 40, prob = 0.3),         hmmBinom(size = 40, prob = 0.3))), control = list(maxit = 10000))
#> 
#> Maximum likelihood estimates
#> 
#> Transition intensities
#>                   Baseline                    
#> State 1 - State 1 -0.09458 (-0.13940,-0.06416)
#> State 1 - State 2  0.09458 ( 0.06416, 0.13940)
#> 
#> Hidden Markov model, 2 states
#> State 1 
#> Outcome 1 - binomial distribution
#>        Estimate       LCL       UCL
#> size 40.0000000        NA        NA
#> prob  0.1054793 0.0948457 0.1171507
#> 
#> Outcome 2 - binomial distribution
#>        Estimate       LCL       UCL
#> size 40.0000000        NA        NA
#> prob  0.1958908 0.1818942 0.2106871
#> 
#> State 2 
#> Outcome 1 - binomial distribution
#>        Estimate       LCL      UCL
#> size 40.0000000        NA       NA
#> prob  0.5954564 0.5780108 0.612664
#> 
#> Outcome 2 - binomial distribution
#>        Estimate       LCL      UCL
#> size 40.0000000        NA       NA
#> prob  0.5110363 0.4933762 0.528669
#> 
#> 
#> -2 * log-likelihood:  1523.652 
```
