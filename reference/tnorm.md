# Truncated Normal distribution

Density, distribution function, quantile function and random generation
for the truncated Normal distribution with mean equal to `mean` and
standard deviation equal to `sd` before truncation, and truncated on the
interval `[lower, upper]`.

## Usage

``` r
dtnorm(x, mean = 0, sd = 1, lower = -Inf, upper = Inf, log = FALSE)

ptnorm(
  q,
  mean = 0,
  sd = 1,
  lower = -Inf,
  upper = Inf,
  lower.tail = TRUE,
  log.p = FALSE
)

qtnorm(
  p,
  mean = 0,
  sd = 1,
  lower = -Inf,
  upper = Inf,
  lower.tail = TRUE,
  log.p = FALSE
)

rtnorm(n, mean = 0, sd = 1, lower = -Inf, upper = Inf)
```

## Arguments

- x, q:

  vector of quantiles.

- mean:

  vector of means.

- sd:

  vector of standard deviations.

- lower:

  lower truncation point.

- upper:

  upper truncation point.

- log:

  logical; if TRUE, return log density or log hazard.

- lower.tail:

  logical; if TRUE (default), probabilities are P\[X \<= x\], otherwise,
  P\[X \> x\].

- log.p:

  logical; if TRUE, probabilities p are given as log(p).

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dtnorm` gives the density, `ptnorm` gives the distribution function,
`qtnorm` gives the quantile function, and `rtnorm` generates random
deviates.

## Details

The truncated normal distribution has density

\$\$ f(x, \mu, \sigma) = \phi(x, \mu, \sigma) / (\Phi(u, \mu, \sigma) -
\Phi(l, \mu, \sigma)) \$\$ for \\l \<= x \<= u\\, and 0 otherwise.

\\\mu\\ is the mean of the original Normal distribution before
truncation,  
\\\sigma\\ is the corresponding standard deviation,  
\\u\\ is the upper truncation point,  
\\l\\ is the lower truncation point,  
\\\phi(x)\\ is the density of the corresponding normal distribution,
and  
\\\Phi(x)\\ is the distribution function of the corresponding normal
distribution.

If `mean` or `sd` are not specified they assume the default values of
`0` and `1`, respectively.

If `lower` or `upper` are not specified they assume the default values
of `-Inf` and `Inf`, respectively, corresponding to no lower or no upper
truncation.

Therefore, for example, `dtnorm(x)`, with no other arguments, is simply
equivalent to `dnorm(x)`.

Only `rtnorm` is used in the `msm` package, to simulate from hidden
Markov models with truncated normal distributions. This uses the
rejection sampling algorithms described by Robert (1995).

These functions are merely provided for completion, and are not
optimized for numerical stability or speed. To fit a hidden Markov model
with a truncated Normal response distribution, use a
[`hmmTNorm`](https://chjackson.github.io/msm/reference/hmm-dists.md)
constructor. See the
[`hmm-dists`](https://chjackson.github.io/msm/reference/hmm-dists.md)
help page for further details.

## References

Robert, C. P. Simulation of truncated normal variables. Statistics and
Computing (1995) 5, 121–125

## See also

[`dnorm`](https://rdrr.io/r/stats/Normal.html)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r
x <- seq(50, 90, by=1)
plot(x, dnorm(x, 70, 10), type="l", ylim=c(0,0.06)) ## standard Normal distribution
lines(x, dtnorm(x, 70, 10, 60, 80), type="l")       ## truncated Normal distribution

```
