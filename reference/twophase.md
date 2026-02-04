# Coxian phase-type distribution with two phases

Density, distribution, quantile functions and other utilities for the
Coxian phase-type distribution with two phases.

## Usage

``` r
d2phase(x, l1, mu1, mu2, log = FALSE)

p2phase(q, l1, mu1, mu2, lower.tail = TRUE, log.p = FALSE)

q2phase(p, l1, mu1, mu2, lower.tail = TRUE, log.p = FALSE)

r2phase(n, l1, mu1, mu2)

h2phase(x, l1, mu1, mu2, log = FALSE)
```

## Arguments

- x, q:

  vector of quantiles.

- l1:

  Intensity for transition between phase 1 and phase 2.

- mu1:

  Intensity for transition from phase 1 to exit.

- mu2:

  Intensity for transition from phase 2 to exit.

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

`d2phase` gives the density, `p2phase` gives the distribution function,
`q2phase` gives the quantile function, `r2phase` generates random
deviates, and `h2phase` gives the hazard.

## Details

This is the distribution of the time to reach state 3 in a
continuous-time Markov model with three states and transitions permitted
from state 1 to state 2 (with intensity \\\lambda_1\\) state 1 to state
3 (intensity \\\mu_1\\) and state 2 to state 3 (intensity \\\mu_2\\).
States 1 and 2 are the two "phases" and state 3 is the "exit" state.

The density is

\$\$f(t \| \lambda_1, \mu_1) = e^{-(\lambda_1+\mu_1)t}(\mu_1 +
(\lambda_1+\mu_1)\lambda_1 t)\$\$

if \\\lambda_1 + \mu_1 = \mu_2\\, and

\$\$f(t \| \lambda_1, \mu_1, \mu_2) =
\frac{(\lambda_1+\mu_1)e^{-(\lambda_1+\mu_1)t}(\mu_2-\mu_1) +
\mu_2\lambda_1e^{-\mu_2t}}{\lambda_1+\mu_1-\mu_2}\$\$

otherwise. The distribution function is

\$\$F(t \| \lambda_1, \mu_1) = 1 - e^{-(\lambda_1+\mu_1) t} (1 +
\lambda_1 t)\$\$

if \\\lambda_1 + \mu_1 = \mu_2\\, and

\$\$F(t \| \lambda_1, \mu_1, \mu_2) = 1 - \frac{e^{-(\lambda_1 +
\mu_1)t} (\mu_2 - \mu_1) + \lambda_1 e^{-\mu_2 t}}{ \lambda_1 + \mu_1 -
\mu_2}\$\$

otherwise. Quantiles are calculated by numerically inverting the
distribution function.

The mean is \\(1 + \lambda_1/\mu_2) / (\lambda_1 + \mu_1)\\.

The variance is \\(2 + 2\lambda_1(\lambda_1+\mu_1+ \mu_2)/\mu_2^2 - (1 +
\lambda_1/\mu_2)^2)/(\lambda_1+\mu_1)^2\\.

If \\\mu_1=\mu_2\\ it reduces to an exponential distribution with rate
\\\mu_1\\, and the parameter \\\lambda_1\\ is redundant. Or also if
\\\lambda_1=0\\.

The hazard at \\x=0\\ is \\\mu_1\\, and smoothly increasing if
\\\mu_1\<\mu_2\\. If \\\lambda_1 + \mu_1 \geq \mu_2\\ it increases to an
asymptote of \\\mu_2\\, and if \\\lambda_1 + \mu_1 \leq \mu_2\\ it
increases to an asymptote of \\\lambda_1 + \mu_1\\. The hazard is
decreasing if \\\mu_1\>\mu_2\\, to an asymptote of \\\mu_2\\.

## Alternative parameterisation

An individual following this distribution can be seen as coming from a
mixture of two populations:

1\) "short stayers" whose mean sojourn time is \\M_1 = \\\\
1/(\lambda_1+\mu_1)\\ and sojourn distribution is exponential with rate
\\\lambda_1 + \mu_1\\.

2\) "long stayers" whose mean sojourn time \\M_2 = \\\\
1/(\lambda_1+\mu_1) + 1/\mu_2\\ and sojourn distribution is the sum of
two exponentials with rate \\\lambda_1 + \\\\ \mu_1\\ and \\\mu_2\\
respectively. The individual is a "long stayer" with probability
\\p=\lambda_1/(\lambda_1 + \mu_1)\\.

Thus a two-phase distribution can be more intuitively parameterised by
the short and long stay means \\M_1 \< M_2\\ and the long stay
probability \\p\\. Given these parameters, the transition intensities
are \\\lambda_1=p/M_1\\, \\\mu_1=(1-p)/M_1\\, and \\\mu_2=1/(M_2-M_1)\\.
This can be useful for choosing intuitively reasonable initial values
for procedures to fit these models to data.

The hazard is increasing at least if \\M_2 \< 2M_1\\, and also only if
\\(M_2 - 2M_1)/(M_2 - M_1) \< p\\.

For increasing hazards with \\\lambda_1 + \mu_1 \leq \mu_2\\, the
maximum hazard ratio between any time \\t\\ and time 0 is \\1/(1-p)\\.

For increasing hazards with \\\lambda_1 + \mu_1 \geq \mu_2\\, the
maximum hazard ratio is \\M_1/((1-p)(M_2 - M_1))\\\\ M_1))\\. This is
the minimum hazard ratio for decreasing hazards.

## References

C. Dutang, V. Goulet and M. Pigeon (2008). actuar: An R Package for
Actuarial Science. Journal of Statistical Software, vol. 25, no. 7,
1-37. URL http://www.jstatsoft.org/v25/i07

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
