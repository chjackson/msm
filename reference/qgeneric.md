# Generic function to find quantiles of a distribution

Generic function to find the quantiles of a distribution, given the
equivalent probability distribution function.

## Usage

``` r
qgeneric(pdist, p, special = NULL, ...)
```

## Arguments

- pdist:

  Probability distribution function, for example,
  [`pnorm`](https://rdrr.io/r/stats/Normal.html) for the normal
  distribution, which must be defined in the current workspace. This
  should accept and return vectorised parameters and values. It should
  also return the correct values for the entire real line, for example a
  positive distribution should have `pdist(x)==0` for \\x\<0\\.

- p:

  Vector of probabilities to find the quantiles for.

- special:

  Vector of character strings naming arguments of the distribution
  function that should not be vectorised over. Used, for example, for
  the `rate` and `t` arguments in
  [`qpexp`](https://chjackson.github.io/msm/reference/pexp.md).

- ...:

  The remaining arguments define parameters of the distribution `pdist`.
  These MUST be named explicitly.

  This may also contain the standard arguments `log.p` (logical; default
  `FALSE`, if `TRUE`, probabilities p are given as log(p)), and
  `lower.tail` (logical; if `TRUE` (default), probabilities are P\[X \<=
  x\] otherwise, P\[X \> x\].).

  If the distribution is bounded above or below, then this should
  contain arguments `lbound` and `ubound` respectively, and these will
  be returned if `p` is 0 or 1 respectively. Defaults to `-Inf` and
  `Inf` respectively.

## Value

Vector of quantiles of the distribution at `p`.

## Details

This function is intended to enable users to define `"q"` functions for
new distributions, in cases where the distribution function `pdist` is
available analytically, but the quantile function is not.

It works by finding the root of the equation \\h(q) = pdist(q) - p =
0\\. Starting from the interval \\(-1, 1)\\, the interval width is
expanded by 50% until \\h()\\ is of opposite sign at either end. The
root is then found using
[`uniroot`](https://rdrr.io/r/stats/uniroot.html).

This assumes a suitably smooth, continuous distribution.

An identical function is provided in the flexsurv package.

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>

## Examples

``` r
qnorm(c(0.025, 0.975), 0, 1)
#> [1] -1.959964  1.959964
qgeneric(pnorm, c(0.025, 0.975), mean=0, sd=1) # must name the arguments
#> [1] -1.959964  1.959964
```
