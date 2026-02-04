# The delta method

Delta method for approximating the standard error of a transformation
\\g(X)\\ of a random variable \\X = (x_1, x_2, \ldots)\\, given
estimates of the mean and covariance matrix of \\X\\.

## Usage

``` r
deltamethod(g, mean, cov, ses = TRUE)
```

## Arguments

- g:

  A formula representing the transformation. The variables must be
  labelled `x1, x2,...{}` For example,

  `~ 1 / (x1 + x2)`

  If the transformation returns a vector, then a list of formulae
  representing (\\g_1, g_2, \ldots\\) can be provided, for example

  `list( ~ x1 + x2, ~ x1 / (x1 + x2) )`

- mean:

  The estimated mean of \\X\\

- cov:

  The estimated covariance matrix of \\X\\

- ses:

  If `TRUE`, then the standard errors of \\g_1(X), g_2(X),\ldots\\ are
  returned. Otherwise the covariance matrix of \\g(X)\\ is returned.

## Value

A vector containing the standard errors of \\g_1(X), g_2(X), \ldots\\ or
a matrix containing the covariance of \\g(X)\\.

## Details

The delta method expands a differentiable function of a random variable
about its mean, usually with a first-order Taylor approximation, and
then takes the variance. For example, an approximation to the covariance
matrix of \\g(X)\\ is given by

\$\$ Cov(g(X)) = g'(\mu) Cov(X) \[g'(\mu)\]^T \$\$

where \\\mu\\ is an estimate of the mean of \\X\\. This function uses
symbolic differentiation via
[`deriv`](https://rdrr.io/r/stats/deriv.html).

A limitation of this function is that variables created by the user are
not visible within the formula `g`. To work around this, it is necessary
to build the formula as a string, using functions such as `sprintf`,
then to convert the string to a formula using `as.formula`. See the
example below.

If you can spare the computational time, bootstrapping is a more
accurate method of calculating confidence intervals or standard errors
for transformations of parameters. See
[`boot.msm`](https://chjackson.github.io/msm/reference/boot.msm.md).
Simulation from the asymptotic distribution of the MLEs (see e.g. Mandel
2013) is also a convenient alternative.

## References

Oehlert, G. W. (1992) *A note on the delta method*. American
Statistician 46(1).

Mandel, M. (2013) *Simulation based confidence intervals for functions
with complicated derivatives.* The American Statistician 67(2):76-81.

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r

## Simple linear regression, E(y) = alpha + beta x 
x <- 1:100
y <- rnorm(100, 4*x, 5)
toy.lm <- lm(y ~ x)
estmean <- coef(toy.lm)
estvar <- summary(toy.lm)$cov.unscaled * summary(toy.lm)$sigma^2

## Estimate of (1 / (alphahat + betahat))
1 / (estmean[1] + estmean[2])
#> (Intercept) 
#>   0.2895236 
## Approximate standard error
deltamethod (~ 1 / (x1 + x2), estmean, estvar) 
#> [1] 0.08478445

## We have a variable z we would like to use within the formula.
z <- 1
## deltamethod (~ z / (x1 + x2), estmean, estvar) will not work.
## Instead, build up the formula as a string, and convert to a formula.
form <- sprintf("~ %f / (x1 + x2)", z)
form
#> [1] "~ 1.000000 / (x1 + x2)"
deltamethod(as.formula(form), estmean, estvar)
#> [1] 0.08478445

```
