\name{deltamethod}
\title{The delta method}
\alias{deltamethod}
\concept{Delta method}
\description{
  Delta method for approximating the standard error of a transformation
  \eqn{g(X)} of a random variable \eqn{X = (x_1, x_2, \ldots)}{X = (x1, x2, \ldots)}, given estimates of the
  mean and covariance matrix of \eqn{X}.
}
\usage{
deltamethod(g, mean, cov, ses=TRUE)
}
\arguments{
  \item{g}{A formula representing the transformation. The variables must
    be labelled \code{x1, x2,\dots} For example,
    
    \code{~ 1 / (x1 + x2)}

    If the transformation returns a vector, then a list of formulae
    representing (\eqn{g_1, g_2, \ldots}{g1, g2, \ldots}) can be
    provided, for example

    \code{list( ~ x1 + x2,    ~ x1 / (x1 + x2) )}}

    \item{mean}{The estimated mean of \eqn{X}}

    \item{cov}{The estimated covariance matrix of \eqn{X}}

    \item{ses}{If \code{TRUE}, then the standard errors of \eqn{g_1(X), g_2(X),\ldots}{g1(X), g2(X),\ldots} are
      returned. Otherwise the covariance matrix of \eqn{g(X)} is returned.}
  }
\value{
  A vector containing the standard errors of \eqn{g_1(X), g_2(X), \ldots}{g1(X), g2(X), \ldots}
  or a matrix containing the covariance of \eqn{g(X)}.
}
\details{
  The delta method expands a differentiable function of a random
  variable about its mean, usually with a first-order Taylor approximation,
  and then takes the variance. For example, an approximation to
  the covariance matrix of \eqn{g(X)} is given by
  
  \deqn{ Cov(g(X)) = g'(\mu) Cov(X) [g'(\mu)]^T }{ Cov(g(X)) = g'(mu) Cov(X) [g'(mu)]^T }

  where \eqn{\mu}{mu} is an estimate of the mean of \eqn{X}.  This
  function uses symbolic differentiation via \code{\link{deriv}}.

  A limitation of this function is that variables created by the user
  are not visible within the formula \code{g}.  To work around this, it
  is necessary to build the formula as a string, using functions such as
  \code{sprintf}, then to convert the string to a formula using
  \code{as.formula}.  See the example below. 

  If you can spare the computational time, bootstrapping is a more  
  accurate method of calculating confidence intervals or standard errors
  for transformations of parameters. See \code{\link{boot.msm}}.
  Simulation from the asymptotic distribution of the MLEs (see e.g. Mandel 2013)
  is also a convenient alternative.
}
\references{
  Oehlert, G. W. (1992) \emph{A note on the delta method}. American
  Statistician 46(1).

  Mandel, M. (2013) \emph{Simulation based confidence intervals for
  functions with complicated derivatives.} The American Statistician 67(2):76-81.
}
\examples{

## Simple linear regression, E(y) = alpha + beta x 
x <- 1:100
y <- rnorm(100, 4*x, 5)
toy.lm <- lm(y ~ x)
estmean <- coef(toy.lm)
estvar <- summary(toy.lm)$cov.unscaled * summary(toy.lm)$sigma^2

## Estimate of (1 / (alphahat + betahat))
1 / (estmean[1] + estmean[2])
## Approximate standard error
deltamethod (~ 1 / (x1 + x2), estmean, estvar) 

## We have a variable z we would like to use within the formula.
z <- 1
## deltamethod (~ z / (x1 + x2), estmean, estvar) will not work.
## Instead, build up the formula as a string, and convert to a formula.
form <- sprintf("~ \%f / (x1 + x2)", z)
form
deltamethod(as.formula(form), estmean, estvar)

}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{math}
