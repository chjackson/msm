\name{tnorm}
\alias{tnorm}
\alias{dtnorm}
\alias{ptnorm}
\alias{qtnorm}
\alias{rtnorm}
\title{Truncated Normal distribution}
\description{
     Density, distribution function, quantile function and random
     generation for the truncated Normal distribution with mean equal to \code{mean}
     and standard deviation equal to \code{sd} before truncation, and
     truncated on the interval \code{[lower, upper]}.
}
\usage{
     dtnorm(x, mean=0, sd=1, lower=-Inf, upper=Inf, log = FALSE)
     ptnorm(q, mean=0, sd=1, lower=-Inf, upper=Inf,
            lower.tail = TRUE, log.p = FALSE)
     qtnorm(p, mean=0, sd=1, lower=-Inf, upper=Inf,
            lower.tail = TRUE, log.p = FALSE)
     rtnorm(n, mean=0, sd=1, lower=-Inf, upper=Inf)
}
\arguments{

  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
  \item{lower}{lower truncation point.}
  \item{upper}{upper truncation point.}
  \item{log}{logical; if TRUE, return log density or log hazard.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x],
    otherwise, P[X > x].}
}
\value{
     \code{dtnorm} gives the density, \code{ptnorm} gives the distribution
     function, \code{qtnorm} gives the quantile function, and \code{rtnorm}
     generates random deviates.
}
\details{
  The truncated normal distribution has density

  \deqn{  f(x, \mu, \sigma) = \phi(x, \mu, \sigma) / (\Phi(u, \mu, \sigma) - \Phi(l, \mu, \sigma))
  }{  f(x, mu, sigma) = phi(x, mu, sigma) / (Phi(upper, mu, sigma) - Phi(lower, mu, sigma))
    }
    for \eqn{l <= x <= u}{lower <= x <= upper}, and 0 otherwise.

  \eqn{\mu}{mean} is the mean of the original Normal distribution before
  truncation, \cr
  \eqn{\sigma}{sd} is the corresponding standard deviation, \cr
  \eqn{u}  is the upper truncation point, \cr
  \eqn{l} is the lower truncation point, \cr
    \eqn{\phi(x)}{phi(x)} is the density of the corresponding normal
    distribution, and \cr
    \eqn{\Phi(x)}{Phi(x)} is the distribution function of the corresponding normal
    distribution.

  If \code{mean} or \code{sd} are not specified they assume the default values
  of \code{0} and \code{1}, respectively.

  If \code{lower} or \code{upper} are not specified they assume the default values
  of \code{-Inf} and \code{Inf}, respectively, corresponding to no
  lower or no upper truncation.

  Therefore, for example, \code{dtnorm(x)}, with no other arguments, is
  simply equivalent to \code{dnorm(x)}.

  Only \code{rtnorm} is used in the \code{msm} package, to simulate
  from hidden Markov models with truncated normal
  distributions. This uses the rejection sampling algorithms described
  by Robert (1995).

  These functions are merely provided for completion,
  and are not optimized for numerical stability or speed.  To fit a hidden Markov
  model with a truncated Normal response distribution, use a
  \code{\link{hmmTNorm}} constructor. See the \code{\link{hmm-dists}}
  help page for further details.

}
\references{
 Robert, C. P. Simulation of truncated normal variables. Statistics and
 Computing (1995) 5, 121--125
}
\seealso{
  \code{\link{dnorm}}
}
\examples{
x <- seq(50, 90, by=1)
plot(x, dnorm(x, 70, 10), type="l", ylim=c(0,0.06)) ## standard Normal distribution
lines(x, dtnorm(x, 70, 10, 60, 80), type="l")       ## truncated Normal distribution
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{distribution}
