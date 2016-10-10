\name{pexp}
\alias{pexp}
\alias{dpexp}
\alias{ppexp}
\alias{qpexp}
\alias{rpexp}
\title{Exponential distribution with piecewise-constant rate}
\description{
     Density, distribution function, quantile function and random
     generation for a generalisation of the exponential distribution, in
     which the rate changes at a series of times.
}
\usage{
     dpexp(x, rate=1, t=0, log = FALSE)
     ppexp(q, rate=1, t=0, lower.tail = TRUE, log.p = FALSE)
     qpexp(p, rate=1, t=0, lower.tail = TRUE, log.p = FALSE)
     rpexp(n, rate=1, t=0)
}
\arguments{

  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{rate}{vector of rates.}
  \item{t}{vector of the same length as \code{rate}, giving the times at
    which the rate changes. The first element of \code{t} should be 0,
    and \code{t} should be in increasing order.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as
    log(p), or log density is returned.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x],
    otherwise, P[X > x].}
}
\value{
     \code{dpexp} gives the density, \code{ppexp} gives the distribution
     function, \code{qpexp} gives the quantile function, and \code{rpexp}
     generates random deviates.
}
\details{
  Consider the exponential distribution with rates \eqn{r_1, \ldots,
    r_n}{r1,\dots, rn} changing at times \eqn{t_1, \ldots, t_n}{t1,
    \dots, tn}, with \eqn{t_1 = 0}{t1 = 0}. Suppose \eqn{t_k}{tk} is the
maximum \eqn{t_i}{ti} such that
  \eqn{t_i < x}{ti < x}.
  The density of this distribution at \eqn{x > 0} is \eqn{f(x)} for \eqn{k = 1},
  and \deqn{\prod_{i=1}^k (1 - F(t_{i} - t_{i-1}, r_i)) f(x - t_{k}, r_{k})}{\prod{i=1 \dots k} (1 - F(ti - t{i-1}, r{i-1})) f(x - tk, rk)}
  for k > 1.

  where \eqn{F()} and \eqn{f()} are the distribution and density
  functions of the standard exponential distribution.

  If \code{rate} is of length 1, this is just the standard exponential
  distribution.  Therefore, for example, \code{dpexp(x)}, with no other
  arguments, is simply equivalent to \code{dexp(x)}.

  Only \code{rpexp} is used in the \code{msm} package, to simulate
  from Markov processes with piecewise-constant intensities depending on
  time-dependent covariates.  These functions are merely provided for
  completion, and are not optimized for numerical stability or speed.
}
\seealso{
  \code{\link{dexp}}, \code{\link{sim.msm}}.
}
\examples{
x <- seq(0.1, 50, by=0.1)
rate <- c(0.1, 0.2, 0.05, 0.3)
t <- c(0, 10, 20, 30)
## standard exponential distribution
plot(x, dexp(x, 0.1), type="l")
## distribution with piecewise constant rate
lines(x, dpexp(x, rate, t), type="l", lty=2)
## standard exponential distribution
plot(x, pexp(x, 0.1), type="l")
## distribution with piecewise constant rate
lines(x, ppexp(x, rate, t), type="l", lty=2)
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{distribution}
