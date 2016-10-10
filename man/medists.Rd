\name{medists}
\alias{medists}
\alias{dmenorm}
\alias{pmenorm}
\alias{qmenorm}
\alias{rmenorm}
\alias{dmeunif}
\alias{pmeunif}
\alias{qmeunif}
\alias{rmeunif}
\title{Measurement error distributions}
\concept{Measurement error}
\description{
  Truncated Normal and Uniform distributions, where the response is also
  subject to a Normally distributed measurement error.
}
\usage{
     dmenorm(x, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0,
             log = FALSE)
     pmenorm(q, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0,
             lower.tail = TRUE, log.p = FALSE)
     qmenorm(p, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0,
             lower.tail = TRUE, log.p = FALSE)
     rmenorm(n, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0)
     dmeunif(x, lower=0, upper=1, sderr=0, meanerr=0, log = FALSE)
     pmeunif(q, lower=0, upper=1, sderr=0, meanerr=0, lower.tail = TRUE,
             log.p = FALSE)
     qmeunif(p, lower=0, upper=1, sderr=0, meanerr=0, lower.tail = TRUE,
             log.p = FALSE)
     rmeunif(n, lower=0, upper=1, sderr=0, meanerr=0)
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
  \item{sderr}{Standard deviation of measurement error distribution.}
  \item{meanerr}{Optional shift for the measurement error distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities \eqn{p} are given as
    \eqn{\log(p)}{log(p)}, or log density is returned.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[X <= x]},
    otherwise, \eqn{P[X > x]}.}
}
\value{
  \code{dmenorm}, \code{dmeunif} give the density, \code{pmenorm}, \code{pmeunif} give the distribution
  function, \code{qmenorm}, \code{qmeunif} give the quantile
  function, and \code{rmenorm}, \code{rmeunif} generate random
  deviates, for the Normal and Uniform versions respectively.
}
\details{
  The normal distribution with measurement error has density

  \deqn{
    \frac{\Phi(u, \mu_2, \sigma_3) - \Phi(l, \mu_2,
      \sigma_3)}{\Phi(u, \mu_0, \sigma_0) - \Phi(l, \mu_0, \sigma_0)}
    \phi(x, \mu_0 + \mu_\epsilon, \sigma_2)
  }{(Phi(upper, mu2, sigma3) - Phi(lower, mu2, sigma3)) / (Phi(upper, mean, sd) - Phi(lower, mean, sd)) * phi(x, mean + meanerr, sigma2)}
  where
  \deqn{\sigma_2^2 = \sigma_0^2 + \sigma_\epsilon^2,}{sigma2*sigma2 = sd*sd + sderr*sderr,}
  \deqn{\sigma_3 = \sigma_0 \sigma_\epsilon / \sigma_2,}{sigma3 = sd*sderr / sigma2,}
  \deqn{\mu_2 = (x - \mu_\epsilon) \sigma_0^2 + \mu_0
    \sigma_\epsilon^2,}{mu2 = (x - meanerr)*sd*sd + mean*sderr*sderr,}

  \eqn{\mu_0}{mean} is the mean of the original Normal distribution before
  truncation, \cr
  \eqn{\sigma_0}{sd} is the corresponding standard deviation, \cr
  \eqn{u}  is the upper truncation point, \cr
  \eqn{l} is the lower truncation point, \cr
  \eqn{\sigma_\epsilon}{sderr} is the standard deviation of the additional
  measurement error, \cr
  \eqn{\mu_\epsilon}{meanerr} is the mean of the measurement error (usually 0). \cr
  \eqn{\phi(x)}{phi(x)} is the density of the corresponding normal
  distribution, and \cr
  \eqn{\Phi(x)}{Phi(x)} is the distribution function of the corresponding normal
  distribution.

  The uniform distribution with measurement error has density

  \deqn{(\Phi(x, \mu_\epsilon+l, \sigma_\epsilon) - \Phi(x, \mu_\epsilon+u, \sigma_\epsilon))
    / (u - l)}{(Phi(x, meanerr+l, sderr) - Phi(x, meanerr+u, sderr)) / (upper - lower)}

  These are calculated from the original truncated Normal or Uniform density
  functions \eqn{f(. | \mu, \sigma, l, u)}{f(. | mu, sd)} as

  \deqn{\int f(y
    | \mu, \sigma, l, u) \phi(x, y
    + \mu_\epsilon, \sigma_\epsilon) dy}{integral f(y | mu, sd, l, u) phi(x, y + meanerr, sderr) dy}

  If \code{sderr} and \code{meanerr} are not specified they assume the
  default values of 0, representing no measurement error variance, and no
  constant shift in the measurement error, respectively.

  Therefore, for example with no other arguments, \code{dmenorm(x)}, is
  simply equivalent to \code{dtnorm(x)}, which in turn is
  equivalent to \code{dnorm(x)}.

  These distributions were used by Satten and Longini (1996) for
  CD4 cell counts conditionally on hidden Markov states of HIV
  infection, and later by Jackson and Sharples (2002) for FEV1
  measurements conditionally on states of chronic lung transplant
  rejection.

  These distribution functions are just provided for convenience, and are not
  optimised for numerical accuracy or speed.  To fit a hidden Markov model with these response distributions, use a
  \code{\link{hmmMETNorm}} or \code{\link{hmmMEUnif}} constructor. See
  the \code{\link{hmm-dists}} help page for further details.

}
\seealso{
  \code{\link{dnorm}}, \code{\link{dunif}}, \code{\link{dtnorm}}
}
\references{
  Satten, G.A. and Longini, I.M.  Markov chains with measurement error:
  estimating the 'true' course of a marker of the progression of human
  immunodeficiency virus disease (with discussion) \emph{Applied
    Statistics} 45(3): 275-309 (1996)

  Jackson, C.H. and Sharples, L.D. Hidden Markov models for the
  onset and progression of bronchiolitis obliterans syndrome in lung
  transplant recipients \emph{Statistics in Medicine}, 21(1): 113--128
  (2002).
}
\examples{
## what does the distribution look like?
x <- seq(50, 90, by=1)
plot(x, dnorm(x, 70, 10), type="l", ylim=c(0,0.06)) ## standard Normal
lines(x, dtnorm(x, 70, 10, 60, 80), type="l")       ## truncated Normal
## truncated Normal with small measurement error
lines(x, dmenorm(x, 70, 10, 60, 80, sderr=3), type="l")
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{distribution}
