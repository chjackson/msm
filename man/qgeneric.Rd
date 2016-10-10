\name{qgeneric}
\alias{qgeneric}
\title{Generic function to find quantiles of a distribution}
\description{
  Generic function to find the quantiles of a distribution,
  given the equivalent probability distribution function.
}
\usage{
qgeneric(pdist, p, ...)
}
\arguments{
  \item{pdist}{ Probability distribution function,
    for example, \code{\link{pnorm}} for the normal distribution,
    which must be defined in the current workspace.  This should
    accept and return vectorised parameters and values.
    It should also return the correct values for the entire real line,
    for example a positive distribution should have \code{pdist(x)==0}
    for \eqn{x<0}.
  }
  \item{p}{ Vector of probabilities to find the quantiles for. }

  \item{...}{ The remaining arguments define parameters of the
    distribution \code{pdist}.  These MUST be named explicitly.

    This may also contain the standard arguments \code{log.p} (logical;
    default \code{FALSE}, if \code{TRUE}, probabilities p are given as
    log(p)), and \code{lower.tail} (logical; if \code{TRUE} (default),
    probabilities are P[X <= x] otherwise, P[X > x].).

    If the distribution is bounded above or below, then this should contain
    arguments \code{lbound} and \code{ubound} respectively, and these
    will be returned if \code{p} is 0 or 1 respectively.  Defaults to
    \code{-Inf} and \code{Inf} respectively.
  }
}
\value{
  Vector of quantiles of the distribution at \code{p}.
}
\details{
  This function is intended to enable users to define \code{"q"} functions for
  new distributions, in cases where the distribution function \code{pdist} is
  available analytically, but the quantile function is not.

  It works by finding the root of the equation \eqn{h(q) = pdist(q) - p = 0}.
  Starting from the interval \eqn{(-1, 1)}, the interval width
  is expanded by 50\% until \eqn{h()} is of opposite sign at either end.
  The root is then found using \code{\link{uniroot}}.

  This assumes a suitably smooth, continuous distribution.

  An identical function is provided in the \pkg{flexsurv} package.
}
\examples{
qnorm(c(0.025, 0.975), 0, 1)
qgeneric(pnorm, c(0.025, 0.975), mean=0, sd=1) # must name the arguments
}
\author{
  Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\keyword{distribution}
