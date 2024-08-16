#' Truncated Normal distribution
#' 
#' Density, distribution function, quantile function and random generation for
#' the truncated Normal distribution with mean equal to \code{mean} and
#' standard deviation equal to \code{sd} before truncation, and truncated on
#' the interval \code{[lower, upper]}.
#' 
#' The truncated normal distribution has density
#' 
#' \deqn{ f(x, \mu, \sigma) = \phi(x, \mu, \sigma) / (\Phi(u, \mu, \sigma) -
#' \Phi(l, \mu, \sigma)) }{ f(x, mu, sigma) = phi(x, mu, sigma) / (Phi(upper,
#' mu, sigma) - Phi(lower, mu, sigma)) } for \eqn{l <= x
#' <= u}{lower <= x <= upper}, and 0 otherwise.
#' 
#' \eqn{\mu}{mean} is the mean of the original Normal distribution before
#' truncation, \cr \eqn{\sigma}{sd} is the corresponding standard deviation,
#' \cr \eqn{u} is the upper truncation point, \cr \eqn{l} is the lower
#' truncation point, \cr \eqn{\phi(x)}{phi(x)} is the density of the
#' corresponding normal distribution, and \cr \eqn{\Phi(x)}{Phi(x)} is the
#' distribution function of the corresponding normal distribution.
#' 
#' If \code{mean} or \code{sd} are not specified they assume the default values
#' of \code{0} and \code{1}, respectively.
#' 
#' If \code{lower} or \code{upper} are not specified they assume the default
#' values of \code{-Inf} and \code{Inf}, respectively, corresponding to no
#' lower or no upper truncation.
#' 
#' Therefore, for example, \code{dtnorm(x)}, with no other arguments, is simply
#' equivalent to \code{dnorm(x)}.
#' 
#' Only \code{rtnorm} is used in the \code{msm} package, to simulate from
#' hidden Markov models with truncated normal distributions. This uses the
#' rejection sampling algorithms described by Robert (1995).
#' 
#' These functions are merely provided for completion, and are not optimized
#' for numerical stability or speed.  To fit a hidden Markov model with a
#' truncated Normal response distribution, use a \code{\link{hmmTNorm}}
#' constructor. See the \code{\link{hmm-dists}} help page for further details.
#'
#' @name tnorm
#' @aliases tnorm dtnorm ptnorm qtnorm rtnorm
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @param lower lower truncation point.
#' @param upper upper truncation point.
#' @param log logical; if TRUE, return log density or log hazard.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x],
#' otherwise, P[X > x].
#' @return \code{dtnorm} gives the density, \code{ptnorm} gives the
#' distribution function, \code{qtnorm} gives the quantile function, and
#' \code{rtnorm} generates random deviates.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{dnorm}}
#' @references Robert, C. P. Simulation of truncated normal variables.
#' Statistics and Computing (1995) 5, 121--125
#' @keywords distribution
#' @examples
#' 
#' x <- seq(50, 90, by=1)
#' plot(x, dnorm(x, 70, 10), type="l", ylim=c(0,0.06)) ## standard Normal distribution
#' lines(x, dtnorm(x, 70, 10, 60, 80), type="l")       ## truncated Normal distribution
#' 
NULL



### Truncated normal distribution

#' @rdname tnorm
#' @export
dtnorm <- function(x, mean=0, sd=1, lower=-Inf, upper=Inf, log=FALSE)
  {
      ret <- numeric(length(x))
      ret[x < lower | x > upper] <- if (log) -Inf else 0
      ret[upper < lower] <- NaN
      ind <- x >=lower & x <=upper
      if (any(ind)) {
          denom <- pnorm(upper, mean, sd) - pnorm(lower, mean, sd)
          xtmp <- dnorm(x, mean, sd, log)
          if (log) xtmp <- xtmp - log(denom) else xtmp <- xtmp/denom
          ret[x >=lower & x <=upper] <- xtmp[ind]
      }
      ret
  }

#' @rdname tnorm
#' @export
ptnorm <- function(q, mean=0, sd=1, lower=-Inf, upper=Inf, lower.tail=TRUE, log.p=FALSE)
  {
      ret <- numeric(length(q))
      if (lower.tail) {
          ret[q < lower] <- 0
          ret[q > upper] <- 1
      }
      else {
          ret[q < lower] <- 1
          ret[q > upper] <- 0
      }
      ret[upper < lower] <- NaN
      ind <- q >=lower & q <=upper
      if (any(ind)) {
          denom <- pnorm(upper, mean, sd) - pnorm(lower, mean, sd)
          if (lower.tail) qtmp <- pnorm(q, mean, sd) - pnorm(lower, mean, sd)
          else qtmp <- pnorm(upper, mean, sd) - pnorm(q, mean, sd)
          if (log.p) qtmp <- log(qtmp) - log(denom) else qtmp <- qtmp/denom
          ret[q >=lower & q <=upper] <- qtmp[ind]
      }
      ret
  }

#' @rdname tnorm
#' @export
qtnorm <- function(p, mean=0, sd=1, lower=-Inf, upper=Inf, lower.tail=TRUE, log.p=FALSE)
{
    qgeneric(ptnorm, p=p, mean=mean, sd=sd, lower=lower, upper=upper, lbound=lower, ubound=upper, lower.tail=lower.tail, log.p=log.p)
}

## Rejection sampling algorithm by Robert (Stat. Comp (1995), 5, 121-5)
## for simulating from the truncated normal distribution.

#' @rdname tnorm
#' @export
rtnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
    if (length(n) > 1)
        n <- length(n)
    mean <- rep(mean, length=n)
    sd <- rep(sd, length=n)
    lower <- rep(lower, length=n)
    upper <- rep(upper, length=n)
    ret <- numeric(n)
    ind <- seq(length.out=n)

    sdzero <- sd < .Machine$double.eps
    ## return the mean, unless mean is outside the range, then return nan 
    sdna <- sdzero & ((mean < lower) | (mean > upper))

    lower <- (lower - mean) / sd ## Algorithm works on mean 0, sd 1 scale
    upper <- (upper - mean) / sd
    nas <- is.na(mean) | is.na(sd) | is.na(lower) | is.na(upper) | sdna
    if (any(nas)) warning("NAs produced")
    ## Different algorithms depending on where upper/lower limits lie.
    alg <- ifelse(
                  ((lower > upper) | nas),
                  -1,# return NaN
                  ifelse(
                         sdzero, 
                         4, # SD zero, so set the sampled value to the mean. 
                         ifelse(
                         ((lower < 0 & upper == Inf) |
                          (lower == -Inf & upper > 0) |
                          (is.finite(lower) & is.finite(upper) & (lower < 0) & (upper > 0) & (upper-lower > sqrt(2*pi)))
                          ),
                         0, # standard "simulate from normal and reject if outside limits" method. Use if bounds are wide.
                         ifelse(
                                (lower >= 0 & (upper > lower + 2*sqrt(exp(1)) /
                                 (lower + sqrt(lower^2 + 4)) * exp((lower*2 - lower*sqrt(lower^2 + 4)) / 4))),
                                1, # rejection sampling with exponential proposal. Use if lower >> mean
                                ifelse(upper <= 0 & (-lower > -upper + 2*sqrt(exp(1)) /
                                       (-upper + sqrt(upper^2 + 4)) * exp((upper*2 - -upper*sqrt(upper^2 + 4)) / 4)),
                                       2, # rejection sampling with exponential proposal. Use if upper << mean.
                                       3))))) # rejection sampling with uniform proposal. Use if bounds are narrow and central.

    ind.nan <- ind[alg==-1]; ind.no <- ind[alg==0]; ind.expl <- ind[alg==1]; ind.expu <- ind[alg==2]; ind.u <- ind[alg==3]
    ind.sd0 <- ind[alg==4]; 
    ret[ind.nan] <- NaN
    ret[ind.sd0] <- 0  # SD zero, so set the sampled value to the mean.
    while (length(ind.no) > 0) {
        y <- rnorm(length(ind.no))
        done <- which(y >= lower[ind.no] & y <= upper[ind.no])
        ret[ind.no[done]] <- y[done]
        ind.no <- setdiff(ind.no, ind.no[done])
    }
    stopifnot(length(ind.no) == 0)
    while (length(ind.expl) > 0) {
        a <- (lower[ind.expl] + sqrt(lower[ind.expl]^2 + 4)) / 2
        z <- rexp(length(ind.expl), a) + lower[ind.expl]
        u <- runif(length(ind.expl))
        done <- which((u <= exp(-(z - a)^2 / 2)) & (z <= upper[ind.expl]))
        ret[ind.expl[done]] <- z[done]
        ind.expl <- setdiff(ind.expl, ind.expl[done])
    }
    stopifnot(length(ind.expl) == 0)
    while (length(ind.expu) > 0) {
        a <- (-upper[ind.expu] + sqrt(upper[ind.expu]^2 +4)) / 2
        z <- rexp(length(ind.expu), a) - upper[ind.expu]
        u <- runif(length(ind.expu))
        done <- which((u <= exp(-(z - a)^2 / 2)) & (z <= -lower[ind.expu]))
        ret[ind.expu[done]] <- -z[done]
        ind.expu <- setdiff(ind.expu, ind.expu[done])
    }
    stopifnot(length(ind.expu) == 0)
    while (length(ind.u) > 0) {
        z <- runif(length(ind.u), lower[ind.u], upper[ind.u])
        rho <- ifelse(lower[ind.u] > 0,
                      exp((lower[ind.u]^2 - z^2) / 2), ifelse(upper[ind.u] < 0,
                                                            exp((upper[ind.u]^2 - z^2) / 2),
                                                            exp(-z^2/2)))
        u <- runif(length(ind.u))
        done <- which(u <= rho)
        ret[ind.u[done]] <- z[done]
        ind.u <- setdiff(ind.u, ind.u[done])
    }
    stopifnot(length(ind.u) == 0)
    ret*sd + mean
}


