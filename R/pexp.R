
#' Exponential distribution with piecewise-constant rate
#' 
#' Density, distribution function, quantile function and random generation for
#' a generalisation of the exponential distribution, in which the rate changes
#' at a series of times.
#' 
#' Consider the exponential distribution with rates \eqn{r_1, \ldots,
#' }{r1,\dots, rn}\eqn{ r_n}{r1,\dots, rn} changing at times \eqn{t_1, \ldots,
#' t_n}{t1, \dots, tn}, with \eqn{t_1 = 0}{t1 = 0}. Suppose \eqn{t_k}{tk} is
#' the maximum \eqn{t_i}{ti} such that \eqn{t_i < x}{ti < x}.  The density of
#' this distribution at \eqn{x > 0} is \eqn{f(x)} for \eqn{k = 1}, and
#' \deqn{\prod_{i=1}^k (1 - F(t_{i} - t_{i-1}, r_i)) f(x - t_{k},
#' r_{k})}{\prod{i=1 \dots k} (1 - F(ti - t{i-1}, r{i-1})) f(x - tk, rk)} for k
#' > 1.
#' 
#' where \eqn{F()} and \eqn{f()} are the distribution and density functions of
#' the standard exponential distribution.
#' 
#' If \code{rate} is of length 1, this is just the standard exponential
#' distribution.  Therefore, for example, \code{dpexp(x)}, with no other
#' arguments, is simply equivalent to \code{dexp(x)}.
#' 
#' Only \code{rpexp} is used in the \code{msm} package, to simulate from Markov
#' processes with piecewise-constant intensities depending on time-dependent
#' covariates.  These functions are merely provided for completion, and are not
#' optimized for numerical stability or speed.
#'
#' @name pexp
#' @aliases pexp dpexp ppexp qpexp rpexp
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param rate vector of rates.
#' @param t vector of the same length as \code{rate}, giving the times at which
#' the rate changes. The values of \code{t} should be in increasing order.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p), or
#' log density is returned.
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x],
#' otherwise, P[X > x].
#' @param start numeric scalar; delayed entry time. The random deviates will be
#' left truncated from this start time.
#' @return \code{dpexp} gives the density, \code{ppexp} gives the distribution
#' function, \code{qpexp} gives the quantile function, and \code{rpexp}
#' generates random deviates.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{dexp}}, \code{\link{sim.msm}}.
#' @keywords distribution
#' @examples
#' 
#' x <- seq(0.1, 50, by=0.1)
#' rate <- c(0.1, 0.2, 0.05, 0.3)
#' t <- c(0, 10, 20, 30)
#' ## standard exponential distribution
#' plot(x, dexp(x, 0.1), type="l")
#' ## distribution with piecewise constant rate
#' lines(x, dpexp(x, rate, t), type="l", lty=2)
#' ## standard exponential distribution
#' plot(x, pexp(x, 0.1), type="l")
#' ## distribution with piecewise constant rate
#' lines(x, ppexp(x, rate, t), type="l", lty=2)
#' 
NULL

## The exponential distribution with piecewise-constant rate.  Vector
## of parameters given by rate, change times given by t (first should
## be 0)

#' @rdname pexp
#' @export
dpexp <- function (x, rate = 1, t = 0, log = FALSE)
  {
      if (length(t) != length(rate)) stop("length of t must be equal to length of rate")
      if (!isTRUE(all.equal(0, t[1]))) stop("first element of t should be 0")
      if (is.unsorted(t)) stop("t should be in increasing order")
      ind <- rowSums(outer(x, t, ">="))
      ret <- dexp(x - t[ind], rate[ind], log)
      if (length(t) > 1) {
          dt <- t[-1] - t[-length(t)]
          if (log) {
              cs <- c(0, cumsum(pexp(dt, rate[-length(rate)], log.p=TRUE, lower.tail=FALSE)))
              ret <- cs[ind] + ret
          }
          else {
              cp <- c(1, cumprod(pexp(dt, rate[-length(rate)], lower.tail=FALSE)))
              ret <- cp[ind] * ret
          }
      }
      ret
  }

#' @rdname pexp
#' @export
ppexp <- function(q, rate = 1, t = 0, lower.tail = TRUE, log.p = FALSE)
  {
      if (length(t) != length(rate))
        stop("length of t must be equal to length of rate")
      if (!isTRUE(all.equal(0, t[1]))) stop("first element of t should be 0")
      if (is.unsorted(t)) stop("t should be in increasing order")
      q[q<0] <- 0
      ind <- rowSums(outer(q, t, ">="))
      ret <- pexp(q - t[ind], rate[ind])
      mi <- min(length(t), max(ind))
      if (length(t) > 1) {
          dt <- t[-1] - t[-mi]
          pe <- pexp(dt, rate[-mi])
          cp <- c(1, cumprod(1 - pe))
          ret <- c(0, cumsum(cp[-length(cp)]*pe))[ind] + ret*cp[ind]
      }
      if (!lower.tail) ret <- 1 - ret
      if (log.p) ret <- log(ret)
      ret
  }

#' @rdname pexp
#' @export
qpexp <- function (p, rate = 1, t = 0, lower.tail = TRUE, log.p = FALSE) {
    qgeneric(ppexp, p=p, rate=rate, t=t, special=c("rate", "t"), lower.tail=lower.tail, log.p=log.p)
}

## Simulate n values from exponential distribution with parameters
## rate changing at t.  Simulate from exponentials in turn, simulated
## value is retained if it is less than the next change time.

#' @rdname pexp
#' @export
rpexp <- function(n=1, rate=1, t=0, start=min(t))
{
    if (length(t) != length(rate))
        stop("length of t must be equal to length of rate")
    if (length(start) > 1)
        stop("current implementation does not allow for length(start) > 1")
    if (start < min(t))
        stop("start is less then min(t)")
    if (is.unsorted(t))
        stop("t should be in increasing order")
    if (n == 0)
        return(numeric(0))
    if (length(n) > 1)
        n <- length(n)
    if (length(rate) == 1 || start > max(t))
        return(start + rexp(n, tail(rate,1)))
    if (length(start)==1 && start > min(t)) {
        index <- which(t > start)
        t <- c(start,t[index])
        rate <- rate[c(index[1]-1,index)]
    }
    H <- c(0,cumsum(head(rate,-1)*diff(t)))
    e <- stats::rexp(n)
    i <- findInterval(e,H)
    return(t[i]+(e-H[i])/rate[i])
}

