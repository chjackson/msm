### msm PACKAGE
### USEFUL FUNCTIONS NOT SPECIFIC TO MULTI-STATE MODELS

### Delta method for approximating the covariance matrix of f(X) given cov(X)



#' The delta method
#' 
#' Delta method for approximating the standard error of a transformation
#' \eqn{g(X)} of a random variable \eqn{X = (x_1, x_2, \ldots)}{X = (x1, x2,
#' \ldots)}, given estimates of the mean and covariance matrix of \eqn{X}.
#' 
#' The delta method expands a differentiable function of a random variable
#' about its mean, usually with a first-order Taylor approximation, and then
#' takes the variance. For example, an approximation to the covariance matrix
#' of \eqn{g(X)} is given by
#' 
#' \deqn{ Cov(g(X)) = g'(\mu) Cov(X) [g'(\mu)]^T }{ Cov(g(X)) = g'(mu) Cov(X)
#' [g'(mu)]^T }
#' 
#' where \eqn{\mu}{mu} is an estimate of the mean of \eqn{X}.  This function
#' uses symbolic differentiation via \code{\link{deriv}}.
#' 
#' A limitation of this function is that variables created by the user are not
#' visible within the formula \code{g}.  To work around this, it is necessary
#' to build the formula as a string, using functions such as \code{sprintf},
#' then to convert the string to a formula using \code{as.formula}.  See the
#' example below.
#' 
#' If you can spare the computational time, bootstrapping is a more accurate
#' method of calculating confidence intervals or standard errors for
#' transformations of parameters. See \code{\link{boot.msm}}.  Simulation from
#' the asymptotic distribution of the MLEs (see e.g. Mandel 2013) is also a
#' convenient alternative.
#' 
#' @param g A formula representing the transformation. The variables must be
#' labelled \code{x1, x2,\dots{}} For example,
#' 
#' \code{~ 1 / (x1 + x2)}
#' 
#' If the transformation returns a vector, then a list of formulae representing
#' (\eqn{g_1, g_2, \ldots}{g1, g2, \ldots}) can be provided, for example
#' 
#' \code{list( ~ x1 + x2, ~ x1 / (x1 + x2) )}
#' @param mean The estimated mean of \eqn{X}
#' @param cov The estimated covariance matrix of \eqn{X}
#' @param ses If \code{TRUE}, then the standard errors of \eqn{g_1(X),
#' g_2(X),\ldots}{g1(X), g2(X),\ldots} are returned. Otherwise the covariance
#' matrix of \eqn{g(X)} is returned.
#' @return A vector containing the standard errors of \eqn{g_1(X), g_2(X),
#' \ldots}{g1(X), g2(X), \ldots} or a matrix containing the covariance of
#' \eqn{g(X)}.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @references Oehlert, G. W. (1992) \emph{A note on the delta method}.
#' American Statistician 46(1).
#' 
#' Mandel, M. (2013) \emph{Simulation based confidence intervals for functions
#' with complicated derivatives.} The American Statistician 67(2):76-81.
#' @keywords math
#' @examples
#' 
#' 
#' ## Simple linear regression, E(y) = alpha + beta x 
#' x <- 1:100
#' y <- rnorm(100, 4*x, 5)
#' toy.lm <- lm(y ~ x)
#' estmean <- coef(toy.lm)
#' estvar <- summary(toy.lm)$cov.unscaled * summary(toy.lm)$sigma^2
#' 
#' ## Estimate of (1 / (alphahat + betahat))
#' 1 / (estmean[1] + estmean[2])
#' ## Approximate standard error
#' deltamethod (~ 1 / (x1 + x2), estmean, estvar) 
#' 
#' ## We have a variable z we would like to use within the formula.
#' z <- 1
#' ## deltamethod (~ z / (x1 + x2), estmean, estvar) will not work.
#' ## Instead, build up the formula as a string, and convert to a formula.
#' form <- sprintf("~ %f / (x1 + x2)", z)
#' form
#' deltamethod(as.formula(form), estmean, estvar)
#' 
#' 
#' @export deltamethod
deltamethod <- function(g,       # a formula or list of formulae (functions) giving the transformation g(x) in terms of x1, x2,  etc
                        mean,    # mean, or maximum likelihood estimate, of x
                        cov,     # covariance matrix of x
                        ses=TRUE # return standard errors, else return covariance matrix
                        )
  {
      ## Var (G(x))  =  G'(mu) Var(X) G'(mu)^T
      cov <- as.matrix(cov)
      n <- length(mean)
      if (!is.list(g))
        g <- list(g)
      if ( (dim(cov)[1] != n) || (dim(cov)[2] != n) )
        stop(paste("Covariances should be a ", n, " by ", n, " matrix"))
      syms <- paste("x",1:n,sep="")
      for (i in 1:n)
        assign(syms[i], mean[i])
      gdashmu <- t(sapply(g,
                          function( form ) {
                              as.numeric(attr(eval(
                                                   ## Differentiate each formula in the list
                                                   deriv(form, syms)
                                                   ## evaluate the results at the mean
                                                   ), "gradient"))
                              ## and build the results row by row into a Jacobian matrix
                          }))
      new.covar <- gdashmu %*% cov %*% t(gdashmu)
      if (ses){
          new.se <- sqrt(diag(new.covar))
          new.se
      }
      else
        new.covar
  }

### Matrix exponential
### If a vector of multipliers t is supplied then a list of matrices is returned.



#' Matrix exponential
#' 
#' Calculates the exponential of a square matrix.
#' 
#' See the \code{\link[expm]{expm}} documentation for details of the algorithms
#' it uses.
#' 
#' Generally the exponential \eqn{E} of a square matrix \eqn{M} can often be
#' calculated as
#' 
#' \deqn{E = U \exp(D) U^{-1}}{E = U exp(D) U^{-1}}
#' 
#' where \eqn{D} is a diagonal matrix with the eigenvalues of \eqn{M} on the
#' diagonal, \eqn{\exp(D)}{exp(D)} is a diagonal matrix with the exponentiated
#' eigenvalues of \eqn{M} on the diagonal, and \eqn{U} is a matrix whose
#' columns are the eigenvectors of \eqn{M}.
#' 
#' This method of calculation is used if \code{"pade"} or \code{"series"} is
#' supplied but \eqn{M} has distinct eigenvalues.  I If \eqn{M} has repeated
#' eigenvalues, then its eigenvector matrix may be non-invertible. In this
#' case, the matrix exponential is calculated using the Pade approximation
#' defined by Moler and van Loan (2003), or the less robust power series
#' approximation,
#' 
#' \deqn{\exp(M) = I + M + M^2/2 + M^3 / 3! + M^4 / 4! + ...}{exp(M) = I + M +
#' M^2/2 + M^3 / 3! + M^4 / 4! + ...}
#' 
#' For a continuous-time homogeneous Markov process with transition intensity
#' matrix \eqn{Q}, the probability of occupying state \eqn{s} at time \eqn{u +
#' t} conditional on occupying state \eqn{r} at time \eqn{u} is given by the
#' \eqn{(r,s)} entry of the matrix \eqn{\exp(tQ)}{exp(tQ)}.
#' 
#' If \code{mat} is a valid transition intensity matrix for a continuous-time
#' Markov model (i.e. diagonal entries non-positive, off-diagonal entries
#' non-negative, rows sum to zero), then for certain simpler model structures,
#' there are analytic formulae for the individual entries of the exponential of
#' \code{mat}.  These structures are listed in the PDF manual and the formulae
#' are coded in the \pkg{msm} source file \code{src/analyticp.c}.  These
#' formulae are only used if \code{method="analytic"}.  This is more efficient,
#' but it is not the default in \code{MatrixExp} because the code is not robust
#' to extreme values.  However it is the default when calculating likelihoods
#' for models fitted by \code{\link{msm}}.
#' 
#' The implementation of the Pade approximation used by \code{method="pade"}
#' was taken from JAGS by Martyn Plummer
#' (\url{https://mcmc-jags.sourceforge.io}).
#' 
#' @param mat A square matrix
#' @param t An optional scaling factor for \code{mat}.
#' @param method Under the default of \code{NULL}, this simply wraps the
#' \code{\link[expm]{expm}} function from the \pkg{expm} package.  This is
#' recommended.  Options to \code{\link{expm}} can be supplied to
#' \code{\link{MatrixExp}}, including \code{method}.
#' 
#' Otherwise, for backwards compatibility, the following options, which use
#' code in the \pkg{msm} package, are available: \code{"pade"} for a Pade
#' approximation method, \code{"series"} for the power series approximation, or
#' \code{"analytic"} for the analytic formulae for simpler Markov model
#' intensity matrices (see below).  These options are only used if \code{mat}
#' has repeated eigenvalues, thus the usual eigen-decomposition method cannot
#' be used.
#' @param ... Arguments to pass to \code{\link{expm}}.
#' @return The exponentiated matrix \eqn{\exp(mat)}{exp(mat)}. Or, if \code{t}
#' is a vector of length 2 or more, an array of exponentiated matrices.
#' @references Cox, D. R. and Miller, H. D. \emph{The theory of stochastic
#' processes}, Chapman and Hall, London (1965)
#' 
#' Moler, C and van Loan, C (2003).  Nineteen dubious ways to compute the
#' exponential of a matrix, twenty-five years later.  \emph{SIAM Review}
#' \bold{45}, 3--49.
#' @keywords math
#' @export MatrixExp
MatrixExp <- function(mat, t = 1, method=NULL,...){
    if (!is.matrix(mat) || (nrow(mat)!= ncol(mat)))
        stop("\"mat\" must be a square matrix")
    qmodel <- if (is.qmatrix(mat) && !is.null(method) && method=="analytic") msm.form.qmodel(mat) else list(iso=0, perm=0, qperm=0)
    if (!is.null(method) && method=="analytic") {
        if (!is.qmatrix(mat))
            warning("Analytic method not available since matrix is not a Markov model intensity matrix. Using \"pade\".")
        else if (qmodel$iso==0) warning("Analytic method not available for this Markov model structure. Using \"pade\".")

    }
    if (length(t) > 1) res <- array(dim=c(dim(mat), length(t)))
    for (i in seq_along(t)) {
        if (is.null(method) || !(method %in% c("pade","series","analytic"))) {
            if (is.null(method)) method <- eval(formals(expm::expm)$method)
            resi <- expm::expm(t[i]*mat, method=method, ...)
        } else {
            ccall <- .C("MatrixExpR", as.double(mat), as.integer(nrow(mat)), res=double(length(mat)), as.double(t[i]),
                        as.integer(match(method, c("pade","series"))), # must match macro constants in pijt.c
                        as.integer(qmodel$iso), as.integer(qmodel$perm), as.integer(qmodel$qperm),
                        as.integer(0), NAOK=TRUE)
            resi <- matrix(ccall$res, nrow=nrow(mat))
        }
        if (length(t)==1) res <- resi
        else res[,,i] <- resi
    }
    res
}

## Tests for a valid continuous-time Markov model transition intensity matrix

is.qmatrix <- function(Q) {
    Q2 <- Q; diag(Q2) <- 0
    isTRUE(all.equal(-diag(Q), rowSums(Q2))) && isTRUE(all(diag(Q)<=0)) && isTRUE(all(Q2>=0))
}

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

### Normal distribution with measurement error and optional truncation

#' @rdname medists
#' @export
dmenorm <- function(x, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0, log = FALSE)
{
    sumsq <- sd*sd + sderr*sderr
    sigtmp <- sd*sderr / sqrt(sumsq)
    mutmp <- ((x - meanerr)*sd*sd + mean*sderr*sderr) / sumsq
    nc <- 1/(pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
    nctmp <- pnorm(upper, mutmp, sigtmp) - pnorm(lower, mutmp, sigtmp)
    if (log)
      log(nc) + log(nctmp) + log(dnorm(x, meanerr + mean, sqrt(sumsq), 0))
    else
      nc * nctmp * dnorm(x, meanerr + mean, sqrt(sumsq), 0)
}

#' @rdname medists
#' @export
pmenorm <- function(q, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0, lower.tail = TRUE, log.p = FALSE)
{
    ret <- numeric(length(q))
    dmenorm2 <- function(x)dmenorm(x, mean=mean, sd=sd, lower=lower, upper=upper, sderr=sderr, meanerr=meanerr)
    for (i in 1:length(q)) {
        ret[i] <- integrate(dmenorm2, -Inf, q[i])$value
    }
    if (!lower.tail) ret <- 1 - ret
    if (log.p) ret <- log(ret)
    ret[upper < lower] <- NaN
    ret
}


#' @rdname medists
#' @export
qmenorm <- function(p, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0, lower.tail=TRUE, log.p=FALSE)
{
    qgeneric(pmenorm, p=p, mean=mean, sd=sd, lower=lower, upper=upper, sderr=sderr, meanerr=meanerr,
             lbound=lower, ubound=upper, lower.tail=lower.tail, log.p=log.p)
}

#' @rdname medists
#' @export
rmenorm <- function(n, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0)
{
    rnorm(n, meanerr + rtnorm(n, mean, sd, lower, upper), sderr)
}

### Uniform distribution with measurement error

#' @rdname medists
#' @export
dmeunif <- function(x, lower=0, upper=1, sderr=0, meanerr=0, log = FALSE)
{
    if (log)
      log( pnorm(x, meanerr + lower, sderr) - pnorm(x, meanerr + upper, sderr) ) - log(upper - lower)
    else
      ( pnorm(x, meanerr + lower, sderr) - pnorm(x, meanerr + upper, sderr) ) / (upper - lower)
}

#' @rdname medists
#' @export
pmeunif <- function(q, lower=0, upper=1, sderr=0, meanerr=0, lower.tail = TRUE, log.p = FALSE)
{
    ret <- numeric(length(q))
    dmeunif2 <- function(x)dmeunif(x, lower=lower, upper=upper, sderr=sderr, meanerr=meanerr)
    for (i in 1:length(q)) {
        ret[i] <- integrate(dmeunif2, -Inf, q[i])$value
    }
    if (!lower.tail) ret <- 1 - ret
    if (log.p) ret <- log(ret)
    ret
}

#' @rdname medists
#' @export
qmeunif <- function(p, lower=0, upper=1, sderr=0, meanerr=0, lower.tail=TRUE, log.p=FALSE)
{
    qgeneric(pmeunif, p=p, lower=lower, upper=upper, sderr=sderr, meanerr=meanerr,
             lbound=lower, ubound=upper, lower.tail=lower.tail, log.p=log.p)
}

#' @rdname medists
#' @export
rmeunif <- function(n, lower=0, upper=1, sderr=0, meanerr=0)
{
    rnorm(n, meanerr + runif(n, lower, upper), sderr)
}


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



#' Generic function to find quantiles of a distribution
#' 
#' Generic function to find the quantiles of a distribution, given the
#' equivalent probability distribution function.
#' 
#' This function is intended to enable users to define \code{"q"} functions for
#' new distributions, in cases where the distribution function \code{pdist} is
#' available analytically, but the quantile function is not.
#' 
#' It works by finding the root of the equation \eqn{h(q) = pdist(q) - p = 0}.
#' Starting from the interval \eqn{(-1, 1)}, the interval width is expanded by
#' 50\% until \eqn{h()} is of opposite sign at either end.  The root is then
#' found using \code{\link{uniroot}}.
#' 
#' This assumes a suitably smooth, continuous distribution.
#' 
#' An identical function is provided in the \pkg{flexsurv} package.
#' 
#' @param pdist Probability distribution function, for example,
#' \code{\link{pnorm}} for the normal distribution, which must be defined in
#' the current workspace.  This should accept and return vectorised parameters
#' and values.  It should also return the correct values for the entire real
#' line, for example a positive distribution should have \code{pdist(x)==0} for
#' \eqn{x<0}.
#' @param p Vector of probabilities to find the quantiles for.
#' @param special Vector of character strings naming arguments of the
#' distribution function that should not be vectorised over. Used, for example,
#' for the \code{rate} and \code{t} arguments in \code{\link{qpexp}}.
#' @param ...  The remaining arguments define parameters of the distribution
#' \code{pdist}.  These MUST be named explicitly.
#' 
#' This may also contain the standard arguments \code{log.p} (logical; default
#' \code{FALSE}, if \code{TRUE}, probabilities p are given as log(p)), and
#' \code{lower.tail} (logical; if \code{TRUE} (default), probabilities are P[X
#' <= x] otherwise, P[X > x].).
#' 
#' If the distribution is bounded above or below, then this should contain
#' arguments \code{lbound} and \code{ubound} respectively, and these will be
#' returned if \code{p} is 0 or 1 respectively.  Defaults to \code{-Inf} and
#' \code{Inf} respectively.
#' @return Vector of quantiles of the distribution at \code{p}.
#' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
#' @keywords distribution
#' @examples
#' 
#' qnorm(c(0.025, 0.975), 0, 1)
#' qgeneric(pnorm, c(0.025, 0.975), mean=0, sd=1) # must name the arguments
#' 
#' @export qgeneric
qgeneric <- function(pdist, p, special=NULL, ...)
{
    args <- list(...)
    args.special <- if (!is.null(special)) args[special] else NULL
    args <- args[!names(args) %in% special]
    if (is.null(args$log.p)) args$log.p <- FALSE
    if (is.null(args$lower.tail)) args$lower.tail <- TRUE
    if (is.null(args$lbound)) args$lbound <- -Inf
    if (is.null(args$ubound)) args$ubound <- Inf
    if (args$log.p) p <- exp(p)
    if (!args$lower.tail) p <- 1 - p
    ret <- numeric(length(p))
    ret[p == 0] <- args$lbound
    ret[p == 1] <- args$ubound
    args[c("lower.tail","log.p","lbound","ubound")] <- NULL
    ## Other args assumed to contain params of the distribution.
    ## Replicate all to their maximum length, along with p 
    maxlen <- max(sapply(c(args, p=list(p)), length))
    for (i in seq_along(args))
        args[[i]] <- rep(args[[i]], length.out=maxlen)
    p <- rep(p, length.out=maxlen)

    ret[p < 0 | p > 1] <- NaN
    ind <- (p > 0 & p < 1)
    if (any(ind)) {
        hind <- seq_along(p)[ind]
        h <- function(y) {
            args <- lapply(args, function(x)x[hind[i]])
            p <- p[hind[i]]
            args$q <- y
            args <- c(args, args.special)
            (do.call(pdist, args) - p)
        }
        ptmp <- numeric(length(p[ind]))
        for (i in 1:length(p[ind])) {
            interval <- c(-1, 1)
            while (h(interval[1])*h(interval[2]) >= 0) {
              interval <- interval + c(-1,1)*0.5*(interval[2]-interval[1])
            }
            ptmp[i] <- uniroot(h, interval, tol=.Machine$double.eps)$root
        }
        ret[ind] <- ptmp
    }
    if (any(is.nan(ret))) warning("NaNs produced")
    ret
}


## Transform vector of parameters constrained on [a, b] to real line.
## Vectorised.  a=-Inf or b=Inf represent unbounded below or above.

glogit <- function(x, a, b) {
    if (is.null(a)) a <- -Inf
    if (is.null(b)) b <- Inf
    ret <- numeric(length(x))
    attributes(ret) <- attributes(x)
    nn <- is.infinite(a) & is.infinite(b)
    nb <- is.infinite(a) & is.finite(b)
    an <- is.finite(a) & is.infinite(b)
    ab <- is.finite(a) & is.finite(b)
    ret[nn] <- x[nn]
    ret[nb] <- log(b[nb] - x[nb])
    ret[an] <- log(x[an] - a[an])
    ret[ab] <- log((x[ab] - a[ab]) / (b[ab] - x[ab]))
    ret
}

dglogit <- function(x, a, b) {
    if (is.null(a)) a <- -Inf
    if (is.null(b)) b <- Inf
    ret <- numeric(length(x))
    attributes(ret) <- attributes(x)
    nn <- is.infinite(a) & is.infinite(b)
    nb <- is.infinite(a) & is.finite(b)
    an <- is.finite(a) & is.infinite(b)
    ab <- is.finite(a) & is.finite(b)
    ret[nn] <- 1
    ret[nb] <- -1 / (b[nb] - x[nb])
    ret[an] <- 1 / (x[an] - a[an])
    ret[ab] <- 1/(x[ab] - a[ab]) + 1/(b[ab] - x[ab])
    ret
}

# d/dx log( (x-a)/(b-x) )    = (b-x)/(x-a) * (1/(b-x) + (x-a)/(b-x)^2)
# = 1/(x-a) + 1/(b-x)
# = d/dx   log(x-a) - log(b-x)

## Inverse transform vector of parameters constrained on [a, b]: back
## from real line to constrained scale.  Vectorised.  a=-Inf or b=Inf
## represent unbounded below or above.

gexpit <- function(x, a, b) {
    if (is.null(a)) a <- -Inf
    if (is.null(b)) b <- Inf
    ret <- numeric(length(x))
    attributes(ret) <- attributes(x)
    nn <- is.infinite(a) & is.infinite(b)
    nb <- is.infinite(a) & is.finite(b)
    an <- is.finite(a) & is.infinite(b)
    ab <- is.finite(a) & is.finite(b)
    ret[nn] <- x[nn]
    ret[nb] <- b[nb] - exp(x[nb])
    ret[an] <- exp(x[an]) + a[an]
    ret[ab] <- (b[ab]*exp(x[ab]) + a[ab]) / (1 + exp(x[ab]))
    ret
}

## Derivative of gexpit w.r.t. x

dgexpit <- function(x, a, b) {
    if (is.null(a)) a <- -Inf
    if (is.null(b)) b <- Inf
    ret <- numeric(length(x))
    attributes(ret) <- attributes(x)
    nn <- is.infinite(a) & is.infinite(b)
    nb <- is.infinite(a) & is.finite(b)
    an <- is.finite(a) & is.infinite(b)
    ab <- is.finite(a) & is.finite(b)
    ret[nn] <- 1
    ret[nb] <- - exp(x[nb])
    ret[an] <- exp(x[an])
    ret[ab] <- (b[ab] - a[ab])*exp(x[ab]) / (1 + exp(x[ab]))^2
    ret
}
