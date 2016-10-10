### msm PACKAGE
### USEFUL FUNCTIONS NOT SPECIFIC TO MULTI-STATE MODELS

### Delta method for approximating the covariance matrix of f(X) given cov(X)

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
    for (i in seq(along=t)) {
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

qtnorm <- function(p, mean=0, sd=1, lower=-Inf, upper=Inf, lower.tail=TRUE, log.p=FALSE)
{
    qgeneric(ptnorm, p=p, mean=mean, sd=sd, lower=lower, upper=upper, lbound=lower, ubound=upper, lower.tail=lower.tail, log.p=log.p)
}

## Rejection sampling algorithm by Robert (Stat. Comp (1995), 5, 121-5)
## for simulating from the truncated normal distribution.

rtnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
    if (length(n) > 1)
        n <- length(n)
    mean <- rep(mean, length=n)
    sd <- rep(sd, length=n)
    lower <- rep(lower, length=n)
    upper <- rep(upper, length=n)
    lower <- (lower - mean) / sd ## Algorithm works on mean 0, sd 1 scale
    upper <- (upper - mean) / sd
    ind <- seq(length=n)
    ret <- numeric(n)
    ## Different algorithms depending on where upper/lower limits lie.
    alg <- ifelse(
                  lower > upper,
                  -1,# return NaN if lower > upper
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
                                       3)))) # rejection sampling with uniform proposal. Use if bounds are narrow and central.

    ind.nan <- ind[alg==-1]; ind.no <- ind[alg==0]; ind.expl <- ind[alg==1]; ind.expu <- ind[alg==2]; ind.u <- ind[alg==3]
    ret[ind.nan] <- NaN
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


qmenorm <- function(p, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0, lower.tail=TRUE, log.p=FALSE)
{
    qgeneric(pmenorm, p=p, mean=mean, sd=sd, lower=lower, upper=upper, sderr=sderr, meanerr=meanerr,
             lbound=lower, ubound=upper, lower.tail=lower.tail, log.p=log.p)
}

rmenorm <- function(n, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0)
{
    rnorm(n, meanerr + rtnorm(n, mean, sd, lower, upper), sderr)
}

### Uniform distribution with measurement error

dmeunif <- function(x, lower=0, upper=1, sderr=0, meanerr=0, log = FALSE)
{
    if (log)
      log( pnorm(x, meanerr + lower, sderr) - pnorm(x, meanerr + upper, sderr) ) - log(upper - lower)
    else
      ( pnorm(x, meanerr + lower, sderr) - pnorm(x, meanerr + upper, sderr) ) / (upper - lower)
}

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

qmeunif <- function(p, lower=0, upper=1, sderr=0, meanerr=0, lower.tail=TRUE, log.p=FALSE)
{
    qgeneric(pmeunif, p=p, lower=lower, upper=upper, sderr=sderr, meanerr=meanerr,
             lbound=lower, ubound=upper, lower.tail=lower.tail, log.p=log.p)
}

rmeunif <- function(n, lower=0, upper=1, sderr=0, meanerr=0)
{
    rnorm(n, meanerr + runif(n, lower, upper), sderr)
}


## The exponential distribution with piecewise-constant rate.  Vector
## of parameters given by rate, change times given by t (first should
## be 0)

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

qpexp <- function (p, rate = 1, t = 0, lower.tail = TRUE, log.p = FALSE) {
    qgeneric(ppexp, p=p, rate=rate, t=t, lower.tail=lower.tail, log.p=log.p)
}

## Simulate n values from exponential distribution with parameters
## rate changing at t.  Simulate from exponentials in turn, simulated
## value is retained if it is less than the next change time.

rpexp <- function(n=1, rate=1, t=0)
  {
      if (length(t) != length(rate)) stop("length of t must be equal to length of rate")
      if (!isTRUE(all.equal(0, t[1]))) stop("first element of t should be 0")
      if (is.unsorted(t)) stop("t should be in increasing order")
      if (length(rate) == 1) return(rexp(n, rate))
      if (n == 0) return(numeric(0))
      if (length(n) > 1) n <- length(n)
      ret <- numeric(n)                          # outcome is a vector length n
      left <- 1:n
      for (i in seq(along=rate)){
          re <- rexp(length(left), rate[i])   # simulate as many exponentials as there are values remaining
          r <- t[i] + re
          success <- if (i == length(rate)) seq(along=left) else which(r < t[i+1])
          ret[left[success]] <- r[success]
          left <- setdiff(left, left[success])  # indices of values in outcome remaining to simulate.
          if (length(left)==0) break;
      }
      ret
  }

qgeneric <- function(pdist, p, ...)
{
    args <- list(...)
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
    for (i in seq(along=args))
        args[[i]] <- rep(args[[i]], length.out=maxlen)
    p <- rep(p, length.out=maxlen)

    ret[p < 0 | p > 1] <- NaN
    ind <- (p > 0 & p < 1)
    if (any(ind)) {
        hind <- seq(along=p)[ind]
        h <- function(y) {
            args <- lapply(args, function(x)x[hind[i]])
            p <- p[hind[i]]
            args$q <- y
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
