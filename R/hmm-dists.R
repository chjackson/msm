
### CONSTRUCTORS FOR VARIOUS DISTRIBUTIONS FOR RESPONSE CONDITIONALLY ON HIDDEN STATE

### Categorical distribution on the set 1,...,n

hmmCat <- function(prob, basecat)
  {
      label <- "categorical"
      prob <- lapply(prob, eval)
      p <- unlist(prob)
      if (any(p < 0)) stop("non-positive probability")
      if (all(p == 0)) stop("insufficient positive probabilities")
      p <- p / sum(p)
      ncats <- length(p)
      link <- "log" # covariates are added to log odds relative to baseline in lik.c(AddCovs)
      cats <- seq(ncats)
      basei <- if (missing(basecat)) which.max(p) else which(cats==basecat)
      r <- function(n, rp=p) sample(cats, size=n, prob=rp, replace=TRUE)
      pars <- c(ncats, basei, p)
      plab <- rep("p", ncats)
      plab[p==0] <- "p0"
      plab[basei] <- "pbase"
      names(pars) <- c("ncats", "basecat", plab)
      hdist <- list(label=label, pars=pars, link=link, r=r) ## probabilities are always pars[c(3,3+pars[0])]
      class(hdist) <- "hmmdist"
      hdist
  }

### Constructor for a standard univariate distribution (i.e. not hmmCat)

hmmDIST <- function(label, link, r, call, ...)
{
    call <- c(as.list(call), list(...))
    miss.pars <- which ( ! (.msm.HMODELPARS[[label]] %in% names(call)[-1]) )
    if (length(miss.pars) > 0) {
        stop("Parameter ", .msm.HMODELPARS[[label]][min(miss.pars)], " for ", call[[1]], " not supplied")
    }
    pars <- unlist(lapply(call[.msm.HMODELPARS[[label]]], eval))
    names(pars) <- .msm.HMODELPARS[[label]]
    hmmCheckInits(pars)
    hdist <- list(label = label,
                  pars = pars,
                  link = link,
                  r = r)
    class(hdist) <- "hmmdist"
    hdist
}

### Multivariate distribution composed of independent univariates

hmmMV <- function(...){
    args <- list(...)
    if (any(sapply(args, class) != "hmmdist")) stop("All arguments of \"hmmMV\" should be HMM distribution objects")
    class(args) <- c("hmmMVdist","hmmdist")
    args
}

hmmCheckInits <- function(pars)
  {
      for (i in names(pars)) {
          if (!is.numeric(pars[i]))
            stop("Expected numeric values for all parameters")
          else if (i %in% .msm.INTEGERPARS) {
              if (!identical(all.equal(pars[i], round(pars[i])), TRUE))
                stop("Value of ", i, " should be integer")
          }
          ## Range check now done in msm.form.hranges
      }
  }

hmmIdent <- function(x)
{
    hmm <- hmmDIST(label = "identity",
                   link = "identity",
                   r = function(n)rep(x, n),
                   match.call())
    hmm$pars <- if (missing(x)) numeric() else x
    names(hmm$pars) <- if(length(hmm$pars)>0) "which" else NULL
    hmm
}

hmmUnif <- function(lower, upper)
  {
      hmmDIST (label = "uniform",
               link = "identity",
               r = function(n) runif(n, lower, upper),
               match.call())
  }

hmmNorm <- function(mean, sd)
  {
      hmmDIST (label = "normal",
               link = "identity",
               r = function(n, rmean=mean) rnorm(n, rmean, sd),
               match.call())
  }

hmmLNorm <- function(meanlog, sdlog)
  {
      hmmDIST (label = "lognormal",
               link = "identity",
               r = function(n, rmeanlog=meanlog) rlnorm(n, rmeanlog, sdlog),
               match.call())
  }

hmmExp <- function(rate)
  {
      hmmDIST (label = "exponential",
               link = "log",
               r = function(n, rrate=rate) rexp(n, rrate),
               match.call())
  }

hmmGamma <- function(shape, rate)
  {
      hmmDIST (label = "gamma",
               link = "log",
               r = function(n, rrate=rate) rgamma(n, shape, rrate),
               match.call())
  }

hmmWeibull <- function(shape, scale)
  {
      hmmDIST (label = "weibull",
               link = "log",
               r = function(n, rscale=scale) rweibull(n, shape, rscale),
               match.call())
  }

hmmPois <- function(rate)
  {
      hmmDIST (label = "poisson",
               link = "log",
               r = function(n, rrate=rate) rpois(n, rrate),
               match.call())
  }

hmmBinom <- function(size, prob)
  {
      hmmDIST (label = "binomial",
               link = "qlogis",
               r = function(n, rprob=prob) rbinom(n, size, rprob),
               match.call())
  }

hmmNBinom <- function(disp, prob)
  {
      hmmDIST (label = "nbinom",
               link = "qlogis",
               r = function(n, rprob=prob) rnbinom(n, disp, rprob),
               match.call())
  }

hmmBeta <- function(shape1, shape2)
  {
      hmmDIST (label = "beta",
               link = "log",
               r = function(n) rbeta(n, shape1, shape2),
               match.call())
  }

hmmTNorm <- function(mean, sd, lower=-Inf, upper=Inf)
  {
      hmmDIST (label = "truncnorm",
               link = "identity",
               r = function(n, rmean=mean) rtnorm(n, rmean, sd, lower, upper),
               match.call(),
               lower=lower,
               upper=upper)
  }

hmmMETNorm <- function(mean, sd, lower, upper, sderr, meanerr=0)
  {
      hmmDIST (label = "metruncnorm",
               link = "identity",
               r = function(n, rmeanerr=meanerr) rnorm(n, rmeanerr + rtnorm(n, mean, sd, lower, upper), sderr),
               match.call(),
               meanerr=meanerr)
  }

hmmMEUnif <- function(lower, upper, sderr, meanerr=0)
  {
      hmmDIST (label = "meuniform",
               link = "identity",
               r = function(n, rmeanerr=meanerr) rnorm(n, rmeanerr + runif(n, lower, upper), sderr),
               match.call(),
               meanerr=meanerr)
  }

hmmT <- function(mean, scale, df)
  {
      hmmDIST(label="t",
              link="identity",
              r = function(n, rmean=mean) { rmean + scale*rt(n,df) },
              match.call())
  }
