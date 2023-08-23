### CONSTRUCTORS FOR VARIOUS DISTRIBUTIONS FOR RESPONSE CONDITIONALLY ON HIDDEN STATE

### Categorical distribution on the set 1,...,n

#' @rdname hmm-dists
#' @export
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



#' Multivariate hidden Markov models
#' 
#' Constructor for a a multivariate hidden Markov model (HMM) where each of the
#' \code{n} variables observed at the same time has a (potentially different)
#' standard univariate distribution conditionally on the underlying state.  The
#' \code{n} outcomes are independent conditionally on the hidden state.
#' 
#' If a particular state in a HMM has such an outcome distribution, then a call
#' to \code{\link{hmmMV}} is supplied as the corresponding element of the
#' \code{hmodel} argument to \code{\link{msm}}.  See Example 2 below.
#' 
#' A multivariate HMM where multiple outcomes at the same time are generated
#' from the \emph{same} distribution is specified in the same way as the
#' corresponding univariate model, so that \code{\link{hmmMV}} is not required.
#' The outcome data are simply supplied as a matrix instead of a vector.  See
#' Example 1 below.
#' 
#' The outcome data for such models are supplied as a matrix, with number of
#' columns equal to the maximum number of arguments supplied to the
#' \code{\link{hmmMV}} calls for each state.  If some but not all of the
#' variables are missing (\code{NA}) at a particular time, then the observed
#' data at that time still contribute to the likelihood.  The missing data are
#' assumed to be missing at random.  The Viterbi algorithm may be used to
#' predict the missing values given the fitted model and the observed data.
#' 
#' Typically the outcome model for each state will be from the same family or
#' set of families, but with different parameters.  Theoretically, different
#' numbers of distributions may be supplied for different states.  If a
#' particular state has fewer outcomes than the maximum, then the data for that
#' state are taken from the first columns of the response data matrix.  However
#' this is not likely to be a useful model, since the number of observations
#' will probably give information about the underlying state, violating the
#' missing at random assumption.
#' 
#' Models with outcomes that are dependent conditionally on the hidden state
#' (e.g. correlated multivariate normal observations) are not currently
#' supported.
#' 
#' 
#' @param ...  The number of arguments supplied should equal the maximum number
#' of observations made at one time.  Each argument represents the univariate
#' distribution of that outcome conditionally on the hidden state, and should
#' be the result of calling a univariate hidden Markov model constructor (see
#' \code{\link{hmm-dists}}).
#' @return A list of objects, each of class \code{hmmdist} as returned by the
#' univariate HMM constructors documented in \code{\link{hmm-dists}}.  The
#' whole list has class \code{hmmMVdist}, which inherits from \code{hmmdist}.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{hmm-dists}},\code{\link{msm}}
#' @references Jackson, C. H., Su, L., Gladman, D. D. and Farewell, V. T.
#' (2015) On modelling minimal disease activity.  Arthritis Care and Research
#' (early view).
#' @keywords distribution
#' @examples
#' 
#' ## Simulate data from a Markov model 
#' nsubj <- 30; nobspt <- 5
#' sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt),
#'                      time = seq(0, 20, length=nobspt))
#' set.seed(1)
#' two.q <- rbind(c(-0.1, 0.1), c(0, 0))
#' dat <- simmulti.msm(sim.df[,1:2], qmatrix=two.q, drop.absorb=FALSE)
#' 
#' ### EXAMPLE 1
#' ## Generate two observations at each time from the same outcome
#' ## distribution:
#' ## Bin(40, 0.1) for state 1, Bin(40, 0.5) for state 2
#' dat$obs1[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.1)
#' dat$obs2[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.1)
#' dat$obs1[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.5)
#' dat$obs2[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.5)
#' dat$obs <- cbind(obs1 = dat$obs1, obs2 = dat$obs2)
#' 
#' ## Fitted model should approximately recover true parameters 
#' msm(obs ~ time, subject=subject, data=dat, qmatrix=two.q,
#'     hmodel = list(hmmBinom(size=40, prob=0.2),
#'                   hmmBinom(size=40, prob=0.2)))
#' 
#' ### EXAMPLE 2
#' ## Generate two observations at each time from different
#' ## outcome distributions:
#' ## Bin(40, 0.1) and Bin(40, 0.2) for state 1, 
#' dat$obs1 <- dat$obs2 <- NA
#' dat$obs1[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.1)
#' dat$obs2[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.2)
#' 
#' ## Bin(40, 0.5) and Bin(40, 0.6) for state 2
#' dat$obs1[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.6)
#' dat$obs2[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.5)
#' dat$obs <- cbind(obs1 = dat$obs1, obs2 = dat$obs2)
#' 
#' ## Fitted model should approximately recover true parameters 
#' msm(obs ~ time, subject=subject, data=dat, qmatrix=two.q,   
#'     hmodel = list(hmmMV(hmmBinom(size=40, prob=0.3),
#'                         hmmBinom(size=40, prob=0.3)),                 
#'                  hmmMV(hmmBinom(size=40, prob=0.3),
#'                        hmmBinom(size=40, prob=0.3))),
#'     control=list(maxit=10000))
#' 
#' @export hmmMV
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

#' @rdname hmm-dists
#' @export
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

#' @rdname hmm-dists
#' @export
hmmUnif <- function(lower, upper)
  {
      hmmDIST (label = "uniform",
               link = "identity",
               r = function(n) runif(n, lower, upper),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmNorm <- function(mean, sd)
  {
      hmmDIST (label = "normal",
               link = "identity",
               r = function(n, rmean=mean) rnorm(n, rmean, sd),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmLNorm <- function(meanlog, sdlog)
  {
      hmmDIST (label = "lognormal",
               link = "identity",
               r = function(n, rmeanlog=meanlog) rlnorm(n, rmeanlog, sdlog),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmExp <- function(rate)
  {
      hmmDIST (label = "exponential",
               link = "log",
               r = function(n, rrate=rate) rexp(n, rrate),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmGamma <- function(shape, rate)
  {
      hmmDIST (label = "gamma",
               link = "log",
               r = function(n, rrate=rate) rgamma(n, shape, rrate),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmWeibull <- function(shape, scale)
  {
      hmmDIST (label = "weibull",
               link = "log",
               r = function(n, rscale=scale) rweibull(n, shape, rscale),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmPois <- function(rate)
  {
      hmmDIST (label = "poisson",
               link = "log",
               r = function(n, rrate=rate) rpois(n, rrate),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmBinom <- function(size, prob)
  {
      hmmDIST (label = "binomial",
               link = "qlogis",
               r = function(n, rprob=prob) rbinom(n, size, rprob),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmBetaBinom <- function(size, meanp, sdp)
  {
      hmmDIST (label = "betabinomial",
               link = "qlogis",
               r = function(n) rbinom(n, size, rbeta(n, meanp/sdp, (1-meanp)/sdp)),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmNBinom <- function(disp, prob)
  {
      hmmDIST (label = "nbinom",
               link = "qlogis",
               r = function(n, rprob=prob) rnbinom(n, disp, rprob),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmBeta <- function(shape1, shape2)
  {
      hmmDIST (label = "beta",
               link = "log",
               r = function(n) rbeta(n, shape1, shape2),
               match.call())
  }

#' @rdname hmm-dists
#' @export
hmmTNorm <- function(mean, sd, lower=-Inf, upper=Inf)
  {
      hmmDIST (label = "truncnorm",
               link = "identity",
               r = function(n, rmean=mean) rtnorm(n, rmean, sd, lower, upper),
               match.call(),
               lower=lower,
               upper=upper)
  }

#' @rdname hmm-dists
#' @export
hmmMETNorm <- function(mean, sd, lower, upper, sderr, meanerr=0)
  {
      hmmDIST (label = "metruncnorm",
               link = "identity",
               r = function(n, rmeanerr=meanerr) rnorm(n, rmeanerr + rtnorm(n, mean, sd, lower, upper), sderr),
               match.call(),
               meanerr=meanerr)
  }

#' @rdname hmm-dists
#' @export
hmmMEUnif <- function(lower, upper, sderr, meanerr=0)
  {
      hmmDIST (label = "meuniform",
               link = "identity",
               r = function(n, rmeanerr=meanerr) rnorm(n, rmeanerr + runif(n, lower, upper), sderr),
               match.call(),
               meanerr=meanerr)
  }

#' @rdname hmm-dists
#' @export
hmmT <- function(mean, scale, df)
  {
      hmmDIST(label="t",
              link="identity",
              r = function(n, rmean=mean) { rmean + scale*rt(n,df) },
              match.call())
  }
