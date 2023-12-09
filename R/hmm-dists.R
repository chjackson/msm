
#' Hidden Markov model constructors
#' 
#' These functions are used to specify the distribution of the response
#' conditionally on the underlying state in a hidden Markov model.  A list of
#' these function calls, with one component for each state, should be used for
#' the \code{hmodel} argument to \code{msm}. The initial values for the
#' parameters of the distribution should be given as arguments. Note the
#' initial values should be supplied as literal values - supplying them as
#' variables is currently not supported.
#' 
#' \code{hmmCat} represents a categorical response distribution on the set
#' \code{1, 2, \dots{}, length(prob)}.  The Markov model with misclassification
#' is an example of this type of model. The categories in this case are (some
#' subset of) the underlying states.
#' 
#' The \code{hmmIdent} distribution is used for underlying states which are
#' observed exactly without error.  For hidden Markov models with multiple
#' outcomes, (see \code{\link{hmmMV}}), the outcome in the data which takes the
#' special \code{hmmIdent} value must be the first of the multiple outcomes.
#' 
#' \code{hmmUnif}, \code{hmmNorm}, \code{hmmLNorm}, \code{hmmExp},
#' \code{hmmGamma}, \code{hmmWeibull}, \code{hmmPois}, \code{hmmBinom},
#' \code{hmmTNorm}, \code{hmmNBinom} and \code{hmmBeta} represent Uniform,
#' Normal, log-Normal, exponential, Gamma, Weibull, Poisson, Binomial,
#' truncated Normal, negative binomial and beta distributions, respectively,
#' with parameterisations the same as the default parameterisations in the
#' corresponding base R distribution functions.
#' 
#' \code{hmmT} is the Student t distribution with general mean \eqn{\mu}{mu},
#' scale \eqn{\sigma}{sigma} and degrees of freedom \code{df}.  The variance is
#' \eqn{\sigma^2 df/(df + 2)}{sigma^2 df/(df + 2)}.  Note the t distribution in
#' base R \code{\link{dt}} is a standardised one with mean 0 and scale 1.
#' These allow any positive (integer or non-integer) \code{df}.  By default,
#' all three parameters, including \code{df}, are estimated when fitting a
#' hidden Markov model, but in practice, \code{df} might need to be fixed for
#' identifiability - this can be done using the \code{fixedpars} argument to
#' \code{\link{msm}}.
#' 
#' The \code{hmmMETNorm} and \code{hmmMEUnif} distributions are truncated
#' Normal and Uniform distributions, but with additional Normal measurement
#' error on the response. These are generalisations of the distributions
#' proposed by Satten and Longini (1996) for modelling the progression of CD4
#' cell counts in monitoring HIV disease.  See \code{\link{medists}} for
#' density, distribution, quantile and random generation functions for these
#' distributions.  See also \code{\link{tnorm}} for density, distribution,
#' quantile and random generation functions for the truncated Normal
#' distribution.
#' 
#' See the PDF manual \file{msm-manual.pdf} in the \file{doc} subdirectory for
#' algebraic definitions of all these distributions.  New hidden Markov model
#' response distributions can be added to \pkg{msm} by following the
#' instructions in Section 2.17.1.
#' 
#' Parameters which can be modelled in terms of covariates, on the scale of a
#' link function, are as follows.
#' 
#' \tabular{ll}{ PARAMETER NAME \tab LINK FUNCTION \cr \code{mean} \tab
#' identity \cr \code{meanlog} \tab identity \cr \code{rate} \tab log \cr
#' \code{scale} \tab log \cr \code{meanerr} \tab identity \cr \code{meanp} \tab
#' logit \cr \code{prob} \tab logit or multinomial logit }
#' 
#' Parameters \code{basecat, lower, upper, size, meanerr} are fixed at their
#' initial values. All other parameters are estimated while fitting the hidden
#' Markov model, unless the appropriate \code{fixedpars} argument is supplied
#' to \code{msm}.
#' 
#' For categorical response distributions \code{(hmmCat)} the outcome
#' probabilities initialized to zero are fixed at zero, and the probability
#' corresponding to \code{basecat} is fixed to one minus the sum of the
#' remaining probabilities.  These remaining probabilities are estimated, and
#' can be modelled in terms of covariates via multinomial logistic regression
#' (relative to \code{basecat}).
#'
#' @name hmm-dists
#' @aliases hmm-dists hmmCat hmmIdent hmmUnif hmmNorm hmmLNorm hmmExp hmmGamma
#' hmmWeibull hmmPois hmmBinom hmmTNorm hmmMETNorm hmmMEUnif hmmNBinom
#' hmmBetaBinom hmmBeta hmmT
#' @param prob (\code{hmmCat}) Vector of probabilities of observing category
#' \code{1, 2, \dots{}, length(prob)} respectively.  Or the probability
#' governing a binomial or negative binomial distribution.
#' @param basecat (\code{hmmCat}) Category which is considered to be the
#' "baseline", so that during estimation, the probabilities are parameterised
#' as probabilities relative to this baseline category. By default, the
#' category with the greatest probability is used as the baseline.
#' @param x (\code{hmmIdent}) Code in the data which denotes the
#' exactly-observed state.
#' @param mean (\code{hmmNorm,hmmLNorm,hmmTNorm}) Mean defining a Normal, or
#' truncated Normal distribution.
#' @param sd (\code{hmmNorm,hmmLNorm,hmmTNorm}) Standard deviation defining a
#' Normal, or truncated Normal distribution.
#' @param meanlog (\code{hmmNorm,hmmLNorm,hmmTNorm}) Mean on the log scale, for
#' a log Normal distribution.
#' @param sdlog (\code{hmmNorm,hmmLNorm,hmmTNorm}) Standard deviation on the
#' log scale, for a log Normal distribution.
#' @param rate (\code{hmmPois,hmmExp,hmmGamma}) Rate of a Poisson, Exponential
#' or Gamma distribution (see \code{\link{dpois}}, \code{\link{dexp}},
#' \code{\link{dgamma}}).
#' @param shape (\code{hmmPois,hmmExp,hmmGamma}) Shape parameter of a Gamma or
#' Weibull distribution (see \code{\link{dgamma}}, \code{\link{dweibull}}).
#' @param shape1,shape2 First and second parameters of a beta distribution (see
#' \code{\link{dbeta}}).
#' @param scale (\code{hmmGamma}) Scale parameter of a Gamma distribution (see
#' \code{\link{dgamma}}), or unstandardised Student t distribution.
#' @param df Degrees of freedom of the Student t distribution.
#' @param size Order of a Binomial distribution (see \code{\link{dbinom}}).
#' @param disp Dispersion parameter of a negative binomial distribution, also
#' called \code{size} or \code{order}.  (see \code{\link{dnbinom}}).
#' @param meanp Mean outcome probability in a beta-binomial distribution
#' @param sdp Standard deviation describing the overdispersion of the outcome
#' probability in a beta-binomial distribution
#' @param lower (\code{hmmUnif,hmmTNorm,hmmMEUnif}) Lower limit for an Uniform
#' or truncated Normal distribution.
#' @param upper (\code{hmmUnif,hmmTNorm,hmmMEUnif}) Upper limit for an Uniform
#' or truncated Normal distribution.
#' @param sderr (\code{hmmMETNorm,hmmUnif}) Standard deviation of the Normal
#' measurement error distribution.
#' @param meanerr (\code{hmmMETNorm,hmmUnif}) Additional shift in the
#' measurement error, fixed to 0 by default.  This may be modelled in terms of
#' covariates.
#' @return Each function returns an object of class \code{hmmdist}, which is a
#' list containing information about the model.  The only component which may
#' be useful to end users is \code{r}, a function of one argument \code{n}
#' which returns a random sample of size \code{n} from the given distribution.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}
#' @references Satten, G.A. and Longini, I.M.  Markov chains with measurement
#' error: estimating the 'true' course of a marker of the progression of human
#' immunodeficiency virus disease (with discussion) \emph{Applied Statistics}
#' 45(3): 275-309 (1996).
#' 
#' Jackson, C.H. and Sharples, L.D. Hidden Markov models for the onset and
#' progresison of bronchiolitis obliterans syndrome in lung transplant
#' recipients \emph{Statistics in Medicine}, 21(1): 113--128 (2002).
#' 
#' Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W. and Couto, E.
#' Multi-state Markov models for disease progression with classification error.
#' \emph{The Statistician}, 52(2): 193--209 (2003).
#' @keywords distribution
NULL


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
