#' Calculate the probabilities of underlying states and the most likely path
#' through them
#' 
#' For a fitted hidden Markov model, or a model with censored state
#' observations, the Viterbi algorithm recursively constructs the path with the
#' highest probability through the underlying states.  The probability of each
#' hidden state is also computed for hidden Markov models, using the
#' forward-backward algorithm.
#' 
#' 
#' @param x A fitted hidden Markov multi-state model, or a model with censored
#' state observations, as produced by \code{\link{msm}}
#' @param normboot If \code{TRUE}, then before running the algorithm, the
#' maximum likelihood estimates of the model parameters are replaced by an
#' alternative set of parameters drawn randomly from the asymptotic
#' multivariate normal distribution of the MLEs.
#' @param newdata An optional data frame containing observations on which to
#' construct the Viterbi path and forward-backward probabilities. It must be in
#' the same format as the data frame used to fit \code{x}.  If \code{NULL}, the
#' data frame used to fit \code{x} is used.
#' @return A data frame with columns:
#' 
#' \code{subject} = subject identification numbers
#' 
#' \code{time} = times of observations
#' 
#' \code{observed} = corresponding observed states
#' 
#' \code{fitted} = corresponding fitted states found by Viterbi recursion. If
#' the model is not a hidden Markov model, and there are no censored state
#' observations, this is just the observed states.
#' 
#' For hidden Markov models, an additional matrix \code{pstate} is also
#' returned inside the data frame, giving the probability of each hidden state
#' at each point, conditionally on all the data.  This is computed by the
#' forward/backward algorithm.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}
#' @references Durbin, R., Eddy, S., Krogh, A. and Mitchison, G.
#' \emph{Biological sequence analysis}, Cambridge University Press, 1998.
#' @keywords models
#' @export viterbi.msm
viterbi.msm <- function(x, normboot=FALSE, newdata=NULL)
{
  if (!inherits(x, "msm")) stop("expected x to be a msm model")
  
  if (!is.null(newdata)){
    ## initialise msm with new data but do not fit (hence fixedpars = TRUE)
    newcall <- x$call
    newcall$data <- substitute(newdata)
    newcall$fixedpars <- TRUE
    xnew <- try(eval(newcall))
  } else {
    xnew <- x
  }
  xdata <- expand.data(xnew)
  
  if (x$cmodel$ncens > 0 && !x$hmodel$hidden) {
    ## If censoring but not HMM, then define an identity HMM with
    ## true state known at every time except censoring times
    hmod <- vector(x$qmodel$nstates, mode="list")
    for (i in 1:x$qmodel$nstates)
      hmod[[i]] <- hmmIdent(i)
    x$hmodel <- msm.form.hmodel(hmod, est.initprobs=FALSE)
    x$hmodel <- c(x$hmodel, list(ncovs=rep(rep(0,x$hmodel$nstates),x$hmodel$npars), 
                                 ncoveffs=0, nicovs=rep(0,x$hmodel$nstates-1), nicoveffs=0))
    xdata$mf$"(obstrue)" <- ifelse(xdata$mf$"(state)" %in% x$cmodel$censor, 
                                   0, (xdata$mf$"(state)"))
    xdata$mm.hcov <- vector(mode="list", length=x$hmodel$nstates) # reqd by msm.add.hmmcovs
    for (i in seq_len(x$hmodel$nstates))
      xdata$mm.hcov[[i]] <- model.matrix(~1, xdata$mf)
    x$paramdata$allinits <- c(x$paramdata$allinits,x$hmodel$pars)
    x$paramdata$constr <- c(x$paramdata$constr,max(x$paramdata$constr)+seq_along(x$hmodel$pars))
    npts <- attr(xdata$mf, "npts")
    initstate <- xdata$mf$"(state)"[!duplicated(xdata$mf$"(subject)")]
    initp <- matrix(0,nrow=npts,ncol=x$hmodel$nstates)
    for (i in 1:npts){
      if (initstate[i] %in% x$cmodel$censor) {
        cs <- x$cmodel$states_list[[as.character(initstate[i])]]
        initp[i,cs] <- 1/length(cs)
      }
      else initp[i,initstate[i]] <- 1
    }
    x$hmodel$initprobs <- initp
  }

  if (x$hmodel$hidden) {
    if (normboot)
      params <- rmvnorm(1, x$paramdata$opt$par, x$covmat[x$paramdata$optpars,x$paramdata$optpars])
    else
      params <- x$paramdata$opt$par
    
    ret <- Ccall.msm(params,
                     do.what="viterbi",
                     xdata,
                     x$qmodel, x$qcmodel, x$cmodel, x$hmodel, x$paramdata
                     )
    fitted <- ret[[1]]
    pstate <- ret[[2]]
    fitted <- fitted + 1
  } else {
    fitted <- xdata$mf$"(state)"
    pstate <- NULL
  }
  
  if (!is.null(x$qmodel$phase.states)){
    fitted <- x$qmodel$phase.labs[fitted]
  }
  ret <- data.frame(
    subject = xdata$mf$"(subject)",
    time = xdata$mf$"(time)",
    observed = xdata$mf$"(state)",
    fitted = fitted
  )
  if (!is.null(pstate))
    ret$pstate <- pstate
  ret
}
