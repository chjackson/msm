deriv.supported <- function(msmdata, hmodel, cmodel){
    (!hmodel$hidden || (hmodel$hidden &&
                        ## Models where derivatives not supported: 
                        ## multiple outcomes 
#                        (is.null(ncol(msmdata$mf$"(state)")) || (ncol(msmdata$mf$"(state)")==1)) &&
                        ## unknown initial state probs
                        !hmodel$est.initprobs &&
                        ## constraints on misclassification / categorical outcome probabilities
                        (!any(duplicated(hmodel$constr[hmodel$plabs=="p"]))) &&
                        (!any(duplicated(hmodel$covconstr[.msm.HMODELS[hmodel$models[hmodel$coveffstate]]=="categorical"]))) &&
                        all(.msm.HMODELS[hmodel$models %in% .msm.HMODELS.DERIV]) 
                        ))
}

info.supported <- function(msmdata, hmodel, cmodel){
    ## only panel data and either non-hidden or misclassification models
    deriv.supported(msmdata, hmodel, cmodel) &&
        all(msmdata$mf$"(obstype)"==1) &&
            (!hmodel$hidden || (hmodel$hidden &&
                                all(.msm.HMODELS[hmodel$models] %in% .msm.HMODELS.INFO)))
}

### GENERIC OPTIMISATION FUNCTION

msm.optim <- function(opt.method, p, hessian, use.deriv, msmdata, qmodel, qcmodel, cmodel, hmodel, ...){
    p$params <- p$allinits
    gr <- if (deriv.supported(msmdata,hmodel,cmodel) && use.deriv) deriv.msm else NULL
    optfn <- paste("msm.optim", opt.method, sep=".")
    if (!exists(optfn)) stop("Unknown optimisation method \"", opt.method, "\"")

    assign("nliks", 0, envir=msm.globals)
    args <- c(list(p=p, gr=gr, hessian=hessian,
                   msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                   cmodel=cmodel, hmodel=hmodel), list(...))
    p <- do.call(optfn, args)
    assign("nliks", 0, envir=msm.globals)

    if (isTRUE(getOption("msm.test.analytic.derivatives"))){
        if (!deriv.supported(msmdata,hmodel,cmodel)) warning("Analytic derivatives not available for this model")
        else p$deriv.test <- deriv.test(msmdata, qmodel, qcmodel, cmodel, hmodel, msm.unfixallparams(p))
    }
    ## Attach derivative and information matrix at the MLE.
    ## If all parameters fixed, do this at the initial values for all.
    pp <- if (opt.method=="fixed") msm.unfixallparams(p) else p
    pi <- pp$params[pp$optpars]
    if (deriv.supported(msmdata,hmodel,cmodel)) {
        p$deriv <- deriv.msm(pi, msmdata, qmodel, qcmodel, cmodel, hmodel, pp)
        if (info.supported(msmdata,hmodel,cmodel)) {
            p$information <- information.msm(pi, msmdata, qmodel, qcmodel, cmodel, hmodel, pp)
        }
    }
    p
}


## SPECIFIC OPTIMISATION METHODS
## All should take same arguments
## All should add $lik, $params, $opt terms to "paramdata" object "p"
### This allows the user to just drop in a new method without touching the rest of the msm package code


## Trivial one where all the parameters are fixed at their initial values

msm.optim.fixed <- function(p, gr, hessian, msmdata, qmodel, qcmodel, cmodel, hmodel, ...) {
    p$lik <- lik.msm(p$inits, msmdata, qmodel, qcmodel, cmodel, hmodel, p)
    p$opt <- list(par = p$allinits[!duplicated(abs(p$constr))])
    p$params <- msm.rep.constraints(p$params, p, hmodel)

    p$allinits[!duplicated(abs(p$constr))][abs(p$constr)]*sign(p$constr)
    p.unfix <- msm.unfixallparams(p)
    if (hessian)
        p$opt$hessian <- optimHess(par=p.unfix$inits, fn=lik.msm, gr=gr,
                                   msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                   cmodel=cmodel, hmodel=hmodel, paramdata=p.unfix)
    p
}

msm.optim.optim <- function(p, gr, hessian, msmdata, qmodel, qcmodel, cmodel, hmodel, ...) {
    optim.args <- list(...)
    if (is.null(optim.args$method))
        optim.args$method <- if (deriv.supported(msmdata, hmodel, cmodel) || (length(p$inits)==1)) "BFGS" else "Nelder-Mead"

    if (is.null(optim.args$control)) optim.args$control <- list()

# this might cause more trouble than it solves.    

#    	* R/optim.R: Pass fnscale to optim() automatically, if not already
#	provided, using the likelihood at the initial values.  Wastes a
#	likelihood calculation but should improve convergence.

#    if (is.null(optim.args$control$fnscale))
#        optim.args$control$fnscale <- lik.msm(p$inits, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
#                                              cmodel=cmodel, hmodel=hmodel, paramdata=p)
    
    optim.args <- c(optim.args, list(par=p$inits, fn=lik.msm, hessian=hessian, gr=gr,
                                     msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                     cmodel=cmodel, hmodel=hmodel, paramdata=p))
    opt <- do.call("optim", optim.args)
    if (opt$convergence==1)
        warning("Iteration limit in optim() reached without convergence. Reported estimates are not the maximum likelihood. Increase \"maxit\" or change optimisation method - see help(optim) and help(msm).")
    else if (opt$convergence==10)
        warning("Not converged: Nelder-Mead simplex is degenerate. Reported estimates are not the maximum likelihood.")
    else if (opt$convergence %in% 51:52)
        warning("Not converged: error in L-BFGS-B, see help(optim). Reported estimates are not the maximum likelihood.")
    ctrace <- !is.null(optim.args$control$trace) && optim.args$control$trace > 0
    if (ctrace){
        cat("Used", opt$counts[1], "function and", opt$counts[2], "gradient evaluations\n")
    }
    if (!is.null(opt$message)) warning("optim() returned a message: ",opt$message)
    p$lik <- opt$value
    p$params[p$optpars] <- opt$par
    p$opt <- opt
    p
}

msm.optim.nlm <- function(p, gr, hessian, msmdata, qmodel, qcmodel, cmodel, hmodel, ...) {
    nlmfn <- function(par) {
        ret <- lik.msm(par, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                       cmodel=cmodel, hmodel=hmodel, paramdata=p)
        if (!is.null(gr))
            attr(ret, "gradient") <- deriv.msm(par, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                               cmodel=cmodel, hmodel=hmodel, paramdata=p)
        ret
    }
    optim.args <- c(list(...), list(f=nlmfn, p=p$inits, hessian=hessian))
    ## suspect this check is excessively precise
    if (!is.null(optim.args$check.analyticals)) warning("Forcing check.analyticals=FALSE in \"nlm\"")
    optim.args$check.analyticals <- FALSE
    opt <- do.call("nlm", optim.args)
    if (opt$code==4)
        warning("Iteration limit in nlm() reached without convergence. Reported estimates are not the maximum likelihood. Increase \"iterlim\" or change optimisation method - see help(nlm) and help(msm).")
    else if (opt$code==3)
        warning("Not converged: local minimum or step size too small.  Reported estimates are not the maximum likelihood.  See help(nlm)")
    else if (opt$code==5)
        warning("Not converged: maximum step size exceeded five consecutive times.  Reported estimates are not the maximum likelihood.  See help(nlm)")
    ctrace <- !is.null(list(...)$print.level) && list(...)$print.level > 0
    if (ctrace) cat("Used", opt$iterations, "iterations\n")
    p$lik <- opt$minimum
    p$params[p$optpars] <- opt$estimate
    p$opt <- opt
    p
}

msm.optim.fisher <- function(p, gr, hessian, msmdata, qmodel, qcmodel, cmodel, hmodel, ...) {
    if (hmodel$hidden)
        stop("Fisher scoring not supported for hidden Markov models or censored states")
    if (cmodel$ncens > 0)
        stop("Fisher scoring not supported with censored states")
    if (any(msmdata$mf$"(obstype)"==2))
        stop("Fisher scoring not supported with exact transition times")
    if (any(msmdata$mf$"(obstype)"==3))
        stop("Fisher scoring not supported with exact death times")
    optim.args <- list(...)
    if (is.null(optim.args$control$reltol)) reltol <- sqrt(.Machine$double.eps)
    damp <- if (is.null(optim.args$control$damp)) 0 else optim.args$control$damp
    theta <- p$inits
    lik.old <- -lik.msm(theta, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                        cmodel=cmodel, hmodel=hmodel, paramdata=p)
    converged <- FALSE
    iterations <- 1
    ctrace <- !is.null(optim.args$control$trace) && optim.args$control$trace > 0
    while(!converged) {
        if (ctrace) cat("-2loglik=",-lik.old,", pars=",theta,"\n")
        Info <- information.msm(theta, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                              cmodel=cmodel, hmodel=hmodel, paramdata=p)
        V <- -deriv.msm(theta, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                              cmodel=cmodel, hmodel=hmodel, paramdata=p)
        Info <- Info + diag(damp, nrow(Info))
        theta <- theta + solve(Info, V)
        lik.new <- -lik.msm(theta, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                            cmodel=cmodel, hmodel=hmodel, paramdata=p)
        iterations <- iterations + 1
        if (abs(lik.old - lik.new) < reltol*(abs(lik.old) + reltol))
            converged <- TRUE
        else lik.old <- lik.new
    }
    if (ctrace) cat("Used", iterations, "evaluations of likelihood and information\n")
    p$lik <- -lik.new
    p$params[p$optpars] <- theta
    opt <- list(minimum=-lik.new, estimate=theta, value=-lik.new, par=theta, iterations=iterations)
    if (hessian)
        opt$hessian <- optimHess(par=theta, fn=lik.msm, gr=gr,
                                 msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                 cmodel=cmodel, hmodel=hmodel, paramdata=p)
    p$opt <- opt
    p
}

msm.optim.bobyqa <- function(p, gr, hessian, msmdata, qmodel, qcmodel, cmodel, hmodel, ...) {
    optim.args <- list(...)
    optim.args <- c(optim.args, list(par=p$inits, fn=lik.msm,
                                     msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                     cmodel=cmodel, hmodel=hmodel, paramdata=p))
    if (requireNamespace("minqa", quietly = TRUE)){
        opt <- do.call("bobyqa", optim.args)
    } else stop("\"minqa\" package not available")
    if (opt$ierr %in% 1:5)
        warning(opt$msg, ". Reported estimates are not the maximum likelihood. See help(bobyqa) and help(msm).")
    if (hessian)
        opt$hessian <- optimHess(par=opt$par, fn=lik.msm, gr=gr,
                                 msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                 cmodel=cmodel, hmodel=hmodel, paramdata=p)
    ctrace <- !is.null(optim.args$control$iprint) && optim.args$control$iprint > 0
    if (ctrace) cat("Used", opt$feval, "function evaluations\n")
    p$lik <- opt$fval
    p$params[p$optpars] <- opt$par
    p$opt <- opt
    p
}




### Test analytic against numeric derivatives

deriv.test <- function(msmdata, qmodel, qcmodel, cmodel, hmodel, p){
    an.d <- deriv.msm(p$inits, msmdata, qmodel, qcmodel, cmodel, hmodel, p)

#    if (0){
    ## fiddly method using stats::numericDeriv
    likwrap <- function(x, ...){
        pars <- list(unlist(list(...)))
        do.call("lik.msm", c(pars, x))
    }
    myenv <- new.env()
    assign("x", list(msmdata, qmodel, qcmodel, cmodel, hmodel, p), envir = myenv)
    for (i in 1:p$nopt)
        assign(paste("p", i, sep=""), p$inits[i], envir = myenv)
    pvec <- paste("p",1:p$nopt,sep="")
    foo <- numericDeriv(as.call(lapply(as.list(c("likwrap", "x", pvec)), as.name)), pvec, myenv)
    num.d <- attr(foo,"gradient")
    err <- mean(abs(an.d - num.d))
#}
#    err <- num.d <- NULL

    ## much cleaner method. appears to be more accurate as well
    if (requireNamespace("numDeriv", quietly = TRUE))
        numd2 <- numDeriv::grad(lik.msm, p$inits, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel, cmodel=cmodel, hmodel=hmodel, paramdata=p)
    else stop("\"numDeriv\" package not available")
    err2 <- mean(abs(an.d - numd2))

    res <- cbind(analytic=an.d, numeric.base=as.vector(num.d), numeric.nd=numd2)
    rownames(res) <- names(p$inits)
    list(res=res, error=c(base=err, nd=err2))
}
