draic.msm <- function(msm.full,msm.coarse,
                  likelihood.only=FALSE,
                  information=c("expected", "observed"),
                  tl=0.95
                  )
{
    if (!inherits(msm.full, "msm")) stop("Expected \"msm.full\" to be a msm model")
    if (!inherits(msm.coarse, "msm")) stop("Expected \"msm.coarse\" to be a msm model")
    if (msm.full$hmodel$hidden || msm.coarse$hmodel$hidden) stop("Hidden Markov models not supported")
    if (msm.full$cmodel$ncens > 0 || msm.coarse$cmodel$ncens >0) stop("Models with censoring not supported")
    ## Reconstruct datasets used to fit each model as dataframes. standard names will be bracketed.
    msm.full.data <- reconstruct.data(model.frame(msm.full))
    msm.coarse.data <- reconstruct.data(model.frame(msm.coarse))
    un <- attr(msm.full.data, "usernames")
    subj <- un["subject"]; st <- un["state"]
    unc <- attr(msm.coarse.data, "usernames")
    subjc <- unc["subject"]; stc <- unc["state"]
    if (nrow(msm.full.data) != nrow(msm.coarse.data)) stop("Full dataset has ",nrow(msm.full.data)," rows, but coarse dataset has ",nrow(msm.coarse.data)," rows, these should be equal.")
    
    ## Determine map from full to coarsened data
    coarsening.ematrix <- matrix(0,nrow=msm.full$qmodel$nstates, ncol=msm.full$qmodel$nstates)
    map <- unique(cbind(full=msm.full.data[,st], coarse=msm.coarse.data[,stc]))
    if (!all(tapply(map[,"coarse"], map[,"full"], function(x)length(unique(x)))==1)){
        message("Observed the following state mappings in the data:")
        print(map); stop("Coarse state space not an aggregation of the the full state space")
    }
    coarsening.ematrix[map] <- 1
    ## Could also check for compatible Q matrices, but would be fiddly, so don't bother
    
    ## Initial state occupancy probabilities, assumed equal across possible states
    initstate <- msm.full.data[,st][!duplicated(msm.full.data[,subj])]
    init.probs <- prop.table(table(factor(initstate, levels=1:msm.full$qmodel$nstates)))
    initp <- matrix(0, nrow=attr(msm.full.data,"npts"), ncol=msm.full$qmodel$nstates)
    for(i in 1:length(initstate)) {
        j<-coarsening.ematrix[initstate[i],]==1
        initp[i,] <- init.probs * coarsening.ematrix[,j]
        initp[i,] <- initp[i,] / sum(initp[i,])
    }

    ## Log-likelihood of complex model on restricted dataset
    lh.restricted.obj<-function(pars, by.subject=FALSE, return.model=FALSE) {
        p <- msm.full$estimates
        p[msm.full$paramdata$optpars] <- pars
        p[names(p)=="qbase"] <- exp(p[names(p)=="qbase"])
        call <- msm.full$call
        ## name of state might be different in coarse data, assumes all other names same
        call$formula <- as.formula(paste(stc, "~", unc["time"]))
        call$data <- substitute(msm.coarse.data)
        call$ematrix <- coarsening.ematrix
        qmat <- t(msm.full$qmodel$imatrix) # Supply parameters as a vector for use in optimHess()
        qmat[qmat==1] <- p[names(p)=="qbase"]
        call$qmatrix <- t(qmat)
        call$covinits <- split(p[names(p)=="qcov"],
                               rep(msm.full$qcmodel$covlabels, each=msm.full$qmodel$npars))
        call$fixedpars <- TRUE
        call$initprobs <- initp
        res <- eval(call)
        if (return.model) res else as.numeric(logLik.msm(res, by.subject=by.subject))
    }

    n <- attr(msm.full.data,"npts")
    mle <- msm.full$estimates[msm.full$paramdata$optpars] # handle models with some parameters fixed
    lsm <- lh.restricted.obj(mle) # Loglik of complex model on restricted data
    lmm <- logLik.msm(msm.full)  # Complex model on full data
    lss <- logLik.msm(msm.coarse) # Simple model on restricted data
    dfit <- (1/n)*(lss - lsm)
    res.lik <- c("complex"=-lsm, "simple"=-lss, "complex-simple"=-(lsm-lss),"(complex-simple)/n"=dfit)

    if (!likelihood.only) {
        information <- match.arg(information)
        if (!all(msm.full$data$mf$"(obstype)" == 1)) information <- "observed"
        if (information=="observed") {
            ## Information matrix of complex model on complex data
            I.complex <- -0.5*(1/n)*msm.full$opt$hessian
            ## Information matrix of complex model on restricted data
            I.restricted <- (1/n)*optimHess(mle, lh.restricted.obj)
        } else {
            I.complex <- -0.5*(1/n)*msm.full$paramdata$information
            res <- lh.restricted.obj(mle, return.model=TRUE)
            ## simply using I <- res$paramdata$information doesn't work for models with some parameters fixed
            ## since call producing res unfixes these to produce derivs/info.
            ## so calculate information separately after keeping these fixed
            res$paramdata$optpars <- msm.full$paramdata$optpars
            res$paramdata$fixedpars <- c(msm.full$paramdata$fixedpars, res$paramdata$auxpars)
            I <- information.msm(res$paramdata$allinits[res$paramdata$optpars], expand.data(res), res$qmodel,
                                 res$qcmodel, res$cmodel, res$hmodel, res$paramdata)
            I.restricted <- -0.5*(1/n)*I
        }
        npars.restricted <- sum(diag(I.restricted %*% solve(I.complex)))
        dcompl <- (1/n)*(npars.restricted - msm.coarse$paramdata$npars)
        DRAIC <- dfit + dcompl

        ## Tracking interval.
        ## Difference in log-likelihoods on restricted data
        dlh <- lh.restricted.obj(mle, by.subject=TRUE) - logLik.msm(msm.coarse, by.subject=TRUE)
        omega.sq <- mean(dlh*dlh) - mean(dlh)^2
        half.width <- qnorm(1 - 0.5*(1 - tl)) * sqrt(omega.sq / n)
        prob.DRAIC <- pnorm(-DRAIC / sqrt(omega.sq / n))
        res.int <- c("2.5%"=DRAIC-half.width,"97.5%"=DRAIC+half.width,"Prob<0"=prob.DRAIC)
        if (!isTRUE(all.equal(tl,0.95))) names(res.int)[1:2] <- c("Lower","Upper")
                          
        liks <- cbind("-LL" = res.lik,
                      npars=c(npars.restricted, msm.coarse$paramdata$npars,
                      npars.restricted - msm.coarse$paramdata$npars, dcompl))
        res <- list(lik.restricted = liks, draic = as.numeric(DRAIC), ti=res.int)
    }
    else res <- res.lik
    
    res
}


drlcv.msm <- function(msm.full,msm.coarse,tl=0.95,cores=NULL,verbose=TRUE,outfile=NULL)
{
    if (!inherits(msm.full, "msm")) stop("Expected \"msm.full\" to be a msm model")
    if (!inherits(msm.coarse, "msm")) stop("Expected \"msm.coarse\" to be a msm model")
    if (msm.full$hmodel$hidden || msm.coarse$hmodel$hidden) stop("Hidden Markov models not supported")
    if (msm.full$cmodel$ncens > 0 || msm.coarse$cmodel$ncens >0) stop("Models with censoring not supported")
    ## Reconstruct datasets used to fit each model as dataframes
    msm.full.data <- reconstruct.data(model.frame(msm.full))
    msm.coarse.data <- reconstruct.data(model.frame(msm.coarse))
    un <- attr(msm.full.data, "usernames")
    subj <- un["subject"]; st <- un["state"]
    unc <- attr(msm.coarse.data, "usernames")
    subjc <- unc["subject"]; stc <- unc["state"]
    if (nrow(msm.full.data) != nrow(msm.coarse.data)) stop("Full dataset has ",nrow(msm.full.data)," rows, but coarse dataset has ",nrow(msm.coarse.data)," rows, these should be equal.")
    
    ## Determine map from full to coarsened data
    coarsening.ematrix <- matrix(0,nrow=msm.full$qmodel$nstates, ncol=msm.full$qmodel$nstates)
    map <- unique(cbind(full=msm.full.data[,st], coarse=msm.coarse.data[,stc]))
    if (!all(tapply(map[,"coarse"], map[,"full"], function(x)length(unique(x)))==1)){
        message("Observed the following state mappings in the data:")
        print(map); stop("Coarse state space not an aggregation of the the full state space")
    }
    coarsening.ematrix[map] <- 1

    ## Initial state occupancy probabilities assumed equal across possible states
    initstate <- msm.full.data[,st][!duplicated(msm.full.data[,subj])]
    init.probs <- prop.table(table(factor(initstate, levels=1:msm.full$qmodel$nstates)))
    initp <- matrix(0, nrow=attr(msm.full.data,"npts"), ncol=msm.full$qmodel$nstates)
    for(i in 1:length(initstate)) {
        j<-coarsening.ematrix[initstate[i],]==1
        initp[i,] <- init.probs * coarsening.ematrix[,j]
        initp[i,] <- initp[i,] / sum(initp[i,])
    }

    ## Log-likelihood of complex model on restricted dataset
    call <- msm.full$call
    ## name of state might be different in coarse data, assumes all other names same
    call$formula <- as.formula(paste(stc, "~", unc["time"]))
    call$data <- substitute(msm.coarse.data)
    call$ematrix <- coarsening.ematrix
    call$qmatrix <- qmatrix.msm(msm.full, ci="none") # start it at the full MLE to help convergence
    call$covinits <- lapply(msm.full$Qmatrices[msm.full$qcmodel$covlabels],
                            function(x){t(x)[msm.full$qmodel$imatrix==1]})
    call$fixedpars <- TRUE
    call$initprobs <- initp
    msm.full.on.restricted <- eval(call)

    n <- attr(msm.full.data,"npts")
    lsm <- logLik.msm(msm.full.on.restricted) # Loglik of complex model on restricted data
    lmm <- logLik.msm(msm.full)  # Complex model on full data
    lss <- logLik.msm(msm.coarse) # Simple model on restricted data

    ## Cross validation iteration for leaving out ith subject
    cv.fn <- function(i) {
        ## Refit small model after leaving one subject out
        call <- msm.coarse$call
        call$data <- msm.coarse.data[msm.coarse.data[,subjc] != unique(msm.coarse.data[,subjc])[i],]
        call$qmatrix <- qmatrix.msm(msm.coarse, ci="none") # start it at the full MLE to help convergence
        call$covinits <- lapply(msm.coarse$Qmatrices[msm.coarse$qcmodel$covlabels],
                                function(x){t(x)[msm.coarse$qmodel$imatrix==1]})
        res.coarse <- try(eval(call))
        lssf <- as.numeric(logLik.msm(res.coarse))

        ## Loglik of refitted model on left-out subject
        call$data <- msm.coarse.data[msm.coarse.data[,subjc] == unique(msm.coarse.data[,subjc])[i],]
        call$qmatrix <- qmatrix.msm(res.coarse, ci="none", covariates=0)
        call$center <- FALSE
        call$covinits <- lapply(res.coarse$Qmatrices[res.coarse$qcmodel$covlabels],
                                function(x){t(x)[res.coarse$qmodel$imatrix==1]})
        ## Workaround to avoid dropping unused levels and breaking model.matrix when factor() used around an existing factor in the data
        if (!is.null(call$covariates)){
            call$covariates <- as.formula(gsub("factor\\((.+)\\)", "factor(\\1,levels=levels(\\1))", call$covariates))
            names(call$covinits) <- gsub("factor\\((.+)\\)", "factor(\\1, levels = levels(\\1))", names(call$covinits))
        }
        if (!is.null(call$constraint))
            names(call$constraint) <- gsub("factor\\((.+)\\)", "factor(\\1, levels = levels(\\1))", names(call$constraint))
        call$fixedpars <- TRUE
        res <- suppressWarnings(eval(call)) ## suppress warnings about state vector not containing observations of all states
        lss <- as.numeric(logLik.msm(res))

        ## Refit big model after leaving one subject out.
        call <- msm.full$call
        call$data <- msm.full.data[msm.full.data[,subj] != unique(msm.full.data[,subj])[i],]
        call$qmatrix <- qmatrix.msm(msm.full, ci="none") # start it at the full MLE to help convergence
        call$covinits <- lapply(msm.full$Qmatrices[msm.full$qcmodel$covlabels],
                                function(x){t(x)[msm.full$qmodel$imatrix==1]})
        res.full <- try(eval(call))
        lsmf <- as.numeric(logLik.msm(res.full))

        ## Loglik of refitted model on left-out subject, using coarsened data
        call$formula <- as.formula(paste(stc, "~", unc["time"]))
        call$data <- msm.coarse.data[msm.coarse.data[,subjc] == unique(msm.coarse.data[,subjc])[i],]
        call$ematrix <- coarsening.ematrix
        call$qmatrix <- qmatrix.msm(res.full, ci="none", covariates=0)
        call$center <- FALSE
        call$covinits <- lapply(res.full$Qmatrices[res.full$qcmodel$covlabels],
                                function(x){t(x)[res.full$qmodel$imatrix==1]})
        if (!is.null(call$covariates)){
            call$covariates <- as.formula(gsub("factor\\((.+)\\)", "factor(\\1,levels=levels(\\1))", call$covariates))
            names(call$covinits) <- gsub("factor\\((.+)\\)", "factor(\\1, levels = levels(\\1))", names(call$covinits))
        }
        if (!is.null(call$constraint))
            names(call$constraint) <- gsub("factor\\((.+)\\)", "factor(\\1, levels = levels(\\1))", names(call$constraint))
        call$initprobs <- initp[i,,drop=FALSE]
        call$fixedpars <- TRUE
        res <- eval(call)
        lsm <- as.numeric(logLik.msm(res))

        conv <- ifelse(res.full$opt$convergence==0, "converged", "NOT CONVERGED")
        resc <- sprintf("i=%d, smallrest = %.12f, smalli = %.12f, bigrest = %.12f, bigi = %.12f, diff = %.12f, %s\n",
                        i,  lssf, lss, lsmf, lsm,  lss - lsm, conv)
        if (verbose) cat(resc)
        if (!is.null(outfile)) cat(resc, file=outfile, append=TRUE)

        c(lsm=lsm, lss=lss, lssf=lssf, lsmf=lsmf)
    }
    if (is.null(cores) || cores==1) parallel <- FALSE else parallel <- TRUE;
    if (parallel){
        if (!is.null(cores) && cores=="default") cores <- NULL
        if (requireNamespace("doParallel", quietly = TRUE)){
### can't get this working separated out into a function like portable.parallel(). Variable exporting / scoping doesnt' work.
### no need to export these?   "msm.coarse", "msm.coarse.data", "msm.full", "msm.full.data","coarsening.ematrix", "initp"
            if (.Platform$OS.type == "windows") {
                cl <- parallel::makeCluster(cores)
                doParallel::registerDoParallel(cl)
            } else doParallel::registerDoParallel(cores=cores)
            cv <- foreach::"%dopar%"(foreach::foreach(i=1:n, .packages="msm", .export=c(ls(.GlobalEnv))),
                                     { cv.fn(i) })
        } else stop("\"parallel\" package not available")
    }
    else {
        cv <- vector(mode="list", n)
        for (i in 1:n){
            cv[[i]] <- cv.fn(i)
        }
    }
    LCVi <- sapply(cv, function(x){x["lss"] - x["lsm"]})
    DRLCV <- mean(LCVi)

    dlh <- logLik.msm(msm.full.on.restricted, by.subject=TRUE) -
        logLik.msm(msm.coarse, by.subject=TRUE)
    omega.sq <- mean(dlh*dlh) - mean(dlh)^2
    half.width <- qnorm(1 - 0.5*(1 - tl)) * sqrt(omega.sq / n)
    prob.DRLCV <- pnorm(-DRLCV / sqrt(omega.sq / n))
    res.int <- c("2.5%"=DRLCV-half.width,"97.5%"=DRLCV+half.width,"Prob<0"=prob.DRLCV)
    if (!isTRUE(all.equal(tl,0.95))) names(res.int)[1:2] <- c("Lower","Upper")
    liks <- cbind("-LL" = c(-lsm, -lss, -(lsm-lss), -(lsm-lss)/n))
    rownames(liks) <- c("complex","simple","complex-simple","(complex-simple)/n")
    res <- list(lik.restricted = liks, drlcv = as.numeric(DRLCV), ti=res.int)
    res
}
