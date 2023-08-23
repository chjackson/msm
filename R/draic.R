#' Criteria for comparing two multi-state models with nested state spaces
#' 
#' A modification of Akaike's information criterion, and a leave-one-out
#' likelihood cross-validation criterion, for comparing the predictive ability
#' of two Markov multi-state models with nested state spaces.  This is
#' evaluated based on the restricted or aggregated data which the models have
#' in common.
#' 
#' Note that standard AIC can be computed for one or more fitted \code{msm}
#' models \code{x,y,...} using \code{\link{AIC}(x,y,...)}, and this can be used
#' to compare models fitted to the same data. \code{draic.msm} and
#' \code{drlcv.msm} are designed for models fitted to data with
#' differently-aggregated state spaces.
#' 
#' 
#' The difference in restricted AIC (Liquet and Commenges, 2011), as computed
#' by this function, is defined as
#' 
#' \deqn{D_{RAIC} = l(\gamma_n |\mathbf{x}'' ) - l(\theta_n |\mathbf{x}'' ) +
#' trace ( J(\theta_n |\mathbf{x}'')J(\theta_n |\mathbf{x})^{-1} - J(\gamma_n
#' |\mathbf{x}'' )J(\gamma_n |\mathbf{x}' )^{-1})}{D_RAIC = l(gamma_n |x'' ) -
#' l(theta_n |x'' ) + trace ( J(theta_n |x'')J(theta_n |x)^{-1} - J(gamma_n
#' |x'' )J(gamma_n |x' )^{-1})}
#' 
#' where \eqn{\gamma}{gamma} and \eqn{\theta}{theta} are the maximum likelihood
#' estimates of the smaller and bigger models, fitted to the smaller and bigger
#' data, respectively.
#' 
#' \eqn{l(\gamma_n |x'')}{l(gamma_n |x'')} represents the likelihood of the
#' simpler model evaluated on the restricted data.
#' 
#' \eqn{l(\theta_n |x'')}{l(theta_n |x'')} represents the likelihood of the
#' complex model evaluated on the restricted data.  This is a hidden Markov
#' model, with a misclassification matrix and initial state occupancy
#' probabilities as described by Thom et al (2014).
#' 
#' \eqn{J()} are the corresponding (expected or observed, as specified by the
#' user) information matrices.
#' 
#' \eqn{\mathbf{x}}{x} is the expanded data, to which the bigger model was
#' originally fitted, and \eqn{\mathbf{x}'}{x'} is the data to which the
#' smaller model was originally fitted.  \eqn{\mathbf{x}''}{x''} is the
#' restricted data which the two models have in common.  \eqn{\mathbf{x}'' =
#' \mathbf{x}'}{x'' = x} in this implementation, so the models are nested.
#' 
#' The difference in likelihood cross-validatory criteria (Liquet and
#' Commenges, 2011) is defined as
#' 
#' \deqn{D_{RLCV} = 1/n \sum_{i=1}^n \log( h_{X''}(x_i'' | \gamma_{-i}) /
#' g_{X''}(x_i''| \theta_{-i}))}{D_{RLCV} = 1/n \sum_{i=1}^n log( h_{X''}(x_i''
#' | gamma_{-i}) / g_{X''}(x_i''| theta_{-i}))}
#' 
#' where \eqn{\gamma_{-i}} and \eqn{\theta_{-i}} are the maximum likelihood
#' estimates from the smaller and bigger models fitted to datasets with subject
#' \eqn{i} left out, \eqn{g()} and \eqn{h()} are the densities of the
#' corresponding models, and \eqn{x_i''} is the restricted data from subject
#' \eqn{i}.
#' 
#' Tracking intervals are analogous to confidence intervals, but not strictly
#' the same, since the quantity which D_RAIC aims to estimate, the difference
#' in expected Kullback-Leibler discrepancy for predicting a replicate dataset,
#' depends on the sample size.  See the references.
#' 
#' Positive values for these criteria indicate the coarse model is preferred,
#' while negative values indicate the full model is preferred.
#' 
#' @aliases draic.msm drlcv.msm
#' @param msm.full Model on the bigger state space.
#' @param msm.coarse Model on the smaller state space.
#' 
#' The two models must both be non-hidden Markov models without censored
#' states.
#' 
#' The two models must be fitted to the same datasets, except that the state
#' space of the coarse model must be an aggregated version of the state space
#' of the full model.  That is, every state in the full dataset must correspond
#' to a unique state in the coarse dataset.  For example, for the full state
#' variable \code{c(1,1,2,2,3,4)}, the corresponding coarse states could be
#' \code{c(1,1,2,2,2,3)}, but not \code{c(1,2,3,4,4,4)}.
#' 
#' The structure of allowed transitions in the coarse model must also be a
#' collapsed version of the big model structure, but no check is currently made
#' for this in the code.
#' 
#' To use these functions, all objects which were used in the calls to fit
#' \code{msm.full} and \code{msm.coarse} must be in the working environment,
#' for example, datasets and definitions of transition matrices.
#' @param likelihood.only Don't calculate Hessians and trace term (DRAIC).
#' @param information Use observed or expected information in the DRAIC trace
#' term.  Expected is the default, and much faster, though is only available
#' for models fitted to pure panel data (all \code{obstype=1} in the call to
#' \code{\link{msm}}, thus not exact transition times or exact death times)
#' @param tl Width of symmetric tracking interval, by default 0.95 for a 95\%
#' interval.
#' @param cores Number of processor cores to use in \code{drlcv} for
#' cross-validation by parallel processing.  Requires the \pkg{doParallel}
#' package to be installed.  If not specified, parallel processing is not used.
#' If \code{cores} is set to the string \code{"default"}, the default methods
#' of \code{\link[parallel]{makeCluster}} (on Windows) or
#' \code{\link[doParallel]{registerDoParallel}} (on Unix-like) are used.
#' @param verbose Print intermediate results of each iteration of
#' cross-validation to the console while running. May not work with parallel
#' processing.
#' @param outfile Output file to print intermediate results of
#' cross-validation.  Useful to track execution speed when using parallel
#' processing, where output to the console may not work.
#' @return A list containing \eqn{D_{RAIC}}{D_RAIC} (\code{draic.msm}) or
#' \eqn{D_{RLCV}}{D_RLCV} (\code{drlcv.msm}), its component terms, and tracking
#' intervals.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}, H. H. Z.
#' Thom \email{howard.thom@@bristol.ac.uk}
#' @seealso \code{\link{logLik.msm}}
#' @references Thom, H. and Jackson, C. and Commenges, D. and Sharples, L.
#' (2015) State selection in multistate models with application to quality of
#' life in psoriatic arthritis.  Statistics In Medicine 34(16) 2381 - 2480.
#' 
#' Liquet, B. and Commenges D. (2011) Choice of estimators based on different
#' observations: Modified AIC and LCV criteria. Scandinavian Journal of
#' Statistics; 38:268-287.
#' @keywords models
#' @export draic.msm
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


#' @rdname draic.msm
#' @export
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
