### Reconstruct data for original msm model fit.  Replace
### generically-named variables in model data frame ("(subject)",
### "(time)", "(state)" etc.) with the names specified by the user for
### the original model call, to allow model refits for
### cross-validation or bootstrapping using the original call.  If not
### found in the usernames attribute, these columns are dropped,
### assuming they're not needed for the refit.  Drop imputed
### observations at piecewise constant intensity change points.
### NOTE assuming timeperiod covariate OK to leave as is

### Factor handling: strips factor() from variable names,
## if original data was factor, will be factor in mf as well, so refit will work
## if not factor in df but factor() in formula, stripping factor() allows refit to work
## TODO can test this with simple refits
## test whether labels are OK, check changelog

reconstruct.data <- function(mf){
    un <- attr(mf, "usernames")
    par.name <- paste0("(",names(un),")")
    ids <- match(par.name, names(mf))
    names(mf) <- replace(names(mf), ids, un)
    if (!is.null(mf$"(pci.imp)")) mf <- mf[!mf$"(pci.imp)",]
    if(is.na(un["subject"])) {
        un["subject"] <- "subject"
        mf$subject <- mf$"(subject)"
        attr(mf, "usernames") <- un
    }
    for (i in grep("^\\(.+\\)$", names(mf), value=TRUE)) mf[,i] <- NULL
    colnames(mf) <- gsub("factor\\((.+)\\)", "\\1", colnames(mf))
    mf
}

### Take a bootstrap sample from the data contained in a fitted msm
### model. Sample pairs of consecutive observations, i.e. independent
### transitions.  Not applicable if model is hidden or some states are
### censored.

bootdata.trans.msm <- function(x) {
    dat <- reconstruct.data(model.frame(x))
    un <- attr(dat, "usernames")
    ## sample random rows of original data, excluding last observation
    inds <- sample(which(duplicated(dat[,un["subject"]], fromLast=TRUE)), replace=TRUE)
    ## make new data by interleaving corresponding "from-state" and "to-state" rows
    ntrans <- length(inds)
    z <- array(c(unlist(dat[inds,]), # "from-state" data
                 unlist(dat[inds+1,])), # "to-state" data
               dim=c(ntrans, ncol(dat), 2))
    data.boot <- matrix(aperm(z, c(3,1,2)), nrow=2*ntrans, ncol=ncol(dat),
                        dimnames=list(NULL,names(dat)))
    data.boot <- as.data.frame(data.boot)
    ## label every transition in new data as from a different subject
    data.boot[,un["subject"]] <- rep(1:ntrans, each=2)
    for (i in which(sapply(dat, is.factor)))
        data.boot[,i] <- factor(data.boot[,i], labels=sort(unique(dat[,i])))
    data.boot
}

### Take a bootstrap sample from the data contained in a fitted msm
### model. Sample subjects. Used for hidden models or models with
### censoring, in which the transitions within a subject are not
### independent.

bootdata.subject.msm <- function(x) {
    dat <- model.frame(x)
    subj.num <- match(dat$"(subject)", unique(dat$"(subject)"))
    subjs <- sample(unique(subj.num), replace=TRUE)
    inds <- new.subj <- NULL
    for (i in seq(along=subjs)) {
        subj.inds <- which(subj.num == subjs[i])
        inds <- c(inds, subj.inds)
        new.subj <- c(new.subj, rep(i, length(subj.inds)))
    }
    data.boot <- dat[inds,]
    data.boot[,"(subject)"] <- new.subj
    data.boot <- reconstruct.data(data.boot)
    data.boot
}

### Given a fitted msm model, draw a bootstrap dataset, refit the
### model, and optionally compute a statistic on the refitted model.
### Repeat B times, store the results in a list.
### msm objects tend to be large, so it is advised to compute a statistic on them by specifying "stat", instead
### of using this function to return a list of refitted msm objects.
### To compute more than one statistic, specify, e.g. stat=function(x)list(stat1(x),stat2(x))

### Some of the arguments to the msm call might be user-defined objects.
### e.g. qmatrix, ematrix, hmodel, ...
### Put in help file that these must be in the working environment.

boot.msm <- function(x, stat=pmatrix.msm, B=1000, file=NULL, cores=NULL){
    boot.fn <- function(dummy){
        boot.data <- if (x$hmodel$hidden || x$cmodel$ncens) bootdata.subject.msm(x) else bootdata.trans.msm(x)
        x$call$data <- substitute(boot.data)
        res <- try(eval(x$call))
        if (!inherits(res, "try-error") && !is.null(stat))
            res <- try(stat(res))
        res
    }
    if (is.null(cores) || cores==1) parallel <- FALSE else parallel <- TRUE;
    if (parallel) {
        if (!is.null(cores) && cores=="default") cores <- NULL
        if (requireNamespace("doParallel", quietly = TRUE)){
### can't get this working separated out into a function like portable.parallel(). Variable exporting / scoping doesnt' work.
            if (.Platform$OS.type == "windows") {
                cl <- parallel::makeCluster(cores)
                doParallel::registerDoParallel(cl)
            } else doParallel::registerDoParallel(cores=cores)
            boot.list <- foreach::"%dopar%"(foreach::foreach(i=1:B, .packages="msm", .export=c("x",ls(.GlobalEnv))),
                                            { boot.fn(i) })
            if (.Platform$OS.type == "windows")
                parallel::stopCluster(cl)
        }
        else stop("\"parallel\" package not available")

    }
    else {
        boot.list <- vector(B, mode="list")
        for (i in 1:B) {
            boot.list[[i]] <- boot.fn()
            if (!is.null(file)) save(boot.list, file=file)
        }
    }
    boot.list
}

### Utilities for calculating bootstrap CIs for particular statistics

qmatrix.ci.msm <- function(x, covariates="mean", sojourn=FALSE, cl=0.95, B=1000, cores=NULL) {
    q.list <- boot.msm(x, function(x)qmatrix.msm(x=x, covariates=covariates)$estimates, B=B, cores=cores)
    q.array <- array(unlist(q.list), dim=c(dim(q.list[[1]]), length(q.list)))
    q.ci <- apply(q.array, c(1,2), function(x)(c(quantile(x, c(0.5 - cl/2, 0.5 + cl/2)), sd(x))))
    q.ci <- aperm(q.ci, c(2,3,1))
    if (sojourn) {
        soj.array <- apply(q.array, 3, function(x) -1/diag(x))
        soj.ci <- apply(soj.array, 1, function(x)(c(quantile(x, c(0.5 - cl/2, 0.5 + cl/2)), sd(x))))
        list(q=q.ci, soj=soj.ci)
    }
    else q.ci
}

ematrix.ci.msm <- function(x, covariates="mean", cl=0.95, B=1000, cores=NULL) {
    e.list <- boot.msm(x, function(x)ematrix.msm(x=x, covariates=covariates)$estimates, B=B, cores=cores)
    e.array <- array(unlist(e.list), dim=c(dim(e.list[[1]]), length(e.list)))
    e.ci <- apply(e.array, c(1,2), function(x)(c(quantile(x, c(0.5 - cl/2, 0.5 + cl/2)), sd(x))))
    aperm(e.ci, c(2,3,1))
}

qratio.ci.msm <- function(x, ind1, ind2, covariates="mean", cl=0.95, B=1000, cores=NULL) {
    q.list <- boot.msm(x, function(x)qratio.msm(x=x, ind1=ind1, ind2=ind2, covariates=covariates)["estimate"], B=B, cores=cores)
    q.vec <- unlist(q.list)
    c(quantile(q.vec, c(0.5 - cl/2, 0.5 + cl/2)), sd(q.vec))
}

pnext.ci.msm <- function(x, covariates="mean", cl=0.95, B=1000, cores=NULL) {
    p.list <- boot.msm(x, function(x)pnext.msm(x=x, covariates=covariates, ci="none")$estimates, B=B, cores=cores)
    p.array <- array(unlist(p.list), dim=c(dim(p.list[[1]]), length(p.list)))
    p.ci <- apply(p.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    aperm(p.ci, c(2,3,1))
}

pmatrix.ci.msm <- function(x, t, t1, covariates="mean", cl=0.95, B=1000, cores=NULL) {
    p.list <- boot.msm(x, function(x)pmatrix.msm(x=x, t=t, t1=t1, covariates=covariates,ci="none"), B=B, cores=cores)
    p.array <- array(unlist(p.list), dim=c(dim(p.list[[1]]), length(p.list)))
    p.ci <- apply(p.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    aperm(p.ci, c(2,3,1))
}

pmatrix.piecewise.ci.msm <- function(x, t1, t2, times, covariates="mean", cl=0.95, B=1000, cores=NULL) {
    p.list <- boot.msm(x, function(x)pmatrix.piecewise.msm(x=x, t1=t1, t2=t2, times=times, covariates=covariates,ci="none"), B=B, cores=cores)
    p.array <- array(unlist(p.list), dim=c(dim(p.list[[1]]), length(p.list)))
    p.ci <- apply(p.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    aperm(p.ci, c(2,3,1))
}

totlos.ci.msm <- function(x, start=1, end=NULL, fromt=0, tot=Inf, covariates="mean", piecewise.times=NULL, piecewise.covariates=NULL, discount=0, env=FALSE, cl=0.95, B=1000, cores=NULL,...) {
    t.list <- boot.msm(x, function(x)totlos.msm(x=x, start=start, end=end, fromt=fromt, tot=tot, covariates=covariates,
                                                piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates,
                                                discount=discount, env=env, ci="none",...), B=B, cores=cores)
    t.array <- do.call("rbind", t.list)
    apply(t.array, 2, function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
}

efpt.ci.msm <- function(x, qmatrix=NULL, tostate, start, covariates="mean", cl=0.95, B=1000, cores=NULL) {
    t.list <- boot.msm(x, function(x)efpt.msm(x=x, qmatrix=qmatrix, start=start, tostate=tostate, covariates=covariates, ci="none"), B=B, cores=cores)
    t.array <- do.call("rbind", t.list)
    apply(t.array, 2, function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
}

ppass.ci.msm <- function(x, qmatrix, tot, start, covariates="mean", piecewise.times=NULL, piecewise.covariates=NULL, cl=0.95, B=1000, cores=NULL,...) {
    t.list <- boot.msm(x, function(x)ppass.msm(x=x, qmatrix=qmatrix, tot=tot, start=start, covariates=covariates,
                                                piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates,
                                                ci="none",...), B=B, cores=cores)
    nst <- ncol(t.list[[1]])
    t.array <- do.call("rbind", lapply(t.list, as.vector))
    ci <- apply(t.array, 2, function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    di <- dimnames(t.list[[1]])
    list(L=matrix(ci[1,],ncol=nst,dimnames=di), U=matrix(ci[2,],ncol=nst, dimnames=di))
}

phasemeans.ci.msm <- function(x, covariates="mean", cl=0.95, B=1000, cores=NULL, ...) {
    p.list <- boot.msm(x, function(x)phasemeans.msm(x=x, covariates=covariates, ci="none", ...), B=B, cores=cores)
    p.array <- array(unlist(p.list), dim=c(dim(p.list[[1]]), length(p.list)))
    p.ci <- apply(p.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    aperm(p.ci, c(2,3,1))
}

expected.ci.msm <- function(x,
                            times=NULL,
                            timezero=NULL,
                            initstates=NULL,
                            covariates="mean",
                            misccovariates="mean",
                            piecewise.times=NULL,
                            piecewise.covariates=NULL,
                            risk=NULL,
                            cl=0.95, B=1000, cores=NULL) {
    if(is.null(risk)) risk <- observed.msm(x)$risk
    e.list <- boot.msm(x, function(x){
        expected.msm(x, times, timezero, initstates, covariates, misccovariates, piecewise.times, piecewise.covariates, risk)
    }, B=B, cores=cores)
    e.tab.array <- array(unlist(lapply(e.list, function(x)x[[1]])), dim=c(dim(e.list[[1]][[1]]), length(e.list)))
    e.perc.array <- array(unlist(lapply(e.list, function(x)x[[2]])), dim=c(dim(e.list[[1]][[2]]), length(e.list)))
    e.tab.ci <- apply(e.tab.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    e.perc.ci <- apply(e.perc.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    res <- list(aperm(e.tab.ci, c(2,3,1)),  aperm(e.perc.ci, c(2,3,1)))
    names(res) <- c("Expected", "Expected percentages")
    res
}

### Compute a CI for a statistic using a sample from the assumed MVN
### distribution of MLEs of log Q, logit E and covariate effects on these
### Not user visible: only support statistics based on Q matrix and E matrix
### i.e. statistics computed as functions of x$Qmatrices, x$Ematrices and x$paramdata$params

normboot.msm <- function(x, stat, B=1000) {
    ## simulate from vector of unreplicated parameters, to avoid numerical problems with rmvnorm when lots of correlations are 1
    if (!x$foundse) stop("Asymptotic standard errors not available in fitted model")
    sim <- rmvnorm(B, x$opt$par, x$covmat[x$paramdata$optpars,x$paramdata$optpars])
    params <- matrix(nrow=B, ncol=x$paramdata$npars)  # replicate constrained parameters.
    params[,x$paramdata$optpars] <- sim
    params[,x$paramdata$fixedpars] <- rep(x$paramdata$params[x$paramdata$fixedpars], each=B)
    params[,x$paramdata$hmmpars] <- rep(msm.mninvlogit.transform(x$paramdata$params[x$paramdata$hmmpars], x$hmodel$plabs, x$hmodel$parstate), each=B)
    params <- params[, !duplicated(abs(x$paramdata$constr)), drop=FALSE][, abs(x$paramdata$constr), drop=FALSE] *
        rep(sign(x$paramdata$constr), each=B)

    sim.stat <- vector(B, mode="list")
    for (i in 1:B) {
        x.rep <- x
        x.rep$paramdata$params <- params[i,]
        output <- msm.form.output(x.rep, "intens")
        x.rep$Qmatrices <- output$Qmatrices
        if (x$emodel$misc) {
            output <- msm.form.output(x.rep, "misc")
            x.rep$Ematrices <- output$Ematrices
            names(x.rep$Ematrices)[1] <- "logitbaseline"
        }
        sim.stat[[i]] <- stat(x.rep)
    }
    sim.stat
}

qmatrix.normci.msm <- function(x, covariates="mean", sojourn=FALSE, cl=0.95, B=1000) {
    q.list <- normboot.msm(x, function(x)qmatrix.msm(x=x, covariates=covariates, ci="none"), B)
    q.array <- array(unlist(q.list), dim=c(dim(q.list[[1]]), length(q.list)))
    q.ci <- apply(q.array, c(1,2), function(x)(c(quantile(x, c(0.5 - cl/2, 0.5 + cl/2)), sd(x))))
    q.ci <- aperm(q.ci, c(2,3,1))
    if (sojourn) {
        soj.array <- apply(q.array, 3, function(x) -1/diag(x))
        soj.ci <- apply(soj.array, 1, function(x)(c(quantile(x, c(0.5 - cl/2, 0.5 + cl/2)), sd(x))))
        list(q=q.ci, soj=soj.ci)
    }
    else q.ci
}

ematrix.normci.msm <- function(x, covariates="mean", cl=0.95, B=1000) {
    e.list <- normboot.msm(x, function(x)ematrix.msm(x=x, covariates=covariates, ci="none"), B)
    e.array <- array(unlist(e.list), dim=c(dim(e.list[[1]]), length(e.list)))
    e.ci <- apply(e.array, c(1,2), function(x)(c(quantile(x, c(0.5 - cl/2, 0.5 + cl/2)), sd(x))))
    aperm(e.ci, c(2,3,1))
}

qratio.normci.msm <- function(x, ind1, ind2, covariates="mean", cl=0.95, B=1000) {
    q.list <- normboot.msm(x, function(x)qratio.msm(x=x, ind1=ind1, ind2=ind2, covariates=covariates, ci="none")["estimate"], B)
    q.vec <- unlist(q.list)
    c(quantile(q.vec, c(0.5 - cl/2, 0.5 + cl/2)), sd(q.vec))
}

pnext.normci.msm <- function(x, covariates="mean", cl=0.95, B=1000) {
    p.list <- normboot.msm(x, function(x)pnext.msm(x=x, covariates=covariates, ci="none")$estimates, B)
    p.array <- array(unlist(p.list), dim=c(dim(p.list[[1]]), length(p.list)))
    p.ci <- apply(p.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    aperm(p.ci, c(2,3,1))
}

pmatrix.normci.msm <- function(x, t, t1, covariates="mean", cl=0.95, B=1000) {
    p.list <- normboot.msm(x, function(x)pmatrix.msm(x=x, t=t, t1=t1, covariates=covariates, ci="none"), B)
    p.array <- array(unlist(p.list), dim=c(dim(p.list[[1]]), length(p.list)))
    p.ci <- apply(p.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    aperm(p.ci, c(2,3,1))
}

pmatrix.piecewise.normci.msm <- function(x, t1, t2, times, covariates="mean", cl=0.95, B=1000) {
    p.list <- normboot.msm(x, function(x)pmatrix.piecewise.msm(x=x, t1=t1, t2=t2, times=times, covariates=covariates, ci="none"), B)
    p.array <- array(unlist(p.list), dim=c(dim(p.list[[1]]), length(p.list)))
    p.ci <- apply(p.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    aperm(p.ci, c(2,3,1))
}

totlos.normci.msm <- function(x, start=1, end=NULL, fromt=0, tot=Inf, covariates="mean", piecewise.times=NULL, piecewise.covariates=NULL, discount=0, env=FALSE, cl=0.95, B=1000, ...) {
    t.list <- normboot.msm(x, function(x)totlos.msm(x=x, start=start, end=end, fromt=fromt, tot=tot, covariates=covariates,
                                                    piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates,
                                                    discount=discount, env=env, ci="none", ...), B)
    t.array <- do.call("rbind", t.list)
    apply(t.array, 2, function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
}

efpt.normci.msm <- function(x, qmatrix=NULL, tostate, start, covariates="mean", cl=0.95, B=1000, ...) {
    t.list <- normboot.msm(x, function(x)efpt.msm(x=x, qmatrix=qmatrix, tostate=tostate, start=start, covariates=covariates, ci="none", ...), B)
    t.array <- do.call("rbind", t.list)
    apply(t.array, 2, function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
}

ppass.normci.msm <- function(x, qmatrix, tot, start, covariates="mean", piecewise.times=NULL, piecewise.covariates=NULL, cl=0.95, B=1000, ...) {
    t.list <- normboot.msm(x, function(x)ppass.msm(x=x, qmatrix=qmatrix, tot=tot, start=start, covariates=covariates,
                                                    piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates,
                                                    ci="none", ...), B)
    nst <- ncol(t.list[[1]])
    t.array <- do.call("rbind", lapply(t.list, as.vector))
    ci <- apply(t.array, 2, function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    di <- dimnames(t.list[[1]])
    list(L=matrix(ci[1,],ncol=nst,dimnames=di), U=matrix(ci[2,],ncol=nst, dimnames=di))
}

expected.normci.msm <- function(x,
                                times=NULL,
                                timezero=NULL,
                                initstates=NULL,
                                covariates="mean",
                                misccovariates="mean",
                                piecewise.times=NULL,
                                piecewise.covariates=NULL,
                                risk=NULL,
                                cl=0.95, B=1000) {
    if(is.null(risk)) risk <- observed.msm(x)$risk
    e.list <- normboot.msm(x, function(x){
        expected.msm(x, times, timezero, initstates, covariates, misccovariates, piecewise.times, piecewise.covariates, risk)
    }, B)
    e.tab.array <- array(unlist(lapply(e.list, function(x)x[[1]])), dim=c(dim(e.list[[1]][[1]]), length(e.list)))
    e.perc.array <- array(unlist(lapply(e.list, function(x)x[[2]])), dim=c(dim(e.list[[1]][[2]]), length(e.list)))
    e.tab.ci <- apply(e.tab.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    e.perc.ci <- apply(e.perc.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    res <- list(aperm(e.tab.ci, c(2,3,1)),  aperm(e.perc.ci, c(2,3,1)))
    names(res) <- c("Expected", "Expected percentages")
    res
}

phasemeans.normci.msm <- function(x, covariates="mean", cl=0.95, B=1000, ...) {
    p.list <- normboot.msm(x, function(x)phasemeans.msm(x=x, covariates=covariates, ci="none", ...), B)
    p.array <- array(unlist(p.list), dim=c(dim(p.list[[1]]), length(p.list)))
    p.ci <- apply(p.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    aperm(p.ci, c(2,3,1))
}
