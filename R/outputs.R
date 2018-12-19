### METHODS FOR MSM OBJECTS

qmatrix.msm <- function(x, covariates="mean", sojourn=FALSE, ci=c("delta","normal","bootstrap","none"), cl=0.95, B=1000, cores=NULL)
{   
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    nst <- x$qmodel$nstates
    ni <- x$qmodel$npars
    covlist <- msm.parse.covariates(x, covariates, x$qcmodel)
    nc <- length(covlist)  # number of effects we need to adjust baseline for
    if ((cl < 0) || (cl > 1)) stop("expected cl in [0,1]")
    se <- lse <- fixed <- numeric(ni)
    logest <- x$Qmatrices$logbaseline
    if (is.null(x$QmatricesFixed)) x <- msm.form.output(x, "intens") # for back-compat with pre 1.4.1 model objects
    fixed <- x$QmatricesFixed$logbaseline
    for (i in seq_len(nc)) {
        logest <- logest + x$Qmatrices[[i+1]] * covlist[[i]]
        fixed <- fixed & x$QmatricesFixed[[i+1]]  # Only true if all coefficients contributing to the estimate are fixed.  Used in print functions
    }
    mat <- exp(logest)
    mat[x$qmodel$imatrix == 0] <- 0
    mat <- msm.fixdiag.qmatrix(mat)

    if (sojourn) soj = -1 / diag(mat)
    ci <- match.arg(ci)
    if (x$foundse && (ci!="none")) {
        if (ci == "delta") {
            ## Work out standard errors.
            ## Transformation for delta method is (intensities)
            ##  exp (x1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... )
            ## expit(sum covs)  / (1 + expit(sum(covs)))    or   1  /  (1  +  expit(sum(covs)))
            ## Use delta method to find approximate SE of the transform on log scale
            ## Work out a CI for this by assuming normal and transforming back
            coefs <- as.numeric(c(1, covlist))
            semat <- lsemat <- lmat <- umat <- matrix(0, nst, nst)
            form <- as.formula(paste("~", expsum(seq(nc + 1), coefs)))
            lform <- as.formula(paste("~", lsum(seq(nc + 1), coefs)))
            ## indices into estimates vector of all intens/miscs, intens covs / misc covs
            inds <- seq(length=x$qmodel$npars + x$qcmodel$npars)
            for (i in 1 : ni){
                ## indices into estimates vector of all intens/miscs, intens covs / misc covs for that particular fromstate-tostate.
                parinds <- inds[seq(i, (nc * ni + i), ni)]
                ests <- x$estimates[parinds]
                cov <- x$covmat[parinds, parinds]
                se[i] <- deltamethod(form, ests, cov)
                lse[i] <- deltamethod(lform, ests, cov)
            }
            ivector <- as.numeric(t(x$qmodel$imatrix))
            semat[ivector == 1] <- se; semat <- t(semat)
            lsemat[ivector == 1] <- lse; lsemat <- t(lsemat)
            lmat <- exp(logest - qnorm(1 - 0.5*(1 - cl))*lsemat)
            umat <- exp(logest + qnorm(1 - 0.5*(1 - cl))*lsemat)
            imatrix <- x$qmodel$imatrix
            lmat[imatrix == 0] <- umat[imatrix == 0] <- 0
            ## SEs of diagonal entries
            diagse <- qmatrix.diagse.msm(x, covlist, sojourn, ni, ivector, nc)
            diag(semat) <- diagse$diagse
            diag(lmat) <- sign(diag(mat)) * (exp(log(abs(diag(mat))) - sign(diag(mat)) * qnorm(1 - 0.5*(1 - cl))*diagse$diaglse))
            diag(umat) <- sign(diag(mat)) * (exp(log(abs(diag(mat))) + sign(diag(mat)) * qnorm(1 - 0.5*(1 - cl))*diagse$diaglse))
            if (sojourn) {
                sojse <- diagse$sojse
                sojl <- exp(log(soj) - qnorm(1 - 0.5*(1 - cl))*diagse$sojlse)
                soju <- exp(log(soj) + qnorm(1 - 0.5*(1 - cl))*diagse$sojlse)
            }
        }
        else if (ci %in% c("normal","bootstrap")) {
            q.ci <- if (ci=="normal")
                        qmatrix.normci.msm(x, covariates, sojourn, cl, B) else qmatrix.ci.msm(x, covariates, sojourn, cl, B, cores)
            if (sojourn) {
                soj.ci <- q.ci$soj
                q.ci <- q.ci$q
                sojl <- soj.ci[1,]; soju <- soj.ci[2,]; sojse <- soj.ci[3,]
            }
            lmat <- q.ci[,,1]; umat <- q.ci[,,2]; semat <- q.ci[,,3]
        }
        dimnames(semat) <- dimnames(lmat) <- dimnames(umat) <- dimnames(x$qmodel$qmatrix)
    }
    else semat <- lmat <- umat <- sojse <- sojl <- soju <- NULL

    dimnames(mat) <-  dimnames(x$qmodel$qmatrix)
    if (ci=="none") res <- if (sojourn) soj else mat
    else {
        if (sojourn)
            res <- list(estimates=mat, SE=semat, L=lmat, U=umat, fixed=fixed, sojourn=soj, sojournSE=sojse, sojournL=sojl, sojournU=soju)
        else
            res <- list(estimates=mat, SE=semat, L=lmat, U=umat, fixed=fixed)
        class(res) <- "msm.est"
    }
    res
}

ematrix.msm <- function(x, covariates="mean", ci=c("delta","normal","bootstrap","none"), cl=0.95, B=1000, cores=NULL)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    if (!x$emodel$misc) return(NULL)
    nst <- x$qmodel$nstates
    ni <- x$emodel$npars
    covlist <- msm.parse.covariates(x, covariates, x$ecmodel)
    nc <- length(covlist)
    if ((cl < 0) || (cl > 1)) stop("expected cl in [0,1]")
    se <- lse <- numeric(ni)
    logest <- x$Ematrices$logitbaseline
    if (is.null(x$EmatricesFixed)) x <- msm.form.output(x, "misc") # for back-compat with pre 1.4.1 model objects
    fixed <- x$EmatricesFixed$logitbaseline
    for (i in seq_len(nc)) {
        logest <- logest + x$Ematrices[[i+1]] * covlist[[i]]
        fixed <- fixed & x$EmatricesFixed[[i+1]]  # Only true if all coefficients contributing to the estimate are fixed.  Used in print functions
    }
    plabs <- x$emodel$imatrix
    plabs[x$emodel$imatrix==1] <- "p"
    diag(plabs)[rowSums(x$emodel$imatrix)>0] <- "pbase"
    hmodeltmp <- list(plabs = as.vector(t(plabs)), parstate = rep(1:nst, each=nst), parout = rep(1, length(plabs)))
    mat <- matrix(msm.mninvlogit.transform(as.vector(t(logest)), hmodeltmp),
                  nrow=nst, byrow=TRUE)
    ## true states with no misclassification or perfect misclassification 
    mat[cbind(which(x$hmodel$model==match("identity", .msm.HMODELS)),
              x$hmodel$pars[names(x$hmodel$pars)=="which"])] <- 1
    
    ci <- match.arg(ci)
    if (x$foundse && (ci!="none")) {
        if (ci == "delta") {
            ## Work out standard errors.
            ## Transformation for delta method is
            ## expit(sum covs)  / (1 + expit(sum(covs)))    or   1  /  (1  +  expit(sum(covs)))
            semat <- lmat <- umat <- matrix(0, nst, nst)
            p.se <- p.se.msm(x, covariates)
            ivector <- as.numeric(t(x$emodel$imatrix))
            if (any(p.se$lab %in% c("p","pbase"))){
                semat[ivector==1] <- p.se$se[p.se$lab=="p"]; semat <- t(semat)
                lmat[ivector==1] <- p.se$LCL[p.se$lab=="p"]; lmat <- t(lmat)
                umat[ivector==1] <- p.se$UCL[p.se$lab=="p"]; umat <- t(umat)
                diag(semat)[rowSums(x$emodel$imatrix)>0] <- p.se$se[p.se$lab=="pbase"]
                diag(lmat)[rowSums(x$emodel$imatrix)>0] <- p.se$LCL[p.se$lab=="pbase"]
                diag(umat)[rowSums(x$emodel$imatrix)>0] <- p.se$UCL[p.se$lab=="pbase"]
            }
            lmat[mat==1] <- umat[mat==1] <- 1
        }
        else if (ci=="normal") { 
            e.ci <- ematrix.normci.msm(x, covariates, cl, B)
            lmat <- e.ci[,,1]; umat <- e.ci[,,2]; semat <- e.ci[,,3]
        }
        else if (ci=="bootstrap") {
            e.ci <- ematrix.ci.msm(x, covariates, cl, B, cores)
            lmat <- e.ci[,,1]; umat <- e.ci[,,2]; semat <- e.ci[,,3]
        }
        dimnames(semat) <- dimnames(lmat) <- dimnames(umat) <- dimnames(x$qmodel$qmatrix)
    }
    else semat <- lmat <- umat <- NULL
    dimnames(mat) <-  dimnames(x$qmodel$qmatrix)
    if (ci=="none") res <- mat
    else {
        res <- list(estimates=mat, SE=semat, L=lmat, U=umat, fixed=fixed)
        class(res) <- "msm.est"
    }
    res
}

## Convert "covariates" argument supplied in output function (like
## qmatrix.msm) to a list of values, one per covariate effect
## (e.g. one per factor contrast).  Handle special arguments 0 and
## "mean".

msm.parse.covariates <- function(x, covariates, mod, consider.center=TRUE){
    nc <- mod$ncovs
    if (nc == 0){
        if (is.list(covariates) && (length(covariates) > 0))
            warning("Ignoring covariates - no covariates in this part of the model")
        return(list())
    }
    if (consider.center) {
        ## no adjustment needed: baseline is what we want
        if (!is.list(covariates) &&
            ((covariates==0 && !x$center) || (covariates=="mean" && x$center)))
            return(list())
    }
    if (!is.list(covariates)) {
        covlist <- list()
        if (covariates == 0) {
            for (i in 1:nc)
                covlist[[mod$covlabels[i]]] <- 0
        }
        else if (covariates == "mean") {
            for (i in 1:nc)
                covlist[[mod$covlabels[i]]] <- mod$covmeans[i]
        }
        else stop("covariates argument must be 0, \"mean\", or a list of values for each named covariate")
    }
    else {
        ## Check supplied list of covariate values, convert factors to numeric contrasts, expand interactions, set unknown values to zero.
        covlist <- factorcov2numeric.msm(covariates, x, mod)
    }
    if (x$center && consider.center)
        for (i in 1:nc)
            covlist[[mod$covlabels[i]]] <- covlist[[mod$covlabels[i]]] - mod$covmeans[i]            
    covlist
}

### Given a "covariates" argument of an extractor function containing
### factor covariates, convert the factor covariates to numeric
### contrasts.  For example, for a categorical covariate "smoke" with
### three levels "NON","CURRENT","EX", with baseline level "NON",
### convert list(smoke="CURRENT") to list(smokeCURRENT=1, smokeEX=0)
### Any unspecified covariate values are set to zero.
### Any unknown covariates are dropped with a warning.

factorcov2numeric.msm <- function(covariates, x, mod=NULL) {
    if (is.null(mod)) mod <- x$qcmodel
    covdata.mf <- x$data$mf[attr(x$data$mf,"covnames")]
    covnames.mm <- mod$covlabels
    covfactor <- sapply(covdata.mf, is.factor)
    covfactorlevels <- lapply(covdata.mf, levels)
    covnames <- names(covdata.mf)

    if (is.null(names(covariates))) {
        if (length(covariates)!=length(covnames)) stop("Expected covariate list of length ",length(covnames))
        names(covariates) <- covnames
    }
    all.covnames <- union(covnames.mm,covnames) # including both factor and contrast names
    miss.covs <- ! names(covariates) %in% all.covnames
    if (any(miss.covs)){
        plural <- if(sum(miss.covs)>1) "s" else ""
        warning("Covariate",plural," \"", paste(names(covariates)[which(!names(covariates)  %in% all.covnames)], collapse=", "), "\" unknown, ignoring")
    }
    cfac <- covariates[names(covariates) %in% covnames[which(covfactor)]]
    cnum <- covariates[! names(covariates) %in% covnames[which(covfactor)]]
    cfac.new <- list()
    for (i in seq_along(cfac)) {
        levs.i <- covfactorlevels[[names(cfac)[[i]]]]
        cfac.i <- rep(0, length(levs.i))
        if (! cfac[[i]] %in% levs.i) stop("Level \"", cfac[[i]], "\" of covariate ", names(cfac)[[i]], " unknown")
        cfac.i[match(cfac[[i]], levs.i)] <- 1
        names(cfac.i) <- paste(names(cfac)[[i]], levs.i, sep="")
        cfac.i <- as.list(cfac.i[-1])
        cfac.new <- c(cfac.new, cfac.i)
    }
    covlabels.noint <- covnames.mm[setdiff(seq_along(covnames.mm), grep(":", covnames.mm))]
    covs.out <- as.list(numeric(length(covlabels.noint)))
    names(covs.out) <- covlabels.noint
    covs.out[names(cnum)] <- cnum
    covs.out[names(cfac.new)] <- cfac.new
    covs.out <- expand.interactions.msm(covs.out, covnames.mm)
    covs.out
}

### Work out SE and CIs of misclassification probabilities for given covariate values.
### Know SEs of beta_rs, use delta method to calculate.
###  SE of p_rs = exp(beta_rs x) / (1 + sum(exp(beta_rs x)))  (non-baseline p) or  1 / (1 + sum(exp(beta_rs x))) (baseline p).
### To calculate symmetric CI for resulting p_rs, assume logit(p_rs) normal, and use SE(logit(p_rs)) calculated by delta method.

p.se.msm <- function(x, covariates)
{
    qmodel <- x$qmodel; emodel <- x$emodel; hmodel <- x$hmodel; qcmodel <- x$qcmodel; ecmodel <- x$ecmodel; paramdata <- x$paramdata
    nst <- qmodel$nstates
    inds <- (qmodel$npars + qcmodel$npars + 1) : (qmodel$npars + qcmodel$npars + sum(hmodel$npars) + hmodel$ncoveffs)
    ni <- emodel$npars
    covlist <- msm.parse.covariates(x, covariates, ecmodel)
    nc <- length(covlist) 
    coefs <- as.numeric(c(1, covlist))
    ppars <- hmodel$plabs %in% c("p","pbase")
    res <- data.frame(lab=hmodel$plabs[ppars])
    hmmallpars <- !(paramdata$plabs %in% c("qbase","qcov","initp","initp0","initpcov"))
    res$est <- msm.mninvlogit.transform(paramdata$params[paramdata$hmmpars], hmodel)[ppars]
    res$parstate <- hmodel$parstate[ppars]
    if (any(ppars)) res$se <- res$lse <- res$LCL <- res$UCL <- res$inds <- res$strs <- 0
    cur.i <- 1
    for (i in unique(res$parstate)) {
        nir <- sum(hmodel$parstate[hmodel$plabs=="p"] == i) # number of independent misc probs for current state
        if (nir > 0){  # might be perfect misclassification
            p.inds <- which(hmodel$plabs=="p" & hmodel$parstate==i) # indices into HMM parameter vector of logit baseline p for that state
            cov.inds <- sum(hmodel$npars) + (cur.i-1)*nc + seq(length=(nc*nir)) # indices into HMM parameter vector of corresp cov effects
            parinds <- numeric(); formstr <- character(nir)
            for (j in 1:nir) {
                formstr[j] <- expsum(1:((nc+1)*nir), coefs) # string of form exp(1*x1+beta1*x2*beta2*x3), exp(x4+beta*x5+...)
                parinds <- c(parinds, p.inds[j], cov.inds[(j-1)*nc + seq(length=nc)]) # indices into HMM par vector corresp to x1,x2,x3,...
            }
            sumstr <- paste(formstr, collapse = " + ") # "parinds" are
            formstr <- paste("1 / (1 + ", sumstr, ")")
            form <- as.formula(paste("~ ",formstr))
            lform <- as.formula(paste("~ log( (", formstr, ") / (1 - ", formstr, "))"))
            ests <- paramdata$params[hmmallpars][parinds]
            cov <- paramdata$covmat[hmmallpars,hmmallpars][parinds, parinds]
            res$se[res$parstate==i & res$lab=="pbase"] <- deltamethod(form, ests, cov)
            res$lse[res$parstate==i & res$lab=="pbase"] <- deltamethod(lform, ests, cov)
            res$strs[res$parstate==i & res$lab=="pbase"] <- paste(as.character(form),collapse="")
            res$inds[res$parstate==i & res$lab=="pbase"] <- paste(parinds,collapse=",")
            for (j in 1:nir){
                istr <- expsum(((j-1)*(nc+1)+1):(j*(nc+1)), coefs)
                formstr <- paste(istr, "/ (1 + ", sumstr, ")")
                form <- as.formula(paste("~ ",formstr))
                lform <- as.formula(paste("~ log( (", formstr, ") / (1 - ", formstr, "))"))
                res$se[res$parstate==i & res$lab=="p"][j] <- deltamethod(form, ests, cov)
                res$lse[res$parstate==i & res$lab=="p"][j] <- deltamethod(lform, ests, cov)
                res$strs[res$parstate==i & res$lab=="p"][j] <- paste(as.character(form), collapse="")
                res$inds[res$parstate==i & res$lab=="p"][j] <- paste(parinds,collapse=",")
            }
            cur.i <- cur.i + nir
    }
    }
    res$LCL <- plogis(qlogis(res$est) - qnorm(0.975)*res$lse)
    res$UCL <- plogis(qlogis(res$est) + qnorm(0.975)*res$lse)
    res
}

### Work out standard error of a ratio of intensities using delta method
### Uuugh.  What a fuss for one little number.

qratio.se.msm <- function(x, ind1, ind2, covariates, cl=0.95)
{
    nst <- x$qmodel$nstates
    ni <- x$qmodel$npars
    covlist <- msm.parse.covariates(x, covariates, x$qcmodel)    
    nc <- length(covlist)
    indmat <- t(x$qmodel$imatrix)
    indmat[indmat == 1] <- seq(length = x$qmodel$npars)
    indmat <- t(indmat) # matrix of indices of estimate vector
    inds <- seq(length = x$qmodel$npars+x$qcmodel$npars) # identifiers for q and beta parameters
    coefs <- as.numeric(c(1, unlist(covlist)))
    parinds <- numeric()
    indmatrow.n <- indmat[ind1[1],-ind1[1]]
    nir.n <- sum(indmatrow.n > 0)
    indmatrow.d <- indmat[ind2[1],-ind2[1]]
    nir.d <- sum(indmatrow.d > 0)
    formstr.n <- character(nir.n)
    formstr.d <- character(nir.d)

    if (ind1[1]!=ind1[2] && ind2[1]!=ind2[2]) { # both intensities are off-diagonal
        parinds <- c(inds[indmat[ind1[1],ind1[2]] - 1 + seq(1, (nc * ni + 1), ni)],
                     inds[indmat[ind2[1],ind2[2]] - 1 + seq(1, (nc * ni + 1), ni)])
        parinds2 <- sort(unique(parinds))
        xinds <- rank(parinds2)[match(parinds, parinds2)]
        formstr.n <- expsum(xinds[1:(nc+1)], coefs)
        formstr.d <- expsum(xinds[1:(nc+1) + nc+1] , coefs)
    }

    else if (ind1[1]!=ind1[2] && ind2[1]==ind2[2]) { # numerator off-diagonal, denom diagonal
        parinds <- inds[indmat[ind1[1],ind1[2]] - 1 + seq(1, (nc * ni + 1), ni)]
        cur.i <- min(indmatrow.d[indmatrow.d>0])
        for (j in 1:nir.d)
            parinds <- c(parinds, inds[cur.i - 1 + seq(j, (nc * ni + j), ni)])
        parinds2 <- sort(unique(parinds))
        xinds <- rank(parinds2)[match(parinds, parinds2)]
        formstr.n <- expsum(xinds[1:(nc+1)], coefs)
        for (j in 1:nir.d)
            formstr.d[j] <- expsum(xinds[1:(nc+1) + j*(nc+1)], coefs)
    }

    else if (ind1[1]==ind1[2] && ind2[1]!=ind2[2]) { # numerator diagonal, denom off-diagonal
        cur.i <- min(indmatrow.n[indmatrow.n>0])
        for (j in 1:nir.n)
            parinds <- c(parinds, inds[cur.i - 1 + seq(j, (nc * ni + j), ni)])
        parinds <- c(parinds, inds[indmat[ind2[1],ind2[2]] - 1 + seq(1, (nc * ni + 1), ni)])
        parinds2 <- sort(unique(parinds))
        xinds <- rank(parinds2)[match(parinds, parinds2)]
        for (j in 1:nir.n)
            formstr.n[j] <- expsum(xinds[1:(nc+1) + (j-1)*(nc+1)], coefs)
        formstr.d <- expsum(xinds[nir.n*(nc+1) + 1:(nc+1)], coefs)
    }

    else if (ind1[1]==ind1[2] && ind2[1]==ind2[2]) { # both intensities diagonal
        cur.i <- min(indmatrow.n[indmatrow.n>0])
        for (j in 1:nir.n)
            parinds <- c(parinds, inds[cur.i - 1 + seq(j, (nc * ni + j), ni)])
        cur.i <- min(indmatrow.d[indmatrow.d>0])
        for (j in 1:nir.d)
            parinds <- c(parinds, inds[cur.i - 1 + seq(j, (nc * ni + j), ni)])
        parinds2 <- sort(unique(parinds))
        xinds <- rank(parinds2)[match(parinds, parinds2)]
        for (j in 1:nir.n)
            formstr.n[j] <- expsum(xinds[1:(nc+1) + (j-1)*(nc+1)], coefs)
        for (j in 1:nir.d)
            formstr.d[j] <- expsum(xinds[nir.n*(nc+1) + 1:(nc+1) + (j-1)*(nc+1)], coefs)
    }

    num <- paste(formstr.n, collapse = " + ")
    denom <- paste(formstr.d, collapse = " + ")
    form <- as.formula(paste("~", "(", num, ") / (", denom, ")"))
    lform <- as.formula(paste("~ ", "log (", num, ") - log (", denom, ")"))
    ests <- x$estimates[parinds2]
    cov <- x$covmat[parinds2,parinds2]
    se <- deltamethod(form, ests, cov)
    lse <- deltamethod(lform, ests, cov)
    list(se=se, lse=lse)
}

### Work out standard errors of diagonal entries of intensity matrix, or sojourn times, using delta method

qmatrix.diagse.msm <- function(x, covlist, sojourn, ni, ivector, nc)
{
    nst <- x$qmodel$nstates
    diagse <- diaglse <- sojse <- sojlse <- numeric(nst)
    indmat <- matrix(ivector, nst, nst)
    indmat[indmat==1] <- seq(length = ni)
    indmat <- t(indmat) # matrix of indices of estimate vector
    inds <- seq(length = ni + ni*nc)
    cur.i <- 1
    coefs <- as.numeric(c(1, unlist(covlist)))
    for (i in 1:nst){
        ## Transformation for delta method is
        ## exp(x1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... ) +
        ##  exp(x4 + x5 (cov1 - covmean1) + x6 (cov2 - covmean2) + ... ) +   (or expit(...))
        nir <- sum(indmat[i,-i] > 0) # number of intens/misc for current state
        if (nir > 0) {
            qf <- expsum.formstr(nir, inds, cur.i, ni, nc, coefs)
            form <- as.formula(paste("~", paste(qf$formstr, collapse = " + ")))
            lform <- as.formula(paste("~ log (", paste(qf$formstr, collapse = " + "), ")"))
            ests <- x$estimates[qf$parinds2]
            cov <- x$covmat[qf$parinds2, qf$parinds2]
            diagse[i] <- deltamethod(form, ests, cov)
            diaglse[i] <- deltamethod(lform, ests, cov)
            if (sojourn){
                ## Mean sojourn times are -1 / diagonal entries of q matrix. Calculate their SEs and CIs.
                form <- as.formula(paste("~ 1 / (", paste(qf$formstr, collapse = " + "), ")"))
                lform <- as.formula(paste("~ log ( 1 / (", paste(qf$formstr, collapse = " + "), ")", ")"))
                sojse[i] <- deltamethod(form, ests, cov)
                sojlse[i] <- deltamethod(lform, ests, cov)
            }
            cur.i <- cur.i + nir
        }
        else diagse[i] <- 0
    }
    list(diagse=diagse, diaglse=diaglse, sojse=sojse, sojlse=sojlse)
}

### Make a list of covariate lists to supply to pmatrix.piecewise.msm for models with "pci" time-dependent intensities.
### One for each time period, with time constant covariates replicated.
### For use in model assessment functions
### Returns factor covariates as contrasts, not factor levels.

msm.fill.pci.covs <- function(x, covariates="mean"){
    nc <- x$qcmodel$ncovs
    ## indices of covariates representing time periods
    ti <- grep("timeperiod\\[.+\\)", x$qcmodel$covlabels)
    ni <- setdiff(1:nc, ti) # indices of other covariates
    covlist <- msm.parse.covariates(x, covariates, x$qcmodel, consider.center=FALSE)
    for (i in names(covariates))
        if (length(grep("^timeperiod",i))==0) {
            if (i %in% union(attr(x$data$mf, "covnames"),x$qcmodel$covlabels))
                covlist[[i]] <- covariates[[i]]
            else warning("Covariate ",i," unknown")
        }
    for (i in ti) covlist[[i]] <- 0
    ## set contrasts for each successive time period to 1
    ncut <- length(x$pci)
    covlistlist <- vector(ncut+1, mode="list")
    names(covlistlist) <- levels(x$data$mf$timeperiod)
    covlistlist[[1]] <- covlist
    for (i in seq(length=ncut)){
        covlistlist[[i+1]] <- covlist
        covlistlist[[i+1]][[ti[i]]] <- 1
    }
    covlistlist
}

printold.msm <- function(x, ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")

    if (!attr(x,"fixed")) {
        cat ("Maximum likelihood estimates: \n")
        covmessage <-
            if (x$qcmodel$ncovs == 0) ""
            else paste("with covariates set to", (if (x$center) "their means" else "0"))
        for (i in c("baseline", x$qcmodel$covlabels)) {
            title <-
                if (i == "baseline") paste("Transition intensity matrix",covmessage,"\n")
                else paste("Log-linear effects of", i, "\n")
            cat (title, "\n")
            print.ci(x$Qmatrices[[i]], x$QmatricesL[[i]], x$QmatricesU[[i]])
            cat("\n")
        }
        if (x$emodel$misc) {
            misccovmessage <-
                if (x$ecmodel$ncovs == 0) ""
                else paste("with covariates set to", (if (x$center) "their means" else "0"))
            for (i in c("baseline", x$ecmodel$covlabels)) {
                title <-
                    if (i == "baseline") paste("Misclassification matrix",misccovmessage,"\n")
                    else paste("Effects of", i, "on log (P(state r)/P(baseline state))\n")
                cat (title, "\n")
                print.ci(x$Ematrices[[i]], x$EmatricesL[[i]], x$EmatricesU[[i]])
                cat("\n")
            }
            if (any(x$paramdata$plabs[x$paramdata$optpars] == "initp")) {
                cat("Initial state occupancy probabilities\n\n")
                print(x$hmodel$initprobs)
                cat("\n")
                if (any(x$hmodel$nicovs > 0)) {
                    cat("Covariates on logit initial state probabilities\n")
                    print(x$hmodel$icoveffect)
                }
                cat("\n")
            }
        }
        else if (x$hmodel$hidden && is.null(x$qmodel$phase.states)) {
            print(x$hmodel); cat("\n")
        }
    }
    cat ("-2 * log-likelihood: ", x$minus2loglik, "\n")
#    cat("[Note: a cleaner summary is available from \"printnew.msm\",\n which will be the default in future versions.]\n")

}

### Convert three-transition-matrices (estimate,lower,upper) format to three-columns format

mattotrans <- function(x, matrix, lower, upper, fixed, keep.diag=FALSE, intmisc="intens"){
    imat <- if (intmisc=="intens") x$qmodel$imatrix else x$emodel$imatrix
    if (keep.diag) diag(imat) <- as.numeric(rowSums(imat) > 0)
    keep <- which(t(imat)==1, arr.ind=TRUE)
    keep <- keep[,2:1,drop=FALSE]  # order by row(from-state), not column(to-state)
    fromlabs <- rownames(imat)[keep[,1]]
    tolabs <- colnames(imat)[keep[,2]]
    res <- matrix(nrow=sum(imat), ncol=4)
    rnames <- if (intmisc=="intens") paste(fromlabs, "-", tolabs) else paste("Obs", tolabs, "|", fromlabs)
    dimnames(res) <- list(rnames, c("Estimate", "L", "U","Fixed"))
    res[,1] <- matrix[keep]
    res[,2] <- lower[keep]
    res[,3] <- upper[keep]
    res[,4] <- fixed[keep]
    res
}

### Format transition intensities and their covariate effects in one tidy matrix

msm.form.qoutput <- function(x, covariates="mean", cl=0.95, digits=4, ...){
    qbase <- qmatrix.msm(x, covariates=covariates, cl=cl)
    if (is.null(x$QmatricesFixed)) x <- msm.form.output(x, "intens") # for back-compat with pre 1.4.1 model objects
    y <- mattotrans(x, qbase$estimates, qbase$L, qbase$U, qbase$fixed, keep.diag=TRUE)
    ret <- data.frame(base=y)
    fres <- matrix("", nrow=nrow(y), ncol=x$qcmodel$ncovs+1)
    colnames(fres) <- c("Baseline", x$qcmodel$covlabels)
    rownames(fres) <- rownames(y)
    fres[,1] <- format.ci(y[,1],y[,2],y[,3],y[,4],digits=digits,...)
    im <- t(x$qmodel$imatrix); diag(im) <- -colSums(im); nd <- which(im[im!=0]==1)
    for (i in seq(length=x$qcmodel$ncovs)){
        nm <- x$qcmodel$covlabels[[i]]
        hrs <- mattotrans(x, x$Qmatrices[[nm]], x$QmatricesL[[nm]], x$QmatricesU[[nm]], x$QmatricesFixed[[nm]], keep.diag=FALSE)
        hrs[,1:3] <- exp(hrs[,1:3])
        ret[nm] <- matrix(ncol=3, nrow=nrow(ret), dimnames=list(NULL,colnames(hrs)[1:3]))
        ret[nd,nm] <- hrs[,1:3,drop=FALSE]
        fres[nd,1+i] <- format.ci(hrs[,1], hrs[,2], hrs[,3], hrs[,4], digits=digits, ...)
    }
    attr(ret, "formatted") <- fres # as strings with formatted CIs instead of numbers
    ret
}

### Format misclassification intensities and their covariate effects in one tidy matrix

msm.form.eoutput <- function(x, covariates="mean", cl=0.95, digits=4, ...){
    ebase <- ematrix.msm(x, covariates=covariates, cl=cl)
    if (is.null(x$EmatricesFixed)) x <- msm.form.output(x, "misc") # for back-compat with pre 1.4.1 model objects
    y <- mattotrans(x, ebase$estimates, ebase$L, ebase$U, ebase$fixed, keep.diag=TRUE, intmisc="misc")
    rete <- data.frame(base=y)
    frese <- matrix("", nrow=nrow(y), ncol=x$ecmodel$ncovs+1)
    colnames(frese) <- c("Baseline", x$ecmodel$covlabels)
    rownames(frese) <- rownames(y)
    frese[,1] <- format.ci(y[,1],y[,2],y[,3],y[,4],digits=digits,...)
    im <- t(x$emodel$imatrix); diag(im) <- -colSums(im); nd <- which(im[im!=0]==1)
    for (i in seq(length=x$ecmodel$ncovs)){
        nm <- x$ecmodel$covlabels[[i]]
        ors <- mattotrans(x, x$Ematrices[[nm]], x$EmatricesL[[nm]], x$EmatricesU[[nm]], x$EmatricesFixed[[nm]], keep.diag=FALSE, intmisc="misc")
        ors[,1:3] <- exp(ors[,1:3])
        rete[nm] <- matrix(ncol=3, nrow=nrow(rete), dimnames=list(NULL,colnames(ors)[1:3]))
        rete[nd,nm] <- ors[,1:3]
        frese[nd,1+i] <- format.ci(ors[,1], ors[,2], ors[,3], ors[,4], digits=digits,...)
    }
    attr(rete, "formatted") <- frese # as strings with formatted CIs instead of numbers
    rete
}

## New more helpful and tidier print output

print.msm <- function(x, covariates=NULL, digits=4, ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    ret <- NULL
    if (!x$foundse & !attr(x, "fixed")) {
        cat("Optimisation probably not converged to the maximum likelihood.\noptim() reported convergence but estimated Hessian not positive-definite.\n")
    }
    else {
        if (is.null(x$cl)) {
            cl <- 0.95
            warning("Found msm object saved before version 1.3. Models will need to be refitted under the newer version for output functions to work")
        } else cl <- x$cl
        if (!attr(x,"fixed")) {
            if (is.null(covariates)) {
                covvalue <- (if (x$center) "their means" else "0")
                covariates <- if (x$center) "mean" else 0
            } else covvalue <- "the values supplied in \"covariates\""
            covmessage <-
                if (attr(x$data$mf, "ncovs")==0) ""
                else paste0("\nBaselines are with covariates set to ", covvalue)
            cat("Maximum likelihood estimates", covmessage, "\n", sep="")
            if (x$qcmodel$ncovs> 0)
                hrmessage <- paste0(" with hazard ratios for each covariate")
            else hrmessage <- ""
            q.header <- paste0("Transition intensities", hrmessage, "\n")
            ret <- msm.form.qoutput(x, covariates, cl=cl, digits=digits, ...)
            fres <- attr(ret, "formatted")
            cat("\n"); cat(q.header)
            print(fres, quote=FALSE)
            if (x$emodel$misc) {
                ormessage <- if (x$ecmodel$ncovs>0) paste0(" with odds ratios for each covariate") else ""
                e.header <- paste0("Misclassification probabilities", ormessage, "\n")
                rete <- msm.form.eoutput(x, covariates, cl=cl, digits=digits, ...)
                frese <- attr(rete, "formatted")
                cat("\n"); cat(e.header)
                print(frese, quote=FALSE)
                if (any(x$paramdata$plabs[x$paramdata$optpars] == "initp")) {
                    i.header <- paste0("Initial state occupancy probabilities\n")
                    cat("\n"); cat(i.header)
                    print(x$hmodel$initprobs)
                    if (any(x$hmodel$nicovs > 0)) {
                        ic.header <- "Covariates on logit initial state probabilities\n"
                        cat("\n"); cat(ic.header)
                        print(x$hmodel$icoveffect)
                    }
                }
            }
            else if (x$hmodel$hidden && (is.null(x$hmodel$phase.only) || !x$hmodel$phase.only)){
                cat("\n")
                print(x$hmodel)
            }
            if (!is.null(x$qmodel$phase.states)) {
                cat("\nPhase-type model\n")
                print(phasemeans.msm(x))
            }
        }
    }
    cat ("\n-2 * log-likelihood: ", x$minus2loglik, "\n")
    cat("[Note, to obtain old print format, use \"printold.msm\"]\n",sep="")
    invisible(ret)
}

printnew.msm <- print.msm

summary.msm <- function(object, # fitted model
                        hazard.scale = 1,
                        ...
                        )
{
    if (!inherits(object, "msm")) stop("expected object to be a msm model")
    prevalences <- prevalence.msm(object, ...)
    if (object$qcmodel$ncovs > 0) {
        if (missing (hazard.scale))
            hazard.scale <- rep(1, object$qcmodel$ncovs)
        hazard <- hazard.msm(object)
    }
    else {hazard <- hazard.scale <- NULL}
    ret <- list(prevalences=prevalences,
                hazard=hazard,
                hazard.scale=hazard.scale)
    class(ret) <- "summary.msm"
    ret
}

print.summary.msm <- function(x,...)
{
    if (!is.null(x$prevalences)) {
        cat("\nObserved numbers of individuals occupying states at each time\n\n")
        print(x$prevalences$Observed)
        cat("\nExpected numbers of individuals occupying states at each time\n\n")
        print(x$prevalences$Expected)
        cat("\nObserved prevalences of states (percentages of population at risk)\n\n")
        print(x$prevalences$"Observed percentages")
        cat("\nExpected prevalences of states (percentages of population at risk)\n\n")
        print(x$prevalences$"Expected percentages")
    }
    i <- 1
    for (cov in names(x$hazard)) {
        cat ("\nTransition hazard ratios corresponding to covariate effects\n\n" )
        cat (cov, " ( unit of",x$hazard.scale[i],")\n")
        print(round(x$hazard[[cov]], 2))
        i <- i+1
    }
    invisible()
}

### Estimated survival probability from each state

plot.msm <- function(x, from=NULL, to=NULL, range=NULL, covariates="mean", legend.pos=NULL, xlab="Time", ylab="Fitted survival probability", lwd=1,...)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    if (is.null(from))
        from <- transient.msm(x)
    else {
        if (!is.numeric(from)) stop("from must be numeric")
        if (any (! (from %in% 1:x$qmodel$nstates ) ) )
            stop("from must be a vector of states in 1, ..., ", x$qmodel$nstates)
    }
    if (is.null(to)){
        if (length(absorbing.msm(x))==0)
            stop("\"to\" not specified, and no absorbing state in the model")
        to <- max(absorbing.msm(x))
    }
    else {
        if (!is.numeric(to)) stop("to must be numeric")
        if (! (to %in% absorbing.msm(x) ) ) stop("to must be an absorbing state")
    }
    if (is.null(range))
        rg <- range(model.extract(x$data$mf, "time"))
    else {
        if (!is.numeric(range) || length(range)!= 2) stop("range must be a numeric vector of two elements")
        rg <- range
    }
    timediff <- (rg[2] - rg[1]) / 50
    times <- seq(rg[1], rg[2], timediff)
    pr <- numeric()
    cols <- rainbow(length(from))
    for (t in times)
        pr <- c(pr, pmatrix.msm(x, t, times[1], covariates)[from[1], to])
    plot(times, 1 - pr, type="l", xlab=xlab, ylab=ylab, lwd=lwd,
         ylim=c(0,1), lty = 1, col=cols[1],...)
    lt <- 2
    for (st in from[-1]){
        pr <- numeric()
        for (t in times)
            pr <- c(pr, pmatrix.msm(x, t, times[1], covariates)[st, to])
        lines(times, 1 - pr, type="l", lty = lt, lwd=lwd, col=cols[lt],...)
        lt <- lt+1
    }
    if (!is.numeric(legend.pos) || length(legend.pos) != 2)
        legend.pos <- c(max(times) - 15*timediff, 1)
    legend(legend.pos[1], legend.pos[2], legend=paste("From state",from), lty = seq(lt-1), col=cols, lwd=lwd)
    invisible()
}

### Plot KM estimate of time to first occurrence of each state

plotprog.msm <- function(formula, subject, data, legend.pos=NULL, xlab="Time", ylab="1 - incidence probability", lwd=1, xlim=NULL,
                         mark.time=TRUE, ...) {
    data <- na.omit(data)
    mf <- model.frame(formula, data=data)
    state <- mf[,1]
    time <- mf[,2]
    if (!is.null(data))
        subject <- eval(substitute(subject), as.list(data), parent.frame())
    subject <- match(subject, unique(subject))
    rg <- range(time)
    if (is.null(xlim)) xlim=rg
    plot(0, xlim=xlim, ylim=c(0,1), type="n", xlab=xlab, ylab=ylab, ...)
    states <- sort(unique(state))[-1]
    cols <- rainbow(length(states))
    for (i in states) {
        dat <- cbind(subject, time, state)
        st <- as.data.frame(
                            do.call("rbind", by(dat, subject, function(x)
                                            {
                                                c(anystate = if(any(x[,"state"]>=i)) 1 else 0,
                                                  mintime = if(any(x[,"state"]>=i)) min(x[x[,"state"] >= i, "time"]) else max(x[,"time"]))
                                            }
                                                ))) # slow
        lines(survfit(Surv(st$mintime,st$anystate) ~ 1),
              col=cols[i-1], lty=i-1, lwd=lwd, mark.time=mark.time, ...)
    }
    timediff <- (rg[2] - rg[1]) / 50
    if (!is.numeric(legend.pos) || length(legend.pos) != 2)
        legend.pos <- c(max(time) - 25*timediff, 1)
    legend(legend.pos[1], legend.pos[2], lty=states-1, lwd=lwd, col=cols,
           legend=paste("To state", states, c(rep("or greater", length(states)-1), "")))
    invisible()
}


### Likelihood surface plots

surface.msm <- function(x, params=c(1,2), np=10, type=c("contour","filled.contour","persp","image"),
                        point=NULL, xrange=NULL, yrange=NULL,...)
{
    type <- match.arg(type)
    if (is.null(point))
        point <- x$paramdata$opt$par
    se <- sqrt(diag(x$covmat[x$paramdata$optpars,x$paramdata$optpars]))
    i1 <- params[1]; i2 <- params[2]
    if (is.null(xrange)) {
        pmin <- point[i1] - 2*se[i1]
        pmax <- point[i1] + 2*se[i1]
        p1 <- seq(pmin, pmax, length=np)
    }
    else p1 <- seq(xrange[1], xrange[2], length=np)
    if (is.null(yrange)){
        pmin <- point[i2] - 2*se[i2]
        pmax <- point[i2] + 2*se[i2]
        p2 <- seq(pmin, pmax, length=np)
    }
    else p2 <- seq(yrange[1], yrange[2], length=np)

    z <- matrix(nrow=np, ncol=np)
    for (i in 1:np) {
        for (j in 1:np) {
            point[i1] <- p1[i]; point[i2] <- p2[j]
            z[i,j] <- -0.5*Ccall.msm(point, "lik", expand.data(x), x$qmodel, x$qcmodel, x$cmodel, x$hmodel, x$paramdata)
        }
    }

    switch(type,
           contour = contour(p1, p2, z, ...),
           filled.contour = filled.contour(p1, p2, z, ...),
           image = image(p1, p2, z, ...),
           persp = persp(p1, p2, z, zlab="Log-likelihood",...)
           )
    invisible()
}

contour.msm <- function(x, ...)
{
    surface.msm(x, type="contour",...)
}

persp.msm <- function(x, ...)
{
    surface.msm(x, type="persp",...)
}

image.msm <- function(x, ...)
{
    surface.msm(x, type="image",...)
}

### Given a "covariates" argument of an extractor function containing
### covariates which interact in the model, form the corresponding
### "covariates" argument with the interactions expanded.

expand.interactions.msm <- function(covariates, covlabels){
    cn.nointer <- names(covariates)
    elist <- strsplit(covlabels, ":")
    elist <- lapply(elist, function(x)covariates[x])
    elist <- lapply(elist, function(x)prod(unlist(x)))
    names(elist) <- covlabels
    elist
}

print.msm.est <- function(x, digits=NULL, ...)
{
    if (is.list(x))
        print.ci(x$estimates, x$L, x$U, x$fixed, digits=digits)
    else print(unclass(x), digits=digits)
}

print.msm.est.cols <- function(x, digits=NULL, diag=TRUE, ...)
{
    inc <- if (diag) (x$estimates>0 | x$estimates<0) else (x$estimates>0)
    res <- cbind(x$estimates[inc], x$L[inc], x$U[inc])
    rn <- rownames(x$estimates)[row(inc)[inc]]
    cn <- colnames(x$estimates)[col(inc)[inc]]
    rownames(res) <- paste(rn, cn, sep="-")
    colnames(res) <- c("Estimate", "LCL", "UCL")
    res
}

"[.msm.est" <- function(x, i, j, drop=FALSE){
    Narg <- nargs() - (!missing(drop)) # number of args including x, excluding drop
    if ((missing(i) && missing(j)))
        res <- x
    else if (!is.list(x))
        res <- unclass(x)[i,j]
    else {
        if (missing(j) && (Narg==2))
            stop("Two dimensions must be supplied, found only one")
        if ("SE" %in% names(x)) {
            x <- array(unlist(x), dim=c(dim(x[[1]]),4))
            dimnames(x) <- list(rownames(x[[1]]), colnames(x[[1]]), c("estimate","SE","lower","upper"))
        }
        else {
            x <- array(unlist(x), dim=c(dim(x[[1]]),3))
            dimnames(x) <- list(rownames(x[[1]]), colnames(x[[1]]), c("estimate","lower","upper"))
        }
        res <- x[i,j,]
    }
    res
}

format.ci <- function(x, l, u, noci=NULL, digits=NULL, ...)
{
    if (is.null(noci)) noci <- rep(FALSE, length(x))
    if (is.null(digits)) digits <- 4
    ## note format() aligns nicely on point, unlike formatC
    est <- format(x, digits=digits, ...)
    res <- est
    if (!is.null(l)) {
        low <- format(l[!noci], digits=digits, ...)
        upp <- format(u[!noci], digits=digits, ...)
        res[!noci] <- paste(res[!noci], " (", low, ",", upp, ")", sep="")
        res[x==0] <- 0
    }
    else res <- est
    dim(res) <- dim(x)
    dimnames(res) <- dimnames(x)
    names(res) <- names(x)
    res
}

print.ci <- function(x, l, u, fixed=NULL, digits=NULL){
    res <- format.ci(x, l, u, fixed, digits)
    print(res, quote=FALSE)
}

### Work out CIs of initial state occupancy probabilities using normal simulation method
initp.ci.msm <- function(paramdata, cl=0.95){
    p <- paramdata
    Sig <- p$covmat[p$plabs%in%c("initp"),p$plabs%in%c("initp"),drop=FALSE]
    mu <- p$params[p$plabs%in%c("initp")]
    rr <- rmvnorm(10000, mu, Sig)
    initp.rep <- exp(rr) / (1 + rowSums(exp(rr)))
    p$ci[p$plabs=="initp",] <- t(apply(initp.rep, 2, quantile, c(0.5*(1-cl), 1 - 0.5*(1-cl)), na.rm=TRUE))
    initp.rep <- 1 / (1 + rowSums(exp(rr)))
    p$ci[p$plabs=="initpbase",] <- quantile(initp.rep, c(0.5*(1-cl), 1 - 0.5*(1-cl)))
    p$ci[p$plabs=="initp0"] <- 0
    p
}

## Form a string,  exp ( x1 1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... )
## to be made into a formula for deltamethod.
## Also return indices of estimates and covariance matrix output from optim() to use
## to inform deltamethod

expsum.formstr <- function(nir, inds, cur.i, ni, nc, coefs)
{
    formstr <- character(nir)
    parinds <- numeric()
    for (j in (cur.i : (cur.i + nir - 1))) {
        indj <- seq(j, (nc * ni + j), ni)
        parinds <- c(parinds, inds[indj])  # first 1 5 7, then 2 5 7
    }
    ## e.g. parinds = 1 5 7 2 5 7 becomes xinds = 1 3 4 2 3 4
    parinds2 <- sort(unique(parinds))
    xinds <- rank(parinds2)[match(parinds, parinds2)]
    for (j in 1:nir)
        formstr[j] <- expsum(xinds[1:(nc+1) + (j-1)*(nc+1)], coefs)
    list(formstr=formstr, parinds=parinds, parinds2=parinds2)
}

## Form a string,  exp ( x1 1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... )
## to be made into a formula for deltamethod

expsum <- function(inds, coefs)
{
    xseq <- paste("x", inds, sep="")
    inprod <- paste(paste(coefs, xseq, sep="*"), collapse=" + ")
    paste("exp(", inprod, ")", sep="")
}

lsum <- function(inds, coefs)
{
    xseq <- paste("x", inds, sep="")
    paste(paste(coefs, xseq, sep="*"), collapse=" + ")
}


### Extract a ratio of transition intensities at given covariate values

qratio.msm <- function(x, ind1, ind2,
                       covariates = "mean",
                       ci=c("delta","normal","bootstrap","none"),
                       cl=0.95, B=1000, cores=NULL)
{
    q <- qmatrix.msm(x, covariates, ci="none")
    if (!is.numeric(ind1) || length(ind1) != 2 || !is.numeric(ind2) || length(ind2) != 2)
        stop("ind1 and ind2 must be numeric vectors of length 2")
    if (any (! (ind1 %in% 1 : x$qmodel$nstates))  |  any (! (ind2 %in% 1 : x$qmodel$nstates) ) )
        stop("ind1 and ind2 must be pairs of states in 1, ..., ", x$qmodel$nstates)
    if ((cl < 0) || (cl > 1)) stop("expected cl in [0,1]")
    if (q[ind2[1], ind2[2]] == 0)
        stop (paste("Denominator q[",ind2[1],",",ind2[2],"", "] is zero\n", sep=""))
    else if (q[ind1[1], ind1[2]] ==  0) {
        warning(paste ("Numerator q[",ind1[1],",",ind1[2],"", "] is zero\n", sep=""))
        estimate <- se <- 0
    }
    else {
        estimate <- q[ind1[1], ind1[2]]  /  q[ind2[1], ind2[2]]
        ci <- match.arg(ci)
        if (x$foundse && (ci != "none")) {
            if (ci == "delta") {
                se <- qratio.se.msm(x, ind1, ind2, covariates, cl)$se
                lse <- qratio.se.msm(x, ind1, ind2, covariates, cl)$lse
                L <- exp ( log(abs(estimate)) - sign(estimate)*qnorm(1 - 0.5*(1 - cl)) * lse ) * sign(estimate)
                U <- exp ( log(abs(estimate)) + sign(estimate)*qnorm(1 - 0.5*(1 - cl)) * lse ) * sign(estimate)
            }
            else if (ci=="normal") {
                q.ci <- qratio.normci.msm(x, ind1, ind2, covariates, cl, B)
                L <- q.ci[1]; U <- q.ci[2]; se=q.ci[3]
            }
            else if (ci=="bootstrap") {
                q.ci <- qratio.ci.msm(x, ind1, ind2, covariates, cl, B, cores)
                L <- q.ci[1]; U <- q.ci[2]; se=q.ci[3]
            }
        }
        else {se <- L <- U <- NULL}
    }
    c(estimate=estimate, se=se, L=L, U=U)
}


### Extract the transition probability matrix at given covariate values

pmatrix.msm <- function(x=NULL, # fitted msm model
                        t = 1, # time interval
                        t1 = 0, # start time for pci models
                        covariates = "mean",  # covariate values to calculate transition matrix for
                        ci=c("none","normal","bootstrap"), # calculate a confidence interval
                                        # using either simulation from asymptotic normal dist of MLEs, or bootstrap
                        cl = 0.95, # width of symmetric confidence interval
                        B = 1000, # number of bootstrap replicates or normal simulations
                        cores=NULL,
                        qmatrix=NULL,
                        ...
                        )
{
    if (!is.numeric(t) || (t < 0)) stop("t must be a positive number")
    if (!is.null(x)) { 
        if (!inherits(x, "msm")) stop("expected x to be a msm model")
        if (is.null(x$pci)) {
            q <- qmatrix.msm(x, covariates, ci="none")
            p <- MatrixExp(q, t, ...)
            colnames(p) <- rownames(p) <- rownames(q)
            ci <- match.arg(ci)
            p.ci <- switch(ci,
                           bootstrap = pmatrix.ci.msm(x=x, t=t, t1=t1, covariates=covariates, cl=cl, B=B, cores=cores),
                           normal = pmatrix.normci.msm(x=x, t=t, t1=t1, covariates=covariates, cl=cl, B=B),
                           none = NULL)
            res <- if (ci=="none") p else list(estimates = p, L=p.ci[,,1], U=p.ci[,,2])
        }
        else {
            piecewise.covariates <- msm.fill.pci.covs(x, covariates)
            res <- pmatrix.piecewise.msm(x, t1, t1 + t, x$pci, piecewise.covariates, ci, cl, B, ...)
        }
    }
    else if (!is.null(qmatrix)){
        res <- MatrixExp(qmatrix, t, ...)
    }
    else stop("Neither a fitted model nor a qmatrix supplied")
    class(res) <- "msm.est"
    res
}

### Extract the transition probability matrix at given covariate values - where the Q matrix is piecewise-constant

pmatrix.piecewise.msm <- function(x=NULL, # fitted msm model
                                  t1, # start time
                                  t2, # stop time
                                  times,  # vector of cut points
                                  covariates, # list of lists of covariates, for (, times1], (times1, times2], ...
                                        # of length one greater than times
                                  ci=c("none","normal","bootstrap"),
                                  cl = 0.95, # width of symmetric confidence interval
                                  B = 1000, # number of bootstrap replicates or normal simulations
                                  cores=NULL,
                                  qlist=NULL,
                                  ... # arguments to pass to MatrixExp
                                  )
{
    if (!is.null(x)) { if (!inherits(x, "msm")) stop("expected x to be a msm model") }
    if (is.null(x) && is.null(qlist))  stop("Neither a fitted model nor a list of Q matrices have been supplied")
    x$pci <- NULL # to avoid infinite recursion when calling pmatrix.msm
    ## Input checks
    if (t2 < t1) stop("Stop time t2 should be greater than or equal to start time t1")
    if (!is.numeric(times) || is.unsorted(times)) stop("times should be a vector of numbers in increasing order")
    if (length(covariates) != length(times) + 1)
        stop("Number of covariate lists must be one greater than the number of cut points")
    if (length(times)==0)
        return(pmatrix.msm(x=x, t=t2-t1, t1=t1, covariates=covariates[[1]], ci=ci, cl=cl, B=B, qmatrix=qlist[[1]], ...))
    ## Locate which intervals t1 and t2 fall in, as indices ind1, ind2 into "times".
    if (t1 <= times[1]) ind1 <- 1
    else if (length(times)==1) ind1 <- 2
    else {
        for (i in 2:length(times))
            if ((t1 > times[i-1]) && (t1 <= times[i]))
            {ind1 <- i; break}
        if (t1 > times[i]) ind1 <- i+1
    }
    if (t2 <= times[1]) ind2 <- 1
    else if (length(times)==1) ind2 <- 2
    else {
        for (i in 2:length(times))
            if ((t2 > times[i-1]) && (t2 <= times[i]))
            {ind2 <- i; break}
        if (t2 > times[i]) ind2 <- i+1
    }

    ## Calculate accumulated pmatrix
    ## Three cases: ind1, ind2 in the same interval
    if (ind1 == ind2) {
        P <- pmatrix.msm(x=x, t = t2 - t1, covariates=covariates[[ind1]], qmatrix=qlist[[ind1]], ...)
    }
    ## ind1, ind2 in successive intervals
    else if (ind2 == ind1 + 1) {
        P.start <- pmatrix.msm(x=x, t = times[ind1] - t1 , covariates=covariates[[ind1]], qmatrix=qlist[[ind1]], ...)
        P.end <- pmatrix.msm(x=x, t = t2 - times[ind2-1], covariates=covariates[[ind2]], qmatrix=qlist[[ind2]], ...)
        P <- P.start %*% P.end
    }
    ## ind1, ind2 separated by one or more whole intervals
    else {
        P.start <- pmatrix.msm(x=x, t = times[ind1] - t1, covariates=covariates[[ind1]], qmatrix=qlist[[ind1]], ...)
        P.end <- pmatrix.msm(x=x, t = t2 - times[ind2-1], covariates=covariates[[ind2]], qmatrix=qlist[[ind2]], ...)
        P.middle <- diag(x$qmodel$nstates)
        for (i in (ind1+1):(ind2-1)) {
            P.middle <- P.middle %*% pmatrix.msm(x=x, t = times[i] - times[i-1], covariates=covariates[[i]], qmatrix=qlist[[i]], ...)
        }
        P <- P.start %*% P.middle %*% P.end
    }

    ci <- if (!is.null(x)) match.arg(ci) else "none"
    P.ci <- switch(ci,
                   bootstrap = pmatrix.piecewise.ci.msm(x=x, t1=t1, t2=t2, times=times, covariates=covariates, cl=cl, B=B, cores=cores),
                   normal = pmatrix.piecewise.normci.msm(x=x, t1=t1, t2=t2, times=times, covariates=covariates, cl=cl, B=B),
                   none = NULL)
    res <- if (ci=="none") P else list(estimates = P, L=P.ci[,,1], U=P.ci[,,2])
    res
}

### Extract the mean sojourn times for given covariate values

sojourn.msm <- function(x, covariates = "mean", ci=c("delta","normal","bootstrap","none"), cl=0.95, B=1000)
{
    ci <- match.arg(ci)
    qmatrix <- qmatrix.msm(x, covariates, sojourn=TRUE, ci=ci, cl=cl, B=B)
    sojstates <- (1 : x$qmodel$nstates) [transient.msm(x)]
    if (x$foundse && (ci != "none")){
        soj <- qmatrix$sojourn[sojstates]
        names (soj) <- rownames(x$qmodel$qmatrix)[sojstates]
        sojse <- qmatrix$sojournSE[sojstates]
        sojl <- qmatrix$sojournL[sojstates]
        soju <- qmatrix$sojournU[sojstates]
        names(sojse) <- names(sojl) <- names(soju) <- names(soj)
        res <- data.frame(estimates=soj, SE=sojse, L=sojl, U=soju)
    }
    else if (ci != "none") {
        res <- list(estimates=qmatrix$sojourn[sojstates])
    }
    else res <- list(estimates=qmatrix[sojstates])
    res
}


### Extract the probabilities of occupying each state next

pnext.msm <- function(x, covariates="mean", ci=c("normal","bootstrap","delta","none"), cl=0.95, B=1000, cores=NULL)
{
    ci <- match.arg(ci)
    Q <- qmatrix.msm(x, covariates, ci="none")
    pnext <- - Q / diag(Q)
    pnext[x$qmodel$imatrix==0] <- 0
    p.ci <- array(0, dim=c(dim(pnext), 2))
    if (x$foundse && (ci != "none")){
        if (ci == "delta") {
            for (i in 1:x$qmodel$nstates) {
                for (j in 1:x$qmodel$nstates) {
                    if (pnext[i,j] > 0) {
                        se <- qratio.se.msm(x, c(i,j), c(i,i), covariates, cl)$se
                        lse <- qratio.se.msm(x, c(i,j), c(i,i), covariates, cl)$lse
                        p.ci[i,j,1] <- exp ( log(pnext[i,j]) - qnorm(1 - 0.5*(1 - cl)) * lse )
                        p.ci[i,j,2] <- exp ( log(pnext[i,j]) + qnorm(1 - 0.5*(1 - cl)) * lse )
                    }
                }
            }
        }
        else if (ci=="normal")
            p.ci <- pnext.normci.msm(x, covariates, cl, B)
        else if (ci=="bootstrap")
            p.ci <- pnext.ci.msm(x, covariates, cl, B, cores)
        res <- list(estimates=pnext, L=p.ci[,,1], U=p.ci[,,2])
    }
    else res <- list(estimates=pnext)
    class(res) <- "msm.est"
    res
}

### Extract the coefficients

coef.msm <- function(object, ...)
{
    if (!inherits(object, "msm")) stop("expected object to be a msm model")
    if (object$emodel$misc)
        object[c("Qmatrices", "Ematrices")]
    else object$Qmatrices
}

### Extract the log-likelihood

logLik.msm <- function(object, by.subject=FALSE, ...)
{
    if (!inherits(object, "msm")) stop("expected object to be a msm model")
    if (by.subject){
        p <- object$paramdata
        val <- -0.5*Ccall.msm(p$opt$par, do.what="lik.subj", expand.data(object), object$qmodel, object$qcmodel, object$cmodel, object$hmodel, p)
        names(val) <- unique(model.extract(object$data$mf, "subject"))
    } else {
        val <- - 0.5*object$minus2loglik
        attr(val, "df") <- object$paramdata$nopt
        class(val) <- "logLik"
    }
    val
}

### Likelihood ratio test between two or more models

lrtest.msm <- function(...){
    mods <- list(...)
    if (length(mods) < 2) stop("Expected 2 or more models as arguments")
    lx <- logLik(mods[[1]])
    res <- matrix(nrow=length(mods)-1, ncol=3)
    colnames(res) <- c("-2 log LR","df","p")
    rownames(res) <- sapply(as.list(match.call())[-(1:2)], deparse)
    for (i in 2:length(mods)) {
        if (!inherits(mods[[i]], "msm"))
            stop("Expected argument",i,"to be a msm object")
        ly <- logLik(mods[[i]])
        lr <- as.numeric(-2 * (lx - ly))
        df <- attr(ly,"df") - attr(lx,"df")
        res[i-1,] <- c(lr, df, 1 - pchisq(lr, df))
    }
    res
}

## Estimate total length of stay in a given state.

totlos.msm <- function(x, start=1, end=NULL, fromt=0, tot=Inf, covariates="mean",
                       piecewise.times=NULL,
                       piecewise.covariates=NULL,
                       num.integ=FALSE, discount=0,
                       env=FALSE,
                       ci=c("none","normal","bootstrap"), # calculate a confidence interval
                       cl = 0.95, # width of symmetric confidence interval
                       B = 1000, # number of bootstrap replicates
                       cores=NULL,
                       ...)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    nst <- x$qmodel$nstates
    if (!is.numeric(start) ||
        ((length(start)==1) && (! start %in% 1 : nst))) stop("start should be a state in 1, ..., ", nst, " or a vector of length ",nst)
    else if (length(start) == 1) {p0 <- rep(0, nst); p0[start] <- 1; start <- p0}
    else if (length(start) > 1) {
        if (length(start) != nst)
            stop("start should be a state in 1, ..., ", nst, " or a vector of length ",nst)
    }
    if (is.null(end)) end <- 1 : nst
    if (! all(end %in% 1 : nst)) stop("end should be a set of states in 1, ..., ", nst)
    if (!is.numeric(fromt) || !is.numeric(tot) || length(fromt) != 1 || length(tot) != 1 || fromt < 0 || tot < 0)
        stop("fromt and tot must be single non-negative numbers")
    if (fromt > tot) stop("tot must be greater than fromt")
    if (length(absorbing.msm(x)) == 0)
        if (tot==Inf) stop("Must specify a finite end time for a model with no absorbing state")

    ncuts <- length(piecewise.times)
    npieces <- length(piecewise.covariates)
    if (!is.null(piecewise.times) && (!is.numeric(piecewise.times) || is.unsorted(piecewise.times)))
        stop("piecewise.times should be a vector of numbers in increasing order")
    if (!is.null(piecewise.covariates) && (npieces != ncuts + 1))
        stop("Number of piecewise.covariate lists must be one greater than the number of cut points")
    if (is.null(piecewise.covariates)) {
        ## define homogeneous model as piecewise with one piece
        npieces <- 1
        covs <- list(covariates)
        ptimes <- c(fromt, tot)
    } else {
        ## ignore all cut points outside [fromt,tot]
        keep <- which((piecewise.times > fromt) & (piecewise.times < tot))
        ## cov value between fromt and min(first cut, tot)
        cov1 <- piecewise.covariates[findInterval(fromt, piecewise.times) + 1]
        covs <- c(cov1, piecewise.covariates[keep+1])
        npieces <- length(covs)
        ptimes <- c(fromt, piecewise.times[keep], tot)
    }

    tmat <- envmat <- matrix(nrow=npieces, ncol=nst)
    if (tot==Inf) {
        tmat[,absorbing.msm(x)] <- Inf # set by hand or else integrate() will fail
        envmat[,absorbing.msm(x)] <- 1
        rem <- setdiff(seq_len(nst), absorbing.msm(x))
    }
    else rem <- seq_len(nst)
    for (i in 1:npieces) {
        from.t <- ptimes[i]
        to.t <- ptimes[i+1]
        Q <- qmatrix.msm(x, covariates=covs[[i]], ci="none")
        if (num.integ || to.t==Inf){
            for (j in rem){
                f <- function(time) {
                    y <- numeric(length(time))
                    for (k in seq_along(y))
                        y[k] <- (start %*% pmatrix.msm(x, time[k], t1=0, covariates=covs[[i]], ci="none")) [j]
                    y
                }
                tmat[i,j] <- integrate(f, from.t, to.t, ...)$value
            }
        } else {
            QQ <- rbind(c(0, start),
                        cbind(rep(0,nst), Q - discount*diag(nst)))
            tmat[i,] <- as.vector(c(1, rep(0, nst)) %*%
                                  (MatrixExp(to.t*QQ) - MatrixExp(from.t*QQ)) %*%
                                  rbind(rep(0, nst), diag(nst)))
        }
        Q0 <- Q; diag(Q0) <- 0
        envmat[i,rem] <- tmat[i,rem] %*% Q0[rem,rem]
    }
    res <- if (env) colSums(envmat) else colSums(tmat)
    names(res) <- rownames(x$qmodel$qmatrix)
    ci <- match.arg(ci)
    t.ci <- switch(ci,
                   bootstrap = totlos.ci.msm(x=x, start=start, end=end, fromt=fromt, tot=tot, covariates=covariates, piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates, discount=discount, env=env, cl=cl, B=B, cores=cores, ...),
                   normal = totlos.normci.msm(x=x, start=start, end=end, fromt=fromt, tot=tot, covariates=covariates, piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates, discount=discount, env=env, cl=cl, B=B, ...),
                   none = NULL)
    if (ci=="none") res[end] else rbind(res, t.ci)[,end]
}

## Expected number of visits

envisits.msm <- function(x=NULL, start=1, end=NULL, fromt=0, tot=Inf, covariates="mean",
                         piecewise.times=NULL,  piecewise.covariates=NULL,
                         num.integ=FALSE, discount=0,
                         ci=c("none","normal","bootstrap"), # calculate a confidence interval
                         cl = 0.95, # width of symmetric confidence interval
                         B = 1000, # number of bootstrap replicates
                         cores=NULL,
                         ...)
{
    totlos.msm(x=x, start=start, end=end, fromt=fromt, tot=tot, covariates=covariates, piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates, num.integ=num.integ, discount=discount, env=TRUE, ci=ci, cl=cl, B=B, cores=cores, ...)
}

## Return indices of transient states (can either call for a fitted model or a qmatrix)

transient.msm <- function(x=NULL, qmatrix=NULL)
{
    if (!is.null(x)) {
        if (!inherits(x, "msm")) stop("expected x to be a msm model")
        qmatrix <- x$Qmatrices[[1]]
        nst <- x$qmodel$nstates
    }
    else if (!is.null(qmatrix)) {
        nst <- nrow(qmatrix)
    }
    else stop("Neither a fitted msm model nor a qmatrix have been supplied")
    which(diag(msm.fixdiag.qmatrix(qmatrix)) != 0)
}

## Return indices of absorbing states (can either call for a fitted model or a qmatrix)

absorbing.msm <- function(x=NULL, qmatrix=NULL)
{
    if (!is.null(x)) {
        if (!inherits(x, "msm")) stop("expected x to be a msm model")
        qmatrix <- x$Qmatrices[[1]]
        nst <- x$qmodel$nstates
    }
    else if (!is.null(qmatrix)) {
        nst <- nrow(qmatrix)
    }
    else stop("Neither a fitted msm model nor a qmatrix have been supplied")
    which(diag(msm.fixdiag.qmatrix(qmatrix)) == 0)
}

## Return two-column matrix containing pairs of states with allowed
## transitions in an interval.  Handles transitions between observed
## states in misclassification models

intervaltrans.msm <- function(x=NULL, qmatrix=NULL, ematrix=NULL, exclude.absabs=FALSE, censor=FALSE) {
    if (!is.null(x)) {
        if (!inherits(x, "msm")) stop("expected x to be a msm model")
        qmatrix <- qmatrix.msm(x, ci="none")
        if (is.null(ematrix) & x$emodel$misc)
            ematrix <- ematrix.msm(x, ci="none") > 0
        abs <- absorbing.msm(x)
    }
    else if (!is.null(qmatrix)) {
        abs <- absorbing.msm(qmatrix=qmatrix)
    }
    else if (is.null(qmatrix))
        stop("Neither a fitted msm model nor a qmatrix have been supplied")
    P <- MatrixExp(qmatrix)
    if (!is.null(ematrix))
        P <- t(ematrix) %*% P %*% ematrix # > 0 iff P(obs state=s | prev obs state = r) > 0
    ## P(obs state = s | obs prev = r)  =  Sum_ij  P(obsst = s | truest = j) P(truest = j | trueprev = i) P(trueprev = i | obsprev = r)
    ##  Sum_ij   Ejs Pij Eir    =  Eir Pij Ejs
    gt0 <- abs(P) > .Machine$double.eps ^ 0.5
    at <- cbind(row(P)[gt0], col(P)[gt0])
    if (exclude.absabs)
        at <- at[!(at[,1] %in% abs & at[,2] %in% abs),]
    if (censor && x$cmodel$ncens > 0) { # consider censoring as separate state
        atcens.all <- numeric()
        for (i in 1:x$cmodel$ncens) {
            truestates <- x$cmodel$states[x$cmodel$index[i] : (x$cmodel$index[i+1] - 1)]
            atcens <- at[at[,2] %in% truestates,]
            atcens[,2] <- 99
            atcens <- unique(atcens)
            atcens.all <- rbind(atcens.all, atcens)
        }
        at <- rbind(at, atcens.all)
    }
    at[order(at[,1],at[,2]),]
}

## Enumerate the distinct subject covariate histories in the data
## Works for time homogeneous and inhomogeneous models
## Used to calculate expected prevalences for a population with those covariates

get.covhist <- function(x, subset=NULL) {
    ## Keep only times where the covariate changes, or first or last obs
    mf <- x$data$mf
    if (x$qcmodel$ncovs > 0) {
        if (!is.null(subset)) {
            subs <- mf$"(subject)" %in% subset
            mf <- mf[subs,,drop=FALSE]
        }
        subj <- match(mf$"(subject)", unique(mf$"(subject)"))
        n <- length(subj)
        apaste <- do.call("paste", mf[,attr(mf,"covnames"),drop=FALSE])
        first <- !duplicated(subj); last <- rev(!duplicated(rev(subj)))
        keep <- (c(0, apaste[1:(n-1)]) != apaste) | first | last
        ## Keep and tabulate unique covariate series
        covseries <- split(apaste[keep], subj[keep]) # as a list of char vectors
        covseries <- sapply(covseries, paste, collapse=" , ") # as one char vector, one series per pt.
        ## also need p matrices for different times as well as different covs.
        ## but only interested in cov change times if there's more than one
        ## transition (at least one times change point)
        change.times <- mf$"(time)"; change.times[first] <- change.times[last] <- 0
        change.times <- split(change.times[keep & (!(first|last))], subj[keep & (!(first|last))])
        change.times <- sapply(change.times, paste, collapse= " , ")
        covseries.t <- paste(covseries, change.times, sep="; ")
        ids <- unique(subj)[!duplicated(covseries.t)] # subj ids, one with each distinct series
        ncombs <- table(covseries.t)[unique(covseries.t)]# how many per series
        covmat <- cbind(subject=subj, time=mf$"(time)", mf[,attr(mf,"covnames"),drop=FALSE])
        covmat <- covmat[(subj %in% ids) & keep,]
        list(example=covmat, # rows of the original data sufficient to define the distinct histories
             hist=covseries.t) # one per subject listing their covariate history as a string
    }
    else NULL
}

### Estimate observed state occupancies in the data at a series of times
### Assume previous observed state is retained until next observation time
### Assumes times are sorted within patient (they are in data in msm objects)

observed.msm <- function(x, times=NULL, interp=c("start","midpoint"), censtime=Inf, subset=NULL)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    ## For general HMMs use the Viterbi estimate of the observed state.
    if (!is.null(x$pci)) {
        state <- x$data$mf$"(state)"[!x$data$mf$"(pci.imp)"]
        time <- x$data$mf$"(time)"[!x$data$mf$"(pci.imp)"]
        subject <- x$data$mf$"(subject)"
        subject <- subject[!x$data$mf$"(pci.imp)"]
    } else {
        if ((x$hmodel$hidden && !x$emodel$misc) ||
            (!x$emodel$misc && x$cmodel$ncens>0) )
            state <- viterbi.msm(x)$fitted
        else if (x$emodel$misc && x$cmodel$ncens>0) {
            vit <- viterbi.msm(x)$fitted
            state <- x$data$mf$"(state)"
            state[state %in% x$cmodel$censor] <- vit[state %in% x$cmodel$censor]
            ## TODO for misc models with censoring, impute only censored obs states from viterbi
        }  else
            state <- x$data$mf$"(state)"

        time <- x$data$mf$"(time)"; subject <- x$data$mf$"(subject)"
    }
    if (is.null(subset)) subset <- unique(subject)
    time <- time[subject %in% subset]
    state <- state[subject %in% subset]
    subject <- subject[subject %in% subset]
    if (is.null(times))
        times <- seq(min(time), max(time), (max(time) - min(time))/10)
    states.expand <- matrix(nrow=length(unique(subject)), ncol=length(times))
    pts <- unique(subject)
    absorb <- absorbing.msm(x)
    interp <- match.arg(interp)
    if (!is.numeric(censtime)) stop("censtime should be numeric")
    if (length(censtime)==1) censtime <- rep(censtime, length(pts))
    else if (length(censtime)!=length(pts)) stop("censtime of length ", length(censtime), ", should be 1 or ", length(pts))
    for (i in seq_along(pts)){
        state.i <- state[(subject==pts[i])]
        time.i <- time[(subject==pts[i])]
        j <- 1
        while(j <= length(times)) {
            if (times[j] < time.i[1]) {
                mtime <- max(which(times-time.i[1] < 0))
                states.expand[i, j:mtime] <- NA
                j <- mtime + 1
                next;
            } else if (times[j] > time.i[length(time.i)]) {
                if (state.i[length(time.i)] %in% absorb && (times[j] <= censtime[i])) {
                    states.expand[i, j:(length(times))] <-  state.i[length(time.i)]
                } else states.expand[i, j:(length(times))] <-  NA
                break;
            } else {
                prevtime.ind <- max(which(time.i <= times[j]))
                prevtime <- time.i[prevtime.ind]
                if (interp=="midpoint") {
                    nexttime.ind <- min(which(time.i >= times[j]))
                    nexttime <- time.i[nexttime.ind]
                    midpoint <- (prevtime + nexttime) / 2
                    states.expand[i,j] <- state.i[if (times[j] <= midpoint) prevtime.ind else nexttime.ind]
                } else
                    states.expand[i,j] <- state.i[prevtime.ind]
            }
            j <- j+1
        }
    }
    obstab <- t(apply(states.expand, 2, function(y) table(factor(y, levels=seq(length=x$qmodel$nstates)))))
    obsperc <- 100*obstab / rep(rowSums(obstab), ncol(obstab))
    dimnames(obstab) <- dimnames(obsperc) <- list(times, paste("State", 1:x$qmodel$nstates))
    obstab <- cbind(obstab, Total=rowSums(obstab))

    covhist <- get.covhist(x, subset)
    covcat <- ## distinct covariate history group each subject falls into (ordinal)
        if (is.null(covhist)) rep(1, length(unique(subject)))
        else match(covhist$hist, unique(covhist$hist))
    risk <- matrix(nrow=length(times), ncol=length(unique(covcat)), dimnames = list(times, unique(covhist$hist)))
    for (i in seq_along(unique(covcat))) {
        obst <- t(apply(states.expand[covcat==unique(covcat)[i],,drop=FALSE], 2,
                        function(y) table(factor(y, levels=seq(length=x$qmodel$nstates)))))
        risk[,i] <- rowSums(obst)
    }

    list(obstab=obstab, obsperc=obsperc, risk=risk)
}


expected.msm <- function(x,
                         times=NULL,
                         timezero=NULL,
                         initstates=NULL,
                         covariates="population",
                         misccovariates="mean",
                         piecewise.times=NULL,
                         piecewise.covariates=NULL,
                         risk=NULL,
                         subset=NULL,
                         ci=c("none","normal","bootstrap"),
                         cl = 0.95,
                         B = 1000,
                         cores = NULL
                         )
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    time <- model.extract(x$data$mf, "time")
    if (is.null(times))
        times <- seq(min(time), max(time), (max(time) - min(time))/10)
    if (is.null(timezero))  timezero <- min(time)
    if (is.null(risk)) risk <- observed.msm(x, times=times, subset=subset)$risk
    exptab <- matrix(0, nrow=length(times), ncol=x$qmodel$nstates)
    start <- min(which(times - timezero >= 0))
    if (x$emodel$misc)
        initprobs <- x$emodel$initprobs
    else {
        if (is.null(initstates))
            initstates <- observed.msm(x, times=timezero)$obstab[1:x$qmodel$nstates]
        initprobs <- initstates / sum(initstates)
    }
    if (length(times) >= start) {
        for (j in start:length(times)) {
            if (x$qcmodel$ncovs>0 && covariates=="population") {
                covmat <- get.covhist(x, subset=subset)$example
                for (i in 1:length(unique(covmat$subject))) {
                    ## sum expected prevalences for each covariate history observed in the data
                    subji <- unique(covmat$subject)[i]
                    ni <- sum(covmat$subject==subji)
                    ctimes <- covmat$time[covmat$subject==subji][-c(1,ni)]
                    covs <- covmat[covmat$subject==subji, attr(x$data$mf,"covnames"),drop=FALSE][-ni,,drop=FALSE]
                    ccovs <- list()
                    for (k in 1:nrow(covs)) ccovs[[k]] <- as.list(covs[k,,drop=FALSE])
                    pmat <-  pmatrix.piecewise.msm(x, t1=timezero, t2=times[j], times=ctimes, covariates=ccovs)
                    expji <- risk[j,i] * initprobs %*% pmat
                    if (x$emodel$misc) { # return expected prev of obs (not true) states
                        if (x$ecmodel$ncovs==0) emat <- ematrix.msm(x, ci="none")
                        else {
                            ecovs <- if(length(ctimes)==0) ccovs else ccovs[[length(ccovs)]]
                            emat <- ematrix.msm(x, covariates=ecovs, ci="none")
                        }
                        expji <- expji %*% emat
                    }
                    exptab[j,] <- exptab[j,] + expji
                }
            }
            else {
                pmat <-
                    if (is.null(piecewise.times))
                        pmatrix.msm(x, t=times[j] - timezero, t1=timezero, covariates=covariates)
                    else
                        pmatrix.piecewise.msm(x, timezero, times[j], piecewise.times, piecewise.covariates)
                expj <- rowSums(risk)[j] * initprobs %*% pmat
                if (x$emodel$misc) # return expected prev of obs (not true) states
                    expj <- expj %*% ematrix.msm(x, covariates=misccovariates, ci="none")
                exptab[j,] <- expj
            }
        }
    }
    exptab <- cbind(exptab, apply(exptab, 1, sum))
    dimnames(exptab) <- list(times, c(rownames(x$qmodel$qmatrix),"Total"))
    expperc <- 100*exptab[,1:x$qmodel$nstates] / exptab[, x$qmodel$nstates+1]

    ci <- match.arg(ci)
    e.ci <- switch(ci,
                   bootstrap = expected.ci.msm(x, times, timezero, initstates, covariates, misccovariates,
                   piecewise.times, piecewise.covariates, risk, cl, B, cores),
                   normal = expected.normci.msm(x, times, timezero, initstates, covariates, misccovariates,
                   piecewise.times, piecewise.covariates, risk, cl, B),
                   none = NULL)
    res <-
        if (ci=="none") list(exptab=exptab, expperc=expperc)
        else list(exptab=list(estimates=exptab, ci=e.ci[[1]]),
                  expperc=list(estimates=expperc, ci=e.ci[[2]]))
    names(res) <- c("Expected","Expected percentages")
    res
}

### Table of observed and expected prevalences (works for misclassification and non-misclassification models)

prevalence.msm <- function(x,
                           times=NULL,
                           timezero=NULL,
                           initstates=NULL,
                           covariates="population",
                           misccovariates="mean",
                           piecewise.times=NULL,
                           piecewise.covariates=NULL,
                           ci=c("none","normal","bootstrap"),
                           cl = 0.95,
                           B = 1000,
                           cores = NULL,
                           interp=c("start","midpoint"),
                           censtime=Inf,
                           subset=NULL,
                           plot = FALSE, ...
                           )
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    ## Estimate observed state occupancies in the data at a series of times
    time <- model.extract(x$data$mf, "time")
    if (is.null(times))
        times <- seq(min(time), max(time), (max(time) - min(time))/10)
    obs <- observed.msm(x, times, interp, censtime, subset)
    ## Work out expected state occupancies by forecasting from transition probabilities
    expec <- expected.msm(x, times, timezero, initstates, covariates, misccovariates,
                          piecewise.times, piecewise.covariates, obs$risk, subset, ci, cl, B, cores)
    res <- list(observed=obs$obstab, expected=expec[[1]], obsperc=obs$obsperc, expperc=expec[[2]])
    names(res) <- c("Observed", "Expected", "Observed percentages", "Expected percentages")
    if (plot) plot.prevalence.msm(x, mintime=min(times), maxtime=max(times), timezero=timezero, initstates=initstates,
                                  interp=interp, censtime=censtime, covariates=covariates, misccovariates=misccovariates,
                                  piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates, ...)
    res
}

plot.prevalence.msm <- function(x, mintime=NULL, maxtime=NULL, timezero=NULL, initstates=NULL,
                                interp=c("start","midpoint"), censtime=Inf, subset=NULL,
                                covariates="population", misccovariates="mean",
                                piecewise.times=NULL, piecewise.covariates=NULL, xlab="Times",ylab="Prevalence (%)",
                                lwd.obs=1, lwd.exp=1, lty.obs=1, lty.exp=2,
                                col.obs="blue", col.exp="red", legend.pos=NULL,...){
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    time <- model.extract(x$data$mf, "time")
    if (is.null(mintime)) mintime <- min(time)
    if (is.null(maxtime)) maxtime <- max(time)
    t <- seq(mintime, maxtime, length=100)
    obs <- observed.msm(x, t, interp, censtime, subset)
    expec <- expected.msm(x, t, timezero=timezero, initstates=initstates, covariates=covariates, misccovariates=misccovariates,
                          piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates, risk=obs$risk, subset=subset, ci="none")[[2]]
    states <- seq(length=x$qmodel$nstates)
    S <- length(states)
    ncols <- ceiling(sqrt(S))
    nrows <- if (floor(sqrt(S))^2 < S && S <= floor(sqrt(S))*ceiling(sqrt(S))) floor(sqrt(S)) else ceiling(sqrt(S))
    par(mfrow=c(nrows, ncols))
    for (i in states) {
        plot(t, obs$obsperc[,i], type="l", ylim=c(0, 100), xlab=xlab, ylab=ylab, lwd=lwd.obs, lty=lty.obs, col=col.obs,
             main=rownames(x$qmodel$qmatrix)[i],...)
        lines(t, expec[,i], lwd=lwd.exp, lty=lty.exp, col=col.exp)
    }
    if (!is.numeric(legend.pos) || length(legend.pos) != 2)
        legend.pos <- c(0.4*maxtime, 40)
    legend(x=legend.pos[1], y=legend.pos[2], legend=c("Observed","Expected"), lty=c(lty.obs,lty.exp), lwd=c(lwd.obs,lwd.exp), col=c(col.obs,col.exp))
    invisible()
}

### Empirical versus fitted survival curve

plot.survfit.msm <- function(x, from=1, to=NULL, range=NULL, covariates="mean",
                             interp=c("start","midpoint"), ci=c("none","normal","bootstrap"), B=100,
                             legend.pos=NULL, xlab="Time", ylab="Survival probability",
                             lty=1, lwd=1, col="red", lty.ci=2, lwd.ci=1, col.ci="red",
                             mark.time=TRUE, col.surv="blue", lty.surv=2, lwd.surv=1,
                             survdata=FALSE, 
                             ...) {
    if (!inherits(x, "msm")) stop("expected \"x\" to be a msm model")
    if (is.null(to))
        to <- max(absorbing.msm(x))
    else {
        if (!is.numeric(to)) stop("\"to\" must be numeric")
        if (! (to %in% absorbing.msm(x) ) ) stop("\"to\" must be an absorbing state")
    }
    if (! (from %in% transient.msm(x) ) ) stop("\"from\" must be a non-absorbing state")
    if (is.null(range))
        rg <- range(model.extract(x$data$mf, "time"))
    else {
        if (!is.numeric(range) || length(range)!= 2) stop("\"range\" must be a numeric vector of two elements")
        rg <- range
    }
    interp <- match.arg(interp)
    ci <- match.arg(ci)
    timediff <- (rg[2] - rg[1]) / 50
    times <- seq(rg[1], rg[2], timediff)
    pr <- lower <- upper <- numeric()

    for (t in times) {
        P <- pmatrix.msm(x, t, t1=times[1], covariates=covariates, ci=ci, B=B)
        if (ci != "none") {
            pr <- c(pr, P$estimates[from, to])
            lower <- c(lower, P$L[from, to])
            upper <- c(upper, P$U[from, to])
        }
        else pr <- c(pr, P[from, to])
    }
    plot(times, 1 - pr, type="l", xlab=xlab, ylab=ylab, lwd=lwd,
         ylim=c(0,1), lty = lty, col=col,...)
    if (ci != "none") {
        lines(times, 1 - lower, lty=lty.ci, col=col.ci, lwd=lwd.ci)
        lines(times, 1 - upper, lty=lty.ci, col=col.ci, lwd=lwd.ci)
    }
    dat <- x$data$mf[,c("(subject)", "(time)", "(state)")]
    dat$"(subject)" <- match(dat$"(subject)", unique(dat$"(subject)")) # guard against factor
    dat$subjstate <- paste(dat$"(subject)", dat$"(state)")

    ## restrict data to subjects with any observations of "from" 
    anyfrom <- tapply(dat$"(state)", dat$"(subject)", function(x)any(x==from))[as.character(dat$"(subject)")]
    dat <- dat[anyfrom,]
    obspt <- sequence(table(dat$"(subject)")) # observation numbers, starting at 1 for each subject
    ## number of first observation of "from" for current subject 
    minfrom <- rep(which(dat$"(state)"==from  & !duplicated(dat$subjstate)) - which(!duplicated(dat$"(subject)")), table(dat$"(subject)")) + 1
    ## restrict data to observations after and including first observation of "from" for each person 
    dat <- dat[obspt>=minfrom,]
    
    first <- !duplicated(dat$"(subject)")
    last <- !duplicated(dat$"(subject)", fromLast=TRUE) 
    ## does each subject have an absorbing state
    anyabs <- tapply(dat$"(state)", dat$"(subject)", function(x)any(x==to))[as.character(dat$"(subject)")]
    subjstate <- paste(dat$"(subject)", dat$"(state)")
    ## index of first occurrence of absorbing state, or last obs if no absorbing state
    minabs <- dat$"(state)"==to  & !duplicated(subjstate)
    dtime <- dat$"(time)" - tapply(dat$"(time)", dat$"(subject)", min)[as.character(dat$"(subject)")]
    if (interp=="midpoint"){
        ## index of state just before first occurrence of abs state
        prevminabs <- c(minabs[-1], FALSE)
        dtime[minabs] <- 0.5*(dtime[minabs] + dtime[prevminabs])
    }
    minabs[!anyabs] <- last[!anyabs]
    survdat <- data.frame(survtime =  dtime[minabs], died = as.numeric(anyabs[minabs]))

    lines(survfit(Surv(survdat$survtime, survdat$died) ~ 1), mark.time=mark.time, col=col.surv, lty=lty.surv, lwd=lwd.surv,...)
    timediff <- (rg[2] - rg[1]) / 50
    if (!is.numeric(legend.pos) || length(legend.pos) != 2)
        legend.pos <- c(max(x$data$mf$"(time)") - 25*timediff, 1)
    if (ci=="none")
        legend(legend.pos[1], legend.pos[2], lty=c(lty, lty.surv), lwd=c(lwd, lwd.surv), col=c(col, col.surv),
               legend=c("Fitted","Empirical"))
    else legend(legend.pos[1], legend.pos[2], lty=c(lty, lty.ci, lty.surv), lwd=c(lwd,lwd.ci, lwd.surv), col=c(col ,col.ci, col.surv),
                legend=c("Fitted","Fitted (confidence interval)", "Empirical"))

    if (survdata) survdat else invisible()
}

### Obtain hazard ratios from estimated effects of covariates on log-transition rates

hazard.msm <- function(x, hazard.scale = 1, cl = 0.95)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    if (length(hazard.scale) == 1) hazard.scale <- rep(hazard.scale, x$qcmodel$ncovs)
    if (length(hazard.scale) != x$qcmodel$ncovs)
        stop ("hazard.scale of length ", length(hazard.scale), ", expected ", x$qcmodel$ncovs)
    keep <- (x$qmodel$imatrix != 0)
    nst <- x$qmodel$nstates
    keepvec <- as.vector(t(keep))
    fromlabs <- rep(rownames(keep), each=nst) [keepvec]
    tolabs <- rep(colnames(keep), nst) [keepvec]
    if (x$qcmodel$ncovs > 0) {
        haz.list <- list()
        if (x$foundse) {
            for (i in 1:x$qcmodel$ncovs) {
                cov <- x$qcmodel$covlabels[i]
                haz.rat <- t(exp(hazard.scale[i]*x$Qmatrices[[cov]]))[keepvec]
                LCL <- t(exp(hazard.scale[i]*(x$Qmatrices[[cov]] - qnorm(1 - 0.5*(1 - cl))*x$QmatricesSE[[cov]]) ))[keepvec]
                UCL <- t(exp(hazard.scale[i]*(x$Qmatrices[[cov]] + qnorm(1 - 0.5*(1 - cl))*x$QmatricesSE[[cov]]) ))[keepvec]
                haz.tab <- cbind(haz.rat, LCL, UCL)
                dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs),
                                          c("HR", "L", "U"))
                haz.list[[cov]] <- haz.tab
            }
        }
        else {
            for (i in 1:x$qcmodel$ncovs) {
                cov <- x$qcmodel$covlabels[i]
                haz.tab <- as.matrix(t(exp(hazard.scale[i]*x$Qmatrices[[cov]]))[keepvec])
                dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs), "HR")
                haz.list[[cov]] <- haz.tab
            }
        }
    }
    else haz.list <- "No covariates on transition intensities"
    haz.list
}


### Obtain odds ratios from estimated effects of covariates on logit-misclassification probabilities
### TODO - equivalent for general HMMs which presents cov effects on natural scale.

odds.msm <- function(x, odds.scale = 1, cl = 0.95)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    if (!x$emodel$misc) stop("Requires a misclassification model specified with ematrix")
    if (length(odds.scale) == 1) odds.scale <- rep(odds.scale, x$ecmodel$ncovs)
    if (length(odds.scale) != x$ecmodel$ncovs)
        stop ("odds.scale of length ", length(odds.scale), ", expected ", x$ecmodel$ncovs)
    keep <- (x$emodel$imatrix != 0)
    nst <- x$qmodel$nstates
    keepvec <- as.vector(t(keep))
    truelabs <- rep(rownames(keep), each=nst) [keepvec]
    obslabs <- rep(colnames(keep), nst) [keepvec]
    if (x$ecmodel$ncovs > 0) {
        odds.list <- list()
        if (x$foundse) {
            for (i in 1:x$ecmodel$ncovs) {
                cov <- x$ecmodel$covlabels[i]
                odds.rat <- t(exp(odds.scale[i]*x$Ematrices[[cov]]))[keepvec]
                LCL <- t(exp(odds.scale[i]*(x$Ematrices[[cov]] - qnorm(1 - 0.5*(1 - cl))*x$EmatricesSE[[cov]]) ))[keepvec]
                UCL <- t(exp(odds.scale[i]*(x$Ematrices[[cov]] + qnorm(1 - 0.5*(1 - cl))*x$EmatricesSE[[cov]]) ))[keepvec]
                odds.tab <- cbind(odds.rat, LCL, UCL)
                dimnames(odds.tab) <- list(paste("Obs", obslabs, "|", truelabs),
                                           c("OR", "L", "U"))
                odds.list[[cov]] <- odds.tab
            }
        }
        else {
            for (i in 1:x$ecmodel$ncovs) {
                cov <- x$ecmodel$covlabels[i]
                odds.tab <- as.matrix(t(exp(odds.scale[i]*x$Ematrices[[cov]]))[keepvec])
                dimnames(odds.tab) <- list(paste("Obs", obslabs, "|", truelabs), "OR")
                odds.list[[cov]] <- odds.tab
            }
        }
    }
    else odds.list <- "No covariates on misclassification probabilities"
    odds.list
}

viterbi.msm <- function(x, normboot=FALSE)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    if (x$cmodel$ncens > 0 && !x$hmodel$hidden) {
        ## If censoring but not HMM, then define an identity HMM with
        ## true state known at every time except censoring times
        hmod <- vector(x$qmodel$nstates, mode="list")
        for (i in 1:x$qmodel$nstates)
            hmod[[i]] <- hmmIdent(i)
        x$hmodel <- msm.form.hmodel(hmod, est.initprobs=FALSE)
        x$hmodel <- c(x$hmodel, list(ncovs=rep(rep(0,x$hmodel$nstates),x$hmodel$npars), ncoveffs=0, nicovs=rep(0,x$hmodel$nstates-1), nicoveffs=0))
        x$data$mf$"(obstrue)" <- ifelse(x$data$mf$"(state)" %in% x$cmodel$censor, 0, (x$data$mf$"(state)"))
        x$data$mm.hcov <- vector(mode="list", length=x$hmodel$nstates) # reqd by msm.add.hmmcovs
        for (i in seq_len(x$hmodel$nstates))
            x$data$mm.hcov[[i]] <- model.matrix(~1, x$data$mf)
        x$paramdata$allinits <- c(x$paramdata$allinits,x$hmodel$pars)
        x$paramdata$constr <- c(x$paramdata$constr,max(x$paramdata$constr)+seq_along(x$hmodel$pars))
    }
    if (x$hmodel$hidden) {        
        params <-
          if (normboot)
              rmvnorm(1, x$paramdata$opt$par, x$covmat[x$paramdata$optpars,x$paramdata$optpars])
          else x$paramdata$opt$par
        ret <- Ccall.msm(params, do.what="viterbi", expand.data(x),
                               x$qmodel, x$qcmodel, x$cmodel, x$hmodel, x$paramdata)
        fitted <- ret[[1]]; pstate <- ret[[2]]
        fitted <- fitted + 1
    }
    else { 
        fitted <- x$data$mf$"(state)"
        pstate <- NULL
    }
    if (!is.null(x$qmodel$phase.states)){
        fitted <- x$qmodel$phase.labs[fitted]
    }
    ret <- data.frame(subject = x$data$mf$"(subject)",
               time = x$data$mf$"(time)",
               observed = x$data$mf$"(state)",
               fitted = fitted)
    if (!is.null(pstate)) ret$pstate <- pstate
    ret
}

scoreresid.msm <- function(x, plot=FALSE){
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    if (!deriv.supported(x$data, x$hmodel, x$cmodel))
        stop("Score residuals not available, since analytic derivatives not implemented for this model")
    derivs <- Ccall.msm(x$paramdata$opt$par, do.what="deriv.subj", expand.data(x), x$qmodel, x$qcmodel, x$cmodel, x$hmodel, x$paramdata)
    cov <- x$paramdata$covmat[x$paramdata$optpars,x$paramdata$optpars]
    sres <- colSums(t(derivs) * cov %*% t(derivs))
    names(sres) <- unique(x$data$mf$"(subject)")
    if (plot) {
        plot(sres, type="n")
        text(seq_along(sres), sres, names(sres))
    }
    sres
}

# Function to calculate expected first passage times for continuous-time Markov chain with arbitrary Q matrix
# Returns vector with EFPT for each "from" state in the state space.
# Could also get CDF simply by making tostate absorbing and calculating pmatrix.
# TODO time-dependent covariates.  Unclear if the expectation has a solution for piecewise-constant rate.

efpt.msm <- function(x=NULL, qmatrix=NULL, tostate, start="all", covariates="mean",
                     ci=c("none","normal","bootstrap"), cl = 0.95, B = 1000, cores=NULL, ...)
{
    ci <- match.arg(ci)
    if (!is.null(x)) {
        if (!inherits(x, "msm")) stop("expected x to be a msm model")
        qmatrix <- qmatrix.msm(x, covariates=covariates, ci="none")
    }
    else if (!is.null(qmatrix)) {
        if (!is.matrix(qmatrix) || (nrow(qmatrix) != ncol(qmatrix)))
            stop("expected qmatrix to be a square matrix")
        if (ci != "none") {warning("No fitted model supplied: not calculating confidence intervals."); ci <- "none"}
    }
    if (is.character(tostate)) {
        if (!tostate %in% rownames(qmatrix)) stop(sprintf("state \"%s\" unknown", tostate))
        tostate <- match(tostate, rownames(qmatrix))
    }
    est <- rep(NA, nrow(qmatrix))
    ## EFPT is zero if we're already in tostate
    est[tostate] <- 0
    abstate <- absorbing.msm(qmatrix=qmatrix)
    ## EFPT is infinite for other absorbing states
    est[setdiff(abstate,tostate)] <- Inf
    fromstate <- setdiff(1:nrow(qmatrix), union(abstate,tostate))
    ## EFPT is infinite if any chance of absorbing elsewhere before
    ## hitting tostate.  To calculate this, form Q matrix with tostate
    ## made absorbing, and look at P matrix in unit time.
    Qred <- qmatrix; Qred[tostate,] <- 0
    Pmat <- MatrixExp(Qred, ...)
    Pmat[Pmat < 1e-16] <- 0
    p.abs <- rowSums(Pmat[fromstate,setdiff(abstate,tostate),drop=FALSE])
    est[fromstate][p.abs>0] <- Inf

    ## Any states left from which EFPT to tostate is nonzero and finite.
    ## Use standard linear equation solution
    ## see, e.g. equation (3) of Harrison and Knottenbelt (2001)
    if (any(is.na(est))){
        fromstate <- which(is.na(est))
        Q <- as.matrix(qmatrix[fromstate, fromstate])
        est[fromstate] <- solve(-Q, rep(1,nrow(Q)))
    }
    if (!is.character(start)) {
        if (!is.numeric(start) || (length(start)!=nrow(qmatrix)))
            stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmatrix))
        start <- start / sum(start)
        est <- est %*% start
    }
    else if (any(start!="all"))
        stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmatrix))
    e.ci <- switch(ci,
                   bootstrap = efpt.ci.msm(x=x, qmatrix=qmatrix, tostate=tostate, start=start, covariates=covariates, cl=cl, B=B, cores=cores),
                   normal = efpt.normci.msm(x=x, qmatrix=qmatrix, tostate=tostate, start=start, covariates=covariates, cl=cl, B=B),
                   none = NULL)
    if (ci=="none") est else rbind(est, e.ci)
}


## TODO time-inhomogeneous models.
## Do msm.fill.pci.covs to get covariate list
## Get list of Qs by calling qmatrix.msm for each member of this list
## Zero out the appropriate rows
## Supply this to pmatrix.piecewise.msm, which will call pmatrix.msm

ppass.msm <- function(x=NULL, qmatrix=NULL, tot, start="all", covariates="mean",
                      piecewise.times=NULL, piecewise.covariates=NULL,
                      ci=c("none","normal","bootstrap"), cl = 0.95, B = 1000, cores=NULL, ...)
{
    ci <- match.arg(ci)
    if (!is.null(x)) {
        if (!inherits(x, "msm")) stop("expected x to be a msm model")
        qmatrix <- qmatrix.msm(x, covariates=covariates, ci="none")
    }
    else if (!is.null(qmatrix)) {
        if (!is.matrix(qmatrix) || (nrow(qmatrix) != ncol(qmatrix)))
            stop("expected qmatrix to be a square matrix")
        if (ci != "none") {warning("No fitted model supplied: not calculating confidence intervals."); ci <- "none"}
    }
    res <- array(dim=dim(qmatrix))
    if (!is.null(dimnames(qmatrix))) {
        dimnames(res) <- dimnames(qmatrix)
        names(dimnames(res)) <- c("from","to")
    }
    states <- 1:nrow(qmatrix)
    for (i in states) {
        Qred <- qmatrix; Qred[states[i],] <- 0
        res[,i] <- MatrixExp(Qred*tot, ...)[,i]
    }
    if (!is.character(start)) {
        if (!is.numeric(start) || (length(start)!=nrow(qmatrix)))
            stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmatrix))
        start <- start / sum(start)
        res <- matrix(start %*% res, nrow=1, dimnames=list(from="start",to=colnames(res)))
    }
    else if (any(start!="all"))
        stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmatrix))
    p.ci <- switch(ci,
                   bootstrap = ppass.ci.msm(x=x, qmatrix=qmatrix, tot=tot, start=start, covariates=covariates, cl=cl, B=B, cores=cores),
                   normal = ppass.normci.msm(x=x, qmatrix=qmatrix, tot=tot, start=start, covariates=covariates, cl=cl, B=B),
                   none = NULL)
    if (ci != "none") {
        res <- list(estimates=res, L=p.ci$L, U=p.ci$U)
        class(res) <- "msm.est"
    }
    res
}
