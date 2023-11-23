### METHODS FOR MSM OBJECTS



#' Transition intensity matrix
#' 
#' Extract the estimated transition intensity matrix, and the corresponding
#' standard errors, from a fitted multi-state model at a given set of covariate
#' values.
#' 
#' Transition intensities and covariate effects are estimated on the log scale
#' by \code{\link{msm}}. A covariance matrix is estimated from the Hessian of
#' the maximised log-likelihood.
#' 
#' A more practically meaningful parameterisation of a continuous-time Markov
#' model with transition intensities \eqn{q_{rs}} is in terms of the mean
#' sojourn times \eqn{-1 / q_{rr}} in each state \eqn{r} and the probabilities
#' that the next move of the process when in state \eqn{r} is to state \eqn{s},
#' \eqn{-q_{rs} / q_{rr}}.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param covariates The covariate values at which to estimate the intensity
#' matrix.  This can either be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula. Or more clearly, a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' 
#' If some covariates are specified but not others, the missing ones default to
#' zero.
#' 
#' With \code{covariates="mean"}, for factor / categorical variables, the mean
#' of the 0/1 dummy variable for each factor level is used, representing an
#' average over all values in the data, rather than a specific factor level.
#' @param sojourn Set to TRUE if the estimated sojourn times and their standard
#' errors should also be returned.
#' @param ci If \code{"delta"} (the default) then confidence intervals are
#' calculated by the delta method, or by simple transformation of the Hessian
#' in the very simplest cases.  Normality on the log scale is assumed.
#' 
#' If \code{"normal"}, then calculate a confidence interval by simulating
#' \code{B} random vectors from the asymptotic multivariate normal distribution
#' implied by the maximum likelihood estimates (and covariance matrix) of the
#' log transition intensities and covariate effects, then transforming.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' @param cl Width of the symmetric confidence interval to present.  Defaults
#' to 0.95.
#' @param B Number of bootstrap replicates, or number of normal simulations
#' from the distribution of the MLEs.
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @return A list with components:
#' 
#' \item{estimate}{Estimated transition intensity matrix.}
#' 
#' \item{SE}{Corresponding approximate standard errors.}
#' 
#' \item{L}{Lower confidence limits}
#' 
#' \item{U}{Upper confidence limits}
#' 
#' Or if \code{ci="none"}, then \code{qmatrix.msm} just returns the estimated
#' transition intensity matrix.
#' 
#' If \code{sojourn} is \code{TRUE}, extra components called \code{sojourn},
#' \code{sojournSE}, \code{sojournL} and \code{sojournU} are included,
#' containing the estimates, standard errors and confidence limits,
#' respectively, of the mean sojourn times in each transient state.
#' 
#' The default print method for objects returned by \code{\link{qmatrix.msm}}
#' presents estimates and confidence limits. To present estimates and standard
#' errors, do something like
#' 
#' \code{qmatrix.msm(x)[c("estimates","SE")]}
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{pmatrix.msm}}, \code{\link{sojourn.msm}},
#' \code{\link{deltamethod}}, \code{\link{ematrix.msm}}
#' @keywords models
#' @export qmatrix.msm
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
            inds <- seq(length.out=x$qmodel$npars + x$qcmodel$npars)
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



#' Misclassification probability matrix
#' 
#' Extract the estimated misclassification probability matrix, and
#' corresponding confidence intervals, from a fitted multi-state model at a
#' given set of covariate values.
#' 
#' Misclassification probabilities and covariate effects are estimated on the
#' multinomial-logit scale by \code{\link{msm}}. A covariance matrix is
#' estimated from the Hessian of the maximised log-likelihood.  From these, the
#' delta method can be used to obtain standard errors of the probabilities on
#' the natural scale at arbitrary covariate values.  Confidence intervals are
#' estimated by assuming normality on the multinomial-logit scale.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}
#' @param covariates
#' 
#' The covariate values for which to estimate the misclassification probability
#' matrix.  This can either be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula, or a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' @param ci If \code{"delta"} (the default) then confidence intervals are
#' calculated by the delta method, or by simple transformation of the Hessian
#' in the very simplest cases.
#' 
#' If \code{"normal"}, then calculate a confidence interval by simulating
#' \code{B} random vectors from the asymptotic multivariate normal distribution
#' implied by the maximum likelihood estimates (and covariance matrix) of the
#' multinomial-logit-transformed misclassification probabilities and covariate
#' effects, then transforming back.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' @param cl Width of the symmetric confidence interval to present.  Defaults
#' to 0.95.
#' @param B Number of bootstrap replicates, or number of normal simulations
#' from the distribution of the MLEs
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @return A list with components:
#' 
#' \item{estimate}{Estimated misclassification probability matrix. The rows
#' correspond to true states, and columns observed states.}
#' \item{SE}{Corresponding approximate standard errors.} \item{L}{Lower
#' confidence limits.} \item{U}{Upper confidence limits.}
#' 
#' Or if \code{ci="none"}, then \code{ematrix.msm} just returns the estimated
#' misclassification probability matrix.
#' 
#' The default print method for objects returned by \code{\link{ematrix.msm}}
#' presents estimates and confidence limits. To present estimates and standard
#' errors, do something like
#' 
#' \code{ematrix.msm(x)[c("estimates","SE")]}
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{qmatrix.msm}}
#' @keywords models
#' @export
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
## fixme when called from bootstrap, hasn't worked. 
## is it because covfactor is all false?  is covdata.mf wrong?  yes both numeric 
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
            cov.inds <- sum(hmodel$npars) + (cur.i-1)*nc + seq(length.out=(nc*nir)) # indices into HMM parameter vector of corresp cov effects
            parinds <- numeric(); formstr <- character(nir)
            for (j in 1:nir) {
                formstr[j] <- expsum(1:((nc+1)*nir), coefs) # string of form exp(1*x1+beta1*x2*beta2*x3), exp(x4+beta*x5+...)
                parinds <- c(parinds, p.inds[j], cov.inds[(j-1)*nc + seq(length.out=nc)]) # indices into HMM par vector corresp to x1,x2,x3,...
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
    indmat[indmat == 1] <- seq(length.out = x$qmodel$npars)
    indmat <- t(indmat) # matrix of indices of estimate vector
    inds <- seq(length.out = x$qmodel$npars+x$qcmodel$npars) # identifiers for q and beta parameters
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
    indmat[indmat==1] <- seq(length.out = ni)
    indmat <- t(indmat) # matrix of indices of estimate vector
    inds <- seq(length.out = ni + ni*nc)
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
    for (i in seq(length.out=ncut)){
        covlistlist[[i+1]] <- covlist
        covlistlist[[i+1]][[ti[i]]] <- 1
    }
    covlistlist
}



#' Print a fitted msm model object
#' 
#' Print a fitted msm model object (in old format, from msm 1.3.1 and earlier)
#' 
#' See \code{\link{print.msm}} for a better and cleaner output format, and an
#' explanation of the change.
#' 
#' @param x Output from \code{\link{msm}}, representing a fitted multi-state
#' model object.
#' @param ... Other arguments to be passed to \code{\link{format}}.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{print.msm}}
#' @keywords models
#' @export
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
            print_ci(x$Qmatrices[[i]], x$QmatricesL[[i]], x$QmatricesU[[i]])
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
                print_ci(x$Ematrices[[i]], x$EmatricesL[[i]], x$EmatricesU[[i]])
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
    if (!is.null(lower)) res[,2] <- lower[keep]
    if (!is.null(upper)) res[,3] <- upper[keep]
    res[,4] <- fixed[keep]
    res
}

### Format transition intensities and their covariate effects in one tidy matrix



#' Extract msm model parameter estimates in compact format
#' 
#' Extract estimates and confidence intervals for transition intensities (or
#' misclassification probabilities), and their covariate effects, in a tidy
#' matrix format with one row per transition.  This is used by the print method
#' (\code{\link{print.msm}}) for \code{msm} objects.  Covariate effects are
#' returned as hazard or odds ratios, not on the log scale.
#' 
#' 
#' @aliases msm.form.qoutput msm.form.eoutput
#' @param x A fitted multi-state model object, as returned by
#' \code{\link{msm}}.
#' @param covariates Covariate values defining the "baseline" parameters (see
#' \code{\link{qmatrix.msm}}).
#' @param cl Width of the symmetric confidence interval to present.  Defaults
#' to 0.95.
#' @param digits Minimum number of significant digits for the formatted
#' character matrix returned as an attribute.  This is passed to
#' \code{\link{format}}. Defaults to 4.
#' @param ... Other arguments to be passed to \code{\link{format}}.
#' @return A numeric matrix with one row per transition, and one column for
#' each estimate or confidence limit.  The \code{"formatted"} attribute
#' contains the same results formatted for pretty printing.
#' \code{msm.form.qoutput} returns the transition intensities and their
#' covariates, and \code{msm.form.eoutput} returns the misclassification
#' probabilities and their covariates.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{print.msm}}
#' @keywords models
#' @export
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
    for (i in seq(length.out=x$qcmodel$ncovs)){
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

#' @rdname msm.form.qoutput
#' @export
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
    for (i in seq(length.out=x$ecmodel$ncovs)){
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



#' Print a fitted msm model object
#' 
#' Print a fitted msm model object
#' 
#' This is the new method of formatting msm objects for printing.  The old
#' method was based on printing lists of matrices. That produced a lot of
#' wasted space for parameters which were zero, and it was difficult to match
#' corresponding numbers between matrices. The new method presents all the
#' transition intensities and covariate effects as a single compact table, and
#' likewise for misclassification matrices.
#' 
#' Also in the old method, covariate effects were presented as log hazard
#' ratios or log odds ratios.  The log scale is more convenient mathematically,
#' but unnatural to interpret.  The new method presents hazard ratios for
#' covariates on transition intensities and odds ratios for misclassification
#' probabilities.
#' 
#' \code{printnew.msm} is an alias for \code{print.msm}.
#' 
#' @aliases print.msm printnew.msm
#' @param x Output from \code{\link{msm}}, representing a fitted multi-state
#' model object.
#' @param covariates Covariates for which to print ``baseline'' transition
#' intensities or misclassification probabilities. See
#' \code{\link{qmatrix.msm}} for more details.
#' @param digits Minimum number of significant digits, passed to
#' \code{\link{format}}. Defaults to 4.
#' @param ... Other arguments to be passed to \code{\link{format}}.
#' @return The object returned by \code{print.msm} is a numeric matrix with one
#' column for each estimate or confidence limit for intensities and their
#' covariates, in the same arrangement as printed, but with the underlying
#' numbers in full precision.  The results formatted for printing are stored in
#' the \code{"formatted"} attribute of the object, as a character matrix.
#' These can alternatively be produced by \code{\link{msm.form.qoutput}}, which
#' has no printing side-effect. \code{\link{msm.form.eoutput}} produces the
#' same arrangement for misclassification probabilities instead of intensities.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}, \code{\link{printold.msm}},
#' \code{\link{msm.form.qoutput}}.
#' @keywords models
#' @export
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
    invisible(ret)
}

#' @rdname print.msm 
#' @export
printnew.msm <- print.msm




#' Summarise a fitted multi-state model
#' 
#' Summary method for fitted \code{\link{msm}} models. This is simply a wrapper
#' around \code{\link{prevalence.msm}} which produces a table of observed and
#' expected state prevalences for each time, and for models with covariates,
#' \code{\link{hazard.msm}} to print hazard ratios with 95\% confidence
#' intervals for covariate effects.
#' 
#' @name summary.msm
#' @aliases summary.msm print.summary.msm
#' @param object A fitted multi-state model object, as returned by
#' \code{\link{msm}}.
#' @param hazard.scale Vector with same elements as number of covariates on
#' transition rates. Corresponds to the increase in each covariate used to
#' calculate its hazard ratio. Defaults to all 1.
#' @param ... Further arguments passed to \code{\link{prevalence.msm}}.
#' @return A list of class \code{summary.msm}, with components:
#' 
#' \item{prevalences}{Output from \code{\link{prevalence.msm}}.}
#' 
#' \item{hazard}{Output from \code{\link{hazard.msm}}.}
#' 
#' \item{hazard.scale}{Value of the \code{hazard.scale} argument.}
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}},\code{\link{prevalence.msm}},
#' \code{\link{hazard.msm}}
#' @keywords models
#' @export
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



#' Plots of multi-state models
#' 
#' This produces a plot of the expected probability of survival against time,
#' from each transient state. Survival is defined as not entering an absorbing
#' state.
#' 
#' Note that while this function is only relevant to models with absorbing
#' states, models in \pkg{msm} can have any transition structure and do not
#' necessarily have to have an absorbing state.
#' 
#' 
#' @param x Output from \code{\link{msm}}, representing a fitted multi-state
#' model object.
#' @param from States from which to consider survival. Defaults to the complete
#' set of transient states.
#' @param to Absorbing state to consider. Defaults to the highest-labelled
#' absorbing state.
#' @param range Vector of two elements, giving the range of times to plot for.
#' @param covariates Covariate values for which to evaluate the expected
#' probabilities.  This can either be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula, or a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' @param legend.pos Vector of the \eqn{x} and \eqn{y} position, respectively,
#' of the legend.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param lwd Line width. See \code{\link{par}}.
#' @param ... Other arguments to be passed to the generic \code{\link{plot}}
#' and \code{\link{lines}} functions.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}
#' @keywords models
#' @export
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
            stop("\"to\" not specified, and no absorbing state. See help(plot.msm)")
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



#' Kaplan Meier estimates of incidence
#' 
#' Compute and plot Kaplan-Meier estimates of the probability that each
#' successive state has not occurred yet.
#' 
#' If the data represent observations of the process at arbitrary times, then
#' the first occurrence of the state in the data will usually be greater than
#' the actual first transition time to that state.  Therefore the probabilities
#' plotted by \code{\link{plotprog.msm}} will be overestimates.
#' 
#' @param formula A formula giving the vectors containing the observed states
#' and the corresponding observation times. For example,
#' 
#' \code{state ~ time}
#' 
#' Observed states should be in the set \code{1, \dots{}, n}, where \code{n} is
#' the number of states.
#' @param subject Vector of subject identification numbers for the data
#' specified by \code{formula}. If missing, then all observations are assumed
#' to be on the same subject. These must be sorted so that all observations on
#' the same subject are adjacent.
#' @param data An optional data frame in which the variables represented by
#' \code{state}, \code{time} and \code{subject} can be found.
#' @param legend.pos Vector of the \eqn{x} and \eqn{y} position, respectively,
#' of the legend.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param lwd Line width. See \code{\link{par}}.
#' @param xlim x axis limits, e.g. c(0,10) for an axis ranging from 0 to 10.
#' Default is the range of observation times.
#' @param mark.time Mark the empirical survival curve at each censoring point,
#' see \code{\link{lines.survfit}}.
#' @param ... Other arguments to be passed to the \code{\link{plot}} and
#' \code{\link[survival]{lines.survfit}} functions.
#' @seealso \code{\link[survival]{survfit}},
#' \code{\link[survival]{plot.survfit}}
#' @keywords models
#' @export
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



#' Explore the likelihood surface
#' 
#' Plot the log-likelihood surface with respect to two parameters.
#' 
#' Draws a contour or perspective plot.  Useful for diagnosing irregularities
#' in the likelihood surface.  If you want to use these plots before running
#' the maximum likelihood estimation, then just run \code{msm} with all
#' estimates fixed at their initial values.
#' 
#' \code{contour.msm} just calls surface.msm with \code{type = "contour"}.
#' 
#' \code{persp.msm} just calls surface.msm with \code{type = "persp"}.
#' 
#' \code{image.msm} just calls surface.msm with \code{type = "image"}.
#' 
#' As these three functions are methods of the generic functions
#' \code{contour}, \code{persp} and \code{image}, they can be invoked as
#' \code{contour(x)}, \code{persp(x)} or \code{image(x)}, where \code{x} is a
#' fitted \code{msm} object.
#' 
#' @aliases surface.msm persp.msm contour.msm image.msm
#' @param x Output from \code{\link{msm}}, representing a fitted msm model.
#' @param params Integer vector with two elements, giving the indices of the
#' parameters to vary. All other parameters will be fixed. Defaults to
#' \code{c(1,2)}, representing the first two log transition intensities. See
#' the \code{fixedpars} argument to \code{msm} for a definition of these
#' indices.
#' @param np Number of grid points to use in each direction, by default 10.  An
#' \code{np x np} grid will be used to evaluate the likelihood surface. If 100
#' likelihood function evaluations is slow, then reduce this.
#' @param type Character string specifying the type of plot to produce.
#' \tabular{ll}{ \code{"contour"} \tab Contour plot, using the R function
#' \code{\link{contour}}. \cr \code{"filled.contour"} \tab Solid-color contour
#' plot, using the R function \code{\link{filled.contour}}. \cr \code{"persp"}
#' \tab Perspective plot, using the R function \code{\link{persp}}. \cr
#' \code{"image"} \tab Grid color plot, using the R function
#' \code{\link{image}}. \cr }
#' @param point Vector of length \code{n}, where \code{n} is the number of
#' parameters in the model, including the parameters that will be varied here.
#' This specifies the point at which to fix the likelihood.  By default, this
#' is the maximum likelihood estimates stored in the fitted model \code{x},
#' \code{x$estimates}.
#' @param xrange Range to plot for the first varied parameter.  Defaults to
#' plus and minus two standard errors, obtained from the Hessian at the maximum
#' likelihood estimate.
#' @param yrange Range to plot for the second varied parameter.  Defaults to
#' plus and minus two standard errors, obtained from the Hessian at the maximum
#' likelihood estimate.
#' @param ... Further arguments to be passed to the plotting function.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}, \code{\link{contour}},
#' \code{\link{filled.contour}}, \code{\link{persp}}, \code{\link{image}}.
#' @keywords models
#' @export
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
        p1 <- seq(pmin, pmax, length.out=np)
    }
    else p1 <- seq(xrange[1], xrange[2], length.out=np)
    if (is.null(yrange)){
        pmin <- point[i2] - 2*se[i2]
        pmax <- point[i2] + 2*se[i2]
        p2 <- seq(pmin, pmax, length.out=np)
    }
    else p2 <- seq(yrange[1], yrange[2], length.out=np)

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

#' @rdname surface.msm
#' @export 
contour.msm <- function(x, ...)
{
    surface.msm(x, type="contour",...)
}

#' @rdname surface.msm
#' @export 
persp.msm <- function(x, ...)
{
    surface.msm(x, type="persp",...)
}

#' @rdname surface.msm
#' @export 
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

#' @export
print.msm.est <- function(x, digits=NULL, ...)
{
    if (is.list(x))
        print_ci(x$estimates, x$L, x$U, x$fixed, digits=digits)
    else print(unclass(x), digits=digits)
}

## Unused, remove eventually
# nocov start
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
# nocov end

#' @export
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

print_ci <- function(x, l, u, fixed=NULL, digits=NULL){
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



#' Estimated ratio of transition intensities
#' 
#' Compute the estimate and approximate standard error of the ratio of two
#' estimated transition intensities from a fitted multi-state model at a given
#' set of covariate values.
#' 
#' 
#' For example, we might want to compute the ratio of the progression rate and
#' recovery rate for a fitted model \code{disease.msm} with a health state
#' (state 1) and a disease state (state 2).  In this case, the progression rate
#' is the (1,2) entry of the intensity matrix, and the recovery rate is the
#' (2,1) entry.  Thus to compute this ratio with covariates set to their means,
#' we call
#' 
#' \code{qratio.msm(disease.msm, c(1,2), c(2,1))} .
#' 
#' Standard errors are estimated by the delta method. Confidence limits are
#' estimated by assuming normality on the log scale.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param ind1 Pair of numbers giving the indices in the intensity matrix of
#' the numerator of the ratio, for example, \code{c(1,2)}.
#' @param ind2 Pair of numbers giving the indices in the intensity matrix of
#' the denominator of the ratio, for example, \code{c(2,1)}.
#' @param covariates The covariate values at which to estimate the intensities.
#' This can either be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula, or a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' @param ci If \code{"delta"} (the default) then confidence intervals are
#' calculated by the delta method.
#' 
#' If \code{"normal"}, then calculate a confidence interval by simulating
#' \code{B} random vectors from the asymptotic multivariate normal distribution
#' implied by the maximum likelihood estimates (and covariance matrix) of the
#' log transition intensities and covariate effects, then transforming.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' @param cl Width of the symmetric confidence interval to present.  Defaults
#' to 0.95.
#' @param B Number of bootstrap replicates, or number of normal simulations
#' from the distribution of the MLEs
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @return A named vector with elements \code{estimate}, \code{se}, \code{L}
#' and \code{U} containing the estimate, standard error, lower and upper
#' confidence limits, respectively, of the ratio of intensities.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{qmatrix.msm}}
#' @keywords models
#' @export
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



#' Transition probability matrix
#' 
#' Extract the estimated transition probability matrix from a fitted
#' continuous-time multi-state model for a given time interval, at a given set
#' of covariate values.
#' 
#' For a continuous-time homogeneous Markov process with transition intensity
#' matrix \eqn{Q}, the probability of occupying state \eqn{s} at time \eqn{u +
#' t} conditionally on occupying state \eqn{r} at time \eqn{u} is given by the
#' \eqn{(r,s)} entry of the matrix \eqn{P(t) = \exp(tQ)}{P(t) = exp(tQ)}, where
#' \eqn{\exp()}{exp()} is the matrix exponential.
#' 
#' For non-homogeneous processes, where covariates and hence the transition
#' intensity matrix \eqn{Q} are piecewise-constant in time, the transition
#' probability matrix is calculated as a product of matrices over a series of
#' intervals, as explained in \code{\link{pmatrix.piecewise.msm}}.
#' 
#' The \code{\link{pmatrix.piecewise.msm}} function is only necessary for
#' models fitted using a time-dependent covariate in the \code{covariates}
#' argument to \code{\link{msm}}. For time-inhomogeneous models fitted using
#' "pci", \code{pmatrix.msm} can be used, with arguments \code{t} and
#' \code{t1}, to calculate transition probabilities over any time period.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param t The time interval to estimate the transition probabilities for, by
#' default one unit.
#' @param t1 The starting time of the interval. Used for models \code{x} with
#' piecewise-constant intensities fitted using the \code{pci} option to
#' \code{\link{msm}}. The probabilities will be computed on the interval [t1,
#' t1+t].
#' @param covariates The covariate values at which to estimate the transition
#' probabilities.  This can either be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula, or a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' 
#' If some covariates are specified but not others, the missing ones default to
#' zero.
#' 
#' For time-inhomogeneous models fitted using the \code{pci} option to
#' \code{\link{msm}}, "covariates" here include only those specified using the
#' \code{covariates} argument to \code{\link{msm}}, and exclude the artificial
#' covariates representing the time period.
#' 
#' For time-inhomogeneous models fitted "by hand" by using a time-dependent
#' covariate in the \code{covariates} argument to \code{\link{msm}}, the
#' function \code{\link{pmatrix.piecewise.msm}} should be used to to calculate
#' transition probabilities.
#' @param ci If \code{"normal"}, then calculate a confidence interval for the
#' transition probabilities by simulating \code{B} random vectors from the
#' asymptotic multivariate normal distribution implied by the maximum
#' likelihood estimates (and covariance matrix) of the log transition
#' intensities and covariate effects, then calculating the resulting transition
#' probability matrix for each replicate. See, e.g. Mandel (2013) for a
#' discussion of this approach.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' 
#' If \code{"none"} (the default) then no confidence interval is calculated.
#' @param cl Width of the symmetric confidence interval, relative to 1.
#' @param B Number of bootstrap replicates, or number of normal simulations
#' from the distribution of the MLEs
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @param qmatrix A transition intensity matrix.  Either this or a fitted model
#' \code{x} must be supplied.  No confidence intervals are available if
#' \code{qmatrix} is supplied.
#' @param ... Optional arguments to be passed to \code{\link{MatrixExp}} to
#' control the method of computing the matrix exponential.
#' @return The matrix of estimated transition probabilities \eqn{P(t)} in the
#' given time.  Rows correspond to "from-state" and columns to "to-state".
#' 
#' Or if \code{ci="normal"} or \code{ci="bootstrap"}, \code{pmatrix.msm}
#' returns a list with components \code{estimates} and \code{ci}, where
#' \code{estimates} is the matrix of estimated transition probabilities, and
#' \code{ci} is a list of two matrices containing the upper and lower
#' confidence limits.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}.
#' @seealso \code{\link{qmatrix.msm}}, \code{\link{pmatrix.piecewise.msm}},
#' \code{\link{boot.msm}}
#' @references Mandel, M. (2013). "Simulation based confidence intervals for
#' functions with complicated derivatives." The American Statistician
#' 67(2):76-81
#' @keywords models
#' @export
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



#' Transition probability matrix for processes with piecewise-constant
#' intensities
#' 
#' Extract the estimated transition probability matrix from a fitted
#' non-time-homogeneous multi-state model for a given time interval.  This is a
#' generalisation of \code{\link{pmatrix.msm}} to models with time-dependent
#' covariates.  Note that \code{\link{pmatrix.msm}} is sufficient to calculate
#' transition probabilities for time-inhomogeneous models fitted using the
#' \code{pci} argument to \code{\link{msm}}.
#' 
#' Suppose a multi-state model has been fitted, in which the transition
#' intensity matrix \eqn{Q(x(t))} is modelled in terms of time-dependent
#' covariates \eqn{x(t)}.  The transition probability matrix \eqn{P(t_1,
#' t_n)}{P(t1, tn)} for the time interval \eqn{(t_1, }{(t1, tn)}\eqn{
#' t_n)}{(t1, tn)} cannot be calculated from the estimated intensity matrix as
#' \eqn{\exp((t_n - t_1) Q)}{exp((tn - t1) Q)}, because \eqn{Q} varies within
#' the interval \eqn{t_1, t_n}{t1, tn}.  However, if the covariates are
#' piecewise-constant, or can be approximated as piecewise-constant, then we
#' can calculate \eqn{P(t_1, t_n)}{P(t1, tn)} by multiplying together
#' individual matrices \eqn{P(t_i, }{P(t_i, t_{i+1}) = exp((t_{i+1} - t_i)
#' Q)}\eqn{ t_{i+1}) = \exp((t_{i+1} - t_i) Q)}{P(t_i, t_{i+1}) = exp((t_{i+1}
#' - t_i) Q)}, calculated over intervals where Q is constant:
#' 
#' \deqn{P(t_1, t_n) = P(t_1, t_2) P(t_2, t_3)\ldots P(t_{n-1}, }{P(t1, tn) =
#' P(t1, t2) P(t2, t3)\ldotsP(tn-1, tn)}\deqn{ t_n)}{P(t1, tn) = P(t1, t2)
#' P(t2, t3)\ldotsP(tn-1, tn)}
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}. This
#' should be a non-homogeneous model, whose transition intensity matrix depends
#' on a time-dependent covariate.
#' @param t1 The start of the time interval to estimate the transition
#' probabilities for.
#' @param t2 The end of the time interval to estimate the transition
#' probabilities for.
#' @param times Cut points at which the transition intensity matrix changes.
#' @param covariates A list with number of components one greater than the
#' length of \code{times}.  Each component of the list is specified in the same
#' way as the \code{covariates} argument to \code{\link{pmatrix.msm}}.  The
#' components correspond to the covariate values in the intervals
#' 
#' \code{(t1, times[1]], (times[1], times[2]], ..., (times[length(times)], t2]}
#' 
#' (assuming that all elements of \code{times} are in the interval \code{(t1,
#' t2)}).
#' @param ci If \code{"normal"}, then calculate a confidence interval for the
#' transition probabilities by simulating \code{B} random vectors from the
#' asymptotic multivariate normal distribution implied by the maximum
#' likelihood estimates (and covariance matrix) of the log transition
#' intensities and covariate effects, then calculating the resulting transition
#' probability matrix for each replicate.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' 
#' If \code{"none"} (the default) then no confidence interval is calculated.
#' @param cl Width of the symmetric confidence interval, relative to 1.
#' @param B Number of bootstrap replicates, or number of normal simulations
#' from the distribution of the MLEs
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @param qlist A list of transition intensity matrices, of length one greater
#' than the length of \code{times}.  Either this or a fitted model \code{x}
#' must be supplied.  No confidence intervals are available if (just)
#' \code{qlist} is supplied.
#' @param ... Optional arguments to be passed to \code{\link{MatrixExp}} to
#' control the method of computing the matrix exponential.
#' @return The matrix of estimated transition probabilities \eqn{P(t)} for the
#' time interval \code{[t1, tn]}.  That is, the probabilities of occupying
#' state \eqn{s} at time \eqn{t_n}{tn} conditionally on occupying state \eqn{r}
#' at time \eqn{t_1}{t1}.  Rows correspond to "from-state" and columns to
#' "to-state".
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{pmatrix.msm}}
#' @keywords models
#' @examples
#' 
#' \dontrun{
#' ## In a clinical study, suppose patients are given a placebo in the
#' ## first 5 weeks, then they begin treatment 1 at 5 weeks, and
#' ## a combination of treatments 1 and 2 from 10 weeks.
#' ## Suppose a multi-state model x has been fitted for the patients'
#' ## progress, with treat1 and treat2 as time dependent covariates.
#' 
#' ## Cut points for when treatment covariate changes
#' times <- c(0, 5, 10)
#' 
#' ## Indicators for which treatments are active in the four intervals
#' ## defined by the three cut points
#' covariates <- list( list (treat1=0, treat2=0), list (treat1=0, treat2=0), list(treat1=1, treat2=0),
#' list(treat1=1, treat2=1) )
#' 
#' ## Calculate transition probabilities from the start of the study to 15 weeks
#' pmatrix.piecewise.msm(x, 0, 15, times, covariates)
#' }
#' 
#' @export
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
        P.middle <- if (!is.null(x)) diag(x$qmodel$nstates) else diag(ncol(qlist[[1]]))
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



#' Mean sojourn times from a multi-state model
#' 
#' Estimate the mean sojourn times in the transient states of a multi-state
#' model and their confidence limits.
#' 
#' The mean sojourn time in a transient state \eqn{r} is estimated by \eqn{- 1
#' / q_{rr}}, where \eqn{q_{rr}} is the \eqn{r}th entry on the diagonal of the
#' estimated transition intensity matrix.
#' 
#' A continuous-time Markov model is fully specified by the mean sojourn times
#' and the probability that each state is next (\code{\link{pnext.msm}}).  This
#' is a more intuitively meaningful description of a model than the transition
#' intensity matrix (\code{\link{qmatrix.msm}}).
#' 
#' Time dependent covariates, or time-inhomogeneous models, are not supported.
#' This would require the mean of a piecewise exponential distribution, and the
#' package author is not aware of any general analytic form for that.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param covariates The covariate values at which to estimate the mean sojourn
#' times. This can either be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' a list of values, with optional names. For example,
#' 
#' \code{list(60, 1)}, where the order of the list follows the order of the
#' covariates originally given in the model formula, or a named list, e.g.
#' 
#' \code{list (age = 60, sex = 1)}
#' 
#' @param ci If \code{"delta"} (the default) then confidence intervals are
#' calculated by the delta method, or by simple transformation of the Hessian
#' in the very simplest cases.
#' 
#' If \code{"normal"}, then calculate a confidence interval by simulating
#' \code{B} random vectors from the asymptotic multivariate normal distribution
#' implied by the maximum likelihood estimates (and covariance matrix) of the
#' log transition intensities and covariate effects, then transforming.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' @param cl Width of the symmetric confidence interval to present.  Defaults
#' to 0.95.
#' @param B Number of bootstrap replicates, or number of normal simulations
#' from the distribution of the MLEs
#' @return A data frame with components:
#' 
#' \item{estimates}{Estimated mean sojourn times in the transient states.}
#' \item{SE}{Corresponding standard errors.} \item{L}{Lower confidence limits.}
#' \item{U}{Upper confidence limits.}
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}, \code{\link{qmatrix.msm}},
#' \code{\link{deltamethod}}
#' @keywords models
#' @export
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



#' Probability of each state being next
#' 
#' Compute a matrix of the probability of each state \eqn{s} being the next
#' state of the process after each state \eqn{r}.  Together with the mean
#' sojourn times in each state (\code{\link{sojourn.msm}}), these fully define
#' a continuous-time Markov model.
#' 
#' For a continuous-time Markov process in state \eqn{r}, the probability that
#' the next state is \eqn{s} is \eqn{-q_{rs} / q_{rr}}, where \eqn{q_{rs}} is
#' the transition intensity (\code{\link{qmatrix.msm}}).
#' 
#' A continuous-time Markov model is fully specified by these probabilities
#' together with the mean sojourn times \eqn{-1/q_{rr}} in each state \eqn{r}.
#' This gives a more intuitively meaningful description of a model than the
#' intensity matrix.
#' 
#' Remember that \pkg{msm} deals with continuous-time, not discrete-time
#' models, so these are \emph{not} the same as the probability of observing
#' state \eqn{s} at a fixed time in the future.  Those probabilities are given
#' by \code{\link{pmatrix.msm}}.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param covariates The covariate values at which to estimate the intensities.
#' This can either be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula, or a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' @param ci If \code{"normal"} (the default) then calculate a confidence
#' interval by simulating \code{B} random vectors from the asymptotic
#' multivariate normal distribution implied by the maximum likelihood estimates
#' (and covariance matrix) of the log transition intensities and covariate
#' effects, then transforming.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' 
#' If \code{"delta"} then confidence intervals are calculated based on the
#' delta method SEs of the log rates, but this is not recommended since it may
#' not respect the constraint that probabilities are less than one.
#' @param cl Width of the symmetric confidence interval to present.  Defaults
#' to 0.95.
#' @param B Number of bootstrap replicates, or number of normal simulations
#' from the distribution of the MLEs.
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @return The matrix of probabilities that the next move of a process in state
#' \eqn{r} (rows) is to state \eqn{s} (columns).
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso
#' \code{\link{qmatrix.msm}},\code{\link{pmatrix.msm}},\code{\link{qratio.msm}}
#' @keywords models
#' @export
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



#' Extract model coefficients
#' 
#' Extract the estimated log transition intensities and the corresponding
#' linear effects of each covariate.
#' 
#' 
#' @param object A fitted multi-state model object, as returned by
#' \code{\link{msm}}.
#' @param ... (unused) further arguments passed to or from other methods.
#' @return If there is no misclassification, \code{coef.msm} returns a list of
#' matrices.  The first component, labelled \code{logbaseline}, is a matrix
#' containing the estimated transition intensities on the log scale with any
#' covariates fixed at their means in the data. Each remaining component is a
#' matrix giving the linear effects of the labelled covariate on the matrix of
#' log intensities. \cr
#' 
#' For misclassification models, \code{coef.msm} returns a list of lists. The
#' first component, \code{Qmatrices}, is a list of matrices as described in the
#' previous paragraph.  The additional component \code{Ematrices} is a list of
#' similar format containing the logit-misclassification probabilities and any
#' estimated covariate effects.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}
#' @keywords models
#' @export
coef.msm <- function(object, ...)
{
    if (!inherits(object, "msm")) stop("expected object to be a msm model")
    if (object$emodel$misc)
        object[c("Qmatrices", "Ematrices")]
    else object$Qmatrices
}

### Extract the log-likelihood



#' Extract model log-likelihood
#' 
#' Extract the log-likelihood and the number of parameters of a model fitted
#' with \code{\link{msm}}.
#' 
#' 
#' @param object A fitted multi-state model object, as returned by
#' \code{\link{msm}}.
#' @param by.subject Return vector of subject-specific log-likelihoods, which
#' should sum to the total log-likelihood.
#' @param ... (unused) further arguments passed to or from other methods.
#' @return The log-likelihood of the model represented by 'object' evaluated at
#' the maximum likelihood estimates.
#' 
#' Akaike's information criterion can also be computed using
#' \code{\link{AIC}(object)}.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}},\code{\link{lrtest.msm}}.
#' @keywords models
#' @export
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



#' Likelihood ratio test
#' 
#' Likelihood ratio test between two or more fitted multi-state models
#' 
#' 
#' @param ... Two or more fitted multi-state models, as returned by
#' \code{\link{msm}}, ordered by increasing numbers of parameters.
#' @return A matrix with three columns, giving the likelihood ratio statistic,
#' difference in degrees of freedom and the chi-squared p-value for a
#' comparison of the first model supplied with each subsequent model.
#' @section Warning: The comparison between models will only be valid if they
#' are fitted to the same dataset. This may be a problem if there are missing
#' values and R's default of 'na.action = na.omit' is used.
#' 
#' The likelihood ratio statistic only has the indicated chi-squared
#' distribution if the models are nested. An alternative for comparing
#' non-nested models is Akaike's information criterion.  This can be computed
#' for one or more fitted \code{msm} models \code{x,y,...} using
#' \code{\link{AIC}(x,y,...)}.
#' @seealso \code{\link{logLik.msm}},\code{\link{msm}}
#' @keywords models
#' @export lrtest.msm
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



#' Total length of stay, or expected number of visits
#' 
#' Estimate the expected total length of stay, or the expected number of
#' visits, in each state, for an individual in a given period of evolution of a
#' multi-state model.
#' 
#' The expected total length of stay in state \eqn{j} between times \eqn{t_1}
#' and \eqn{t_2}, from the point of view of an individual in state \eqn{i} at
#' time 0, is defined by the integral from \eqn{t_1} to \eqn{t_2} of the
#' \eqn{i,j} entry of the transition probability matrix \eqn{P(t) = Exp(tQ)},
#' where \eqn{Q} is the transition intensity matrix.
#' 
#' The corresponding expected number of visits to state \eqn{j} (excluding the
#' stay in the current state at time 0) is \eqn{\sum_{i!=j} T_i Q_{i,j}}, where
#' \eqn{T_i} is the expected amount of time spent in state \eqn{i}.
#' 
#' More generally, suppose that \eqn{\pi_0}{pi_0} is the vector of
#' probabilities of being in each state at time 0, supplied in \code{start},
#' and we want the vector \eqn{\mathbf{x}}{x} giving the expected lengths of
#' stay in each state.  The corresponding integral has the following solution
#' (van Loan 1978; van Rosmalen et al. 2013)
#' 
#' \deqn{\mathbf{x} =
#' \left[
#' \begin{array}{ll}
#' 1  &  \mathbf{0}_K
#' \end{array}
#' \right]
#' Exp(t Q')
#' \left[
#' \begin{array}{l} \mathbf{0}_K\\I_K
#' \end{array}
#' \right]
#' }{x = [1, 0_K] Exp(t Q') [0_K, I_K]'}
#'
#' where \deqn{Q' = \left[
#' \begin{array}{ll} 0 & \mathbf{\pi}_0\\
#' \mathbf{0}_K  &  Q - rI_K
#' \end{array}
#' \right]
#' }{Q' = rbind(c(0, pi_0), cbind(0_K, Q - r I_K))}
#' 
#' \eqn{\pi_0}{pi_0} is the row vector of initial state probabilities supplied
#' in \code{start}, \eqn{\mathbf{0}_K}{0_K} is the row vector of K zeros,
#' \eqn{r} is the discount rate, \eqn{I_K}{I_K} is the K x K identity matrix,
#' and \eqn{Exp} is the matrix exponential.
#' 
#' Alternatively, the integrals can be calculated numerically, using the
#' \code{\link{integrate}} function.  This may take a long time for models with
#' many states where \eqn{P(t)} is expensive to calculate.  This is required
#' where \code{tot = Inf}, since the package author is not aware of any
#' analytic expression for the limit of the above formula as \eqn{t} goes to
#' infinity.
#' 
#' With the argument \code{num.integ=TRUE}, numerical integration is used even
#' where the analytic solution is available. This facility is just provided for
#' checking results against versions 1.2.4 and earlier, and will be removed
#' eventually. Please let the package maintainer know if any results are
#' different.
#' 
#' For a model where the individual has only one place to go from each state,
#' and each state is visited only once, for example a progressive disease model
#' with no recovery or death, these are equal to the mean sojourn time in each
#' state.  However, consider a three-state health-disease-death model with
#' transitions from health to disease, health to death, and disease to death,
#' where everybody starts healthy.  In this case the mean sojourn time in the
#' disease state will be greater than the expected length of stay in the
#' disease state.  This is because the mean sojourn time in a state is
#' conditional on entering the state, whereas the expected total time diseased
#' is a forecast for a healthy individual, who may die before getting the
#' disease.
#' 
#' In the above formulae, \eqn{Q} is assumed to be constant over time, but the
#' results generalise easily to piecewise-constant intensities.  This function
#' automatically handles models fitted using the \code{pci} option to
#' \code{\link{msm}}. For any other inhomogeneous models, the user must specify
#' \code{piecewise.times} and \code{piecewise.covariates} arguments to
#' \code{\link{totlos.msm}}.
#' 
#' @aliases totlos.msm envisits.msm
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param start Either a single number giving the state at the beginning of the
#' period, or a vector of probabilities of being in each state at this time.
#' @param end States to estimate the total length of stay (or number of visits)
#' in. Defaults to all states.  This is deprecated, since with the analytic
#' solution (see "Details") it doesn't save any computation to only estimate
#' for a subset of states.
#' @param fromt Time from which to estimate.  Defaults to 0, the beginning of
#' the process.
#' @param tot Time up to which the estimate is made.  Defaults to infinity,
#' giving the expected time spent in or number of visits to the state until
#' absorption. However, the calculation will be much more efficient if a finite
#' (potentially large) time is specified: see the "Details" section.  For
#' models without an absorbing state, \code{t} must be specified.
#' @param covariates The covariate values to estimate for.  This can either
#' be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula, or a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' 
#' @param piecewise.times Times at which piecewise-constant intensities change.
#' See \code{\link{pmatrix.piecewise.msm}} for how to specify this. This is
#' only required for time-inhomogeneous models specified using explicit
#' time-dependent covariates, and should not be used for models specified using
#' "pci".
#' @param piecewise.covariates Covariates on which the piecewise-constant
#' intensities depend. See \code{\link{pmatrix.piecewise.msm}} for how to
#' specify this.
#' @param num.integ Use numerical integration instead of analytic solution (see
#' below).
#' @param discount Discount rate in continuous time.
#' @param env Supplied to \code{\link{totlos.msm}}.  If \code{TRUE}, return the
#' expected number of visits to each state. If \code{FALSE}, return the total
#' length of stay in each state. \code{\link{envisits.msm}} simply calls
#' \code{\link{totlos.msm}} with \code{env=TRUE}.
#' @param ci If \code{"normal"}, then calculate a confidence interval by
#' simulating \code{B} random vectors from the asymptotic multivariate normal
#' distribution implied by the maximum likelihood estimates (and covariance
#' matrix) of the log transition intensities and covariate effects, then
#' calculating the total length of stay for each replicate.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' 
#' If \code{"none"} (the default) then no confidence interval is calculated.
#' @param cl Width of the symmetric confidence interval, relative to 1
#' @param B Number of bootstrap replicates
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @param ... Further arguments to be passed to the \code{\link{integrate}}
#' function to control the numerical integration.
#' @return A vector of expected total lengths of stay
#' (\code{\link{totlos.msm}}), or expected number of visits
#' (\code{\link{envisits.msm}}), for each transient state.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{sojourn.msm}}, \code{\link{pmatrix.msm}},
#' \code{\link{integrate}}, \code{\link{boot.msm}}.
#' @references C. van Loan (1978). Computing integrals involving the matrix
#' exponential. IEEE Transactions on Automatic Control 23(3)395-404.
#' 
#' J. van Rosmalen, M. Toy and J.F. O'Mahony (2013). A mathematical approach
#' for evaluating Markov models in continuous time without discrete-event
#' simulation.  Medical Decision Making 33:767-779.
#' @keywords models
#' @export totlos.msm
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

    if (!is.null(x$pci)){
        piecewise.times <- x$pci
        piecewise.covariates <- msm.fill.pci.covs(x, covariates)
    }
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

#' @rdname totlos.msm
#' @export
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



#' Transient and absorbing states
#' 
#' Returns the transient and absorbing states of either a fitted model or a
#' transition intensity matrix.
#' 
#' 
#' @aliases transient.msm absorbing.msm
#' @param x A fitted multi-state model as returned by \code{\link{msm}}.
#' @param qmatrix A transition intensity matrix. The diagonal is ignored and
#' taken to be minus the sum of the rest of the row.
#' @return A vector of the ordinal indices of the transient or absorbing
#' states.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @keywords models
#' @export transient.msm
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

#' @rdname transient.msm
#' @export
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
    obstab <- t(apply(states.expand, 2, function(y) table(factor(y, levels=seq(length.out=x$qmodel$nstates)))))
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
                        function(y) table(factor(y, levels=seq(length.out=x$qmodel$nstates)))))
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
            if (x$qcmodel$ncovs>0 && isTRUE(covariates=="population")) {
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



#' Tables of observed and expected prevalences
#' 
#' This provides a rough indication of the goodness of fit of a multi-state
#' model, by estimating the observed numbers of individuals occupying each
#' state at a series of times, and comparing these with forecasts from the
#' fitted model.
#' 
#' The fitted transition probability matrix is used to forecast expected
#' prevalences from the state occupancy at the initial time.  To produce the
#' expected number in state \eqn{j} at time \eqn{t} after the start, the number
#' of individuals under observation at time \eqn{t} (including those who have
#' died, but not those lost to follow-up) is multiplied by the product of the
#' proportion of individuals in each state at the initial time and the
#' transition probability matrix in the time interval \eqn{t}.  The proportion
#' of individuals in each state at the "initial" time is estimated, if
#' necessary, in the same way as the observed prevalences.
#' 
#' For misclassification models (fitted using an \code{ematrix}), this aims to
#' assess the fit of the full model for the \emph{observed} states.  That is,
#' the combined Markov progression model for the true states and the
#' misclassification model. Thus, expected prevalences of \emph{true} states
#' are estimated from the assumed proportion occupying each state at the
#' initial time using the fitted transition probabiliy matrix. The vector of
#' expected prevalences of true states is then multiplied by the fitted
#' misclassification probability matrix to obtain the expected prevalences of
#' \emph{observed} states.
#' 
#' For general hidden Markov models, the observed state is taken to be the
#' predicted underlying state from the Viterbi algorithm
#' (\code{\link{viterbi.msm}}).  The goodness of fit of these states to the
#' underlying Markov model is tested.
#' 
#' In any model, if there are censored states, then these are replaced by
#' imputed values of highest probability from the Viterbi algorithm in order to
#' calculate the observed state prevalences.
#' 
#' For an example of this approach, see Gentleman \emph{et al.} (1994).
#' 
#' @param x A fitted multi-state model produced by \code{\link{msm}}.
#' @param times Series of times at which to compute the observed and expected
#' prevalences of states.
#' @param timezero Initial time of the Markov process. Expected values are
#' forecasted from here. Defaults to the minimum of the observation times given
#' in the data.
#' @param initstates Optional vector of the same length as the number of
#' states. Gives the numbers of individuals occupying each state at the initial
#' time, to be used for forecasting expected prevalences.  The default is those
#' observed in the data.  These should add up to the actual number of people in
#' the study at the start.
#' @param covariates Covariate values for which to forecast expected state
#' occupancy.  With the default \code{covariates="population"}, expected
#' prevalences are produced by summing model predictions over the covariates
#' observed in the original data, for a fair comparison with the observed
#' prevalences.  This may be slow, particularly with continuous covariates.
#' 
#' Predictions for fixed covariates can be obtained by supplying covariate
#' values in the standard way, as in \code{\link{qmatrix.msm}}. Therefore if
#' \code{covariates="population"} is too slow, using the mean observed values
#' through \code{covariates="mean"} may give a reasonable approximation.
#' 
#' This argument is ignored if \code{piecewise.times} is specified. If there
#' are a mixture of time-constant and time-dependent covariates, then the
#' values for all covariates should be supplied in \code{piecewise.covariates}.
#' @param misccovariates (Misclassification models only) Values of covariates
#' on the misclassification probability matrix for converting expected true to
#' expected misclassified states.  Ignored if \code{covariates="population"},
#' otherwise defaults to the mean values of the covariates in the data set.
#' @param piecewise.times Times at which piecewise-constant intensities change.
#' See \code{\link{pmatrix.piecewise.msm}} for how to specify this.  Ignored if
#' \code{covariates="population"}.  This is only required for
#' time-inhomogeneous models specified using explicit time-dependent
#' covariates, and should not be used for models specified using "pci".
#' @param piecewise.covariates Covariates on which the piecewise-constant
#' intensities depend. See \code{\link{pmatrix.piecewise.msm}} for how to
#' specify this. Ignored if \code{covariates="population"}.
#' @param ci If \code{"normal"}, then calculate a confidence interval for the
#' expected prevalences by simulating \code{B} random vectors from the
#' asymptotic multivariate normal distribution implied by the maximum
#' likelihood estimates (and covariance matrix) of the log transition
#' intensities and covariate effects, then calculating the expected prevalences
#' for each replicate.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' 
#' If \code{"none"} (the default) then no confidence interval is calculated.
#' @param cl Width of the symmetric confidence interval, relative to 1
#' @param B Number of bootstrap replicates
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @param interp Suppose an individual was observed in states \eqn{S_{r-1}} and
#' \eqn{S_r} at two consecutive times \eqn{t_{r-1}} and \eqn{t_r}, and we want
#' to estimate 'observed' prevalences at a time \eqn{t} between \eqn{t_{r-1}}
#' and \eqn{t_r}.
#' 
#' If \code{interp="start"}, then individuals are assumed to be in state
#' \eqn{S_{r-1}} at time \eqn{t}, the same state as they were at \eqn{t_{r-1}}.
#' 
#' If \code{interp="midpoint"} then if \eqn{t <= (t_{r-1} + t_r) / 2}, the
#' midpoint of \eqn{t_{r-1}} and \eqn{t_r}, the state at \eqn{t} is assumed to
#' be \eqn{S_{r-1}}, otherwise \eqn{S_{r}}. This is generally more reasonable
#' for "progressive" models.
#' @param censtime Adjustment to the observed prevalences to account for
#' limited follow-up in the data.
#' 
#' If the time is greater than \code{censtime} and the patient has reached an
#' absorbing state, then that subject will be removed from the risk set.  For
#' example, if patients have died but would only have been observed up to this
#' time, then this avoids overestimating the proportion of people who are dead
#' at later times.
#' 
#' This can be supplied as a single value, or as a vector with one element per
#' subject (after any \code{subset} has been taken), in the same order as the
#' original data.  This vector also only includes subjects with complete data,
#' thus it excludes for example subjects with only one observation (thus no
#' observed transitions), and subjects for whom every observation has missing
#' values.  (Note, to help construct this, the complete data used for the model
#' fit can be accessed with \code{model.frame(x)}, where \code{x} is the fitted
#' model object)
#' 
#' This is ignored if it is less than the subject's maximum observation time.
#' @param subset Subset of subjects to calculate observed prevalences for.
#' @param plot Generate a plot of observed against expected prevalences. See
#' \code{\link{plot.prevalence.msm}}
#' @param ... Further arguments to pass to \code{\link{plot.prevalence.msm}}.
#' @return A list of matrices, with components:
#' 
#' \item{Observed}{Table of observed numbers of individuals in each state at
#' each time}
#' 
#' \item{Observed percentages}{Corresponding percentage of the individuals at
#' risk at each time.}
#' 
#' \item{Expected}{Table of corresponding expected numbers.}
#' 
#' \item{Expected percentages}{Corresponding percentage of the individuals at
#' risk at each time.}
#' 
#' Or if \code{ci.boot = TRUE}, the component \code{Expected} is a list with
#' components \code{estimates} and \code{ci}.\cr \code{estimates} is a matrix
#' of the expected prevalences, and \code{ci} is a list of two matrices,
#' containing the confidence limits. The component \code{Expected percentages}
#' has a similar format.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}, \code{\link{summary.msm}}
#' @references Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P.
#' Multi-state Markov models for analysing incomplete disease history data with
#' illustrations for HIV disease.  \emph{Statistics in Medicine} (1994) 13(3):
#' 805--821.
#' 
#' Titman, A.C., Sharples, L. D.  Model diagnostics for multi-state models.
#' \emph{Statistical Methods in Medical Research} (2010) 19(6):621-651.
#' @keywords models
#' @export prevalence.msm
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



#' Plot of observed and expected prevalences
#' 
#' Provides a rough indication of goodness of fit of a multi-state model, by
#' estimating the observed numbers of individuals occupying a state at a series
#' of times, and plotting these against forecasts from the fitted model, for
#' each state.  Observed prevalences are indicated as solid lines, expected
#' prevalences as dashed lines.
#' 
#' See \code{\link{prevalence.msm}} for details of the assumptions underlying
#' this method.
#' 
#' Observed prevalences are plotted with a solid line, and expected prevalences
#' with a dotted line.
#' 
#' @param x A fitted multi-state model produced by \code{\link{msm}}.
#' @param mintime Minimum time at which to compute the observed and expected
#' prevalences of states.
#' @param maxtime Maximum time at which to compute the observed and expected
#' prevalences of states.
#' @param timezero Initial time of the Markov process. Expected values are
#' forecasted from here. Defaults to the minimum of the observation times given
#' in the data.
#' @param initstates Optional vector of the same length as the number of
#' states. Gives the numbers of individuals occupying each state at the initial
#' time, to be used for forecasting expected prevalences.  The default is those
#' observed in the data.  These should add up to the actual number of people in
#' the study at the start.
#' @param interp Interpolation method for observed states, see
#' \code{\link{prevalence.msm}}.
#' @param censtime Subject-specific maximum follow-up times, see
#' \code{\link{prevalence.msm}}.
#' @param subset Vector of the subject identifiers to calculated observed
#' prevalences for.
#' @param covariates Covariate values for which to forecast expected state
#' occupancy.  See \code{\link{prevalence.msm}} --- if this function runs too
#' slowly, as it may if there are continuous covariates, replace
#' \code{covariates="population"} with \code{covariates="mean"}.
#' @param misccovariates (Misclassification models only) Values of covariates
#' on the misclassification probability matrix. See
#' \code{\link{prevalence.msm}}.
#' @param piecewise.times Times at which piecewise-constant intensities change.
#' See \code{\link{prevalence.msm}}.
#' @param piecewise.covariates Covariates on which the piecewise-constant
#' intensities depend. See \code{\link{prevalence.msm}}.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param lwd.obs Line width for observed prevalences. See \code{\link{par}}.
#' @param lwd.exp Line width for expected prevalences. See \code{\link{par}}.
#' @param lty.obs Line type for observed prevalences. See \code{\link{par}}.
#' @param lty.exp Line type for expected prevalences. See \code{\link{par}}.
#' @param col.obs Line colour for observed prevalences. See \code{\link{par}}.
#' @param col.exp Line colour for expected prevalences. See \code{\link{par}}.
#' @param legend.pos Vector of the \eqn{x} and \eqn{y} position, respectively,
#' of the legend.
#' @param ... Further arguments to be passed to the generic \code{\link{plot}}
#' function.
#' @seealso \code{\link{prevalence.msm}}
#' @references Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P.
#' Multi-state Markov models for analysing incomplete disease history data with
#' illustrations for HIV disease.  \emph{Statistics in Medicine} (1994) 13(3):
#' 805--821.
#' @keywords models
#' @export plot.prevalence.msm
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
    t <- seq(mintime, maxtime, length.out=100)
    obs <- observed.msm(x, t, interp, censtime, subset)
    expec <- expected.msm(x, t, timezero=timezero, initstates=initstates, covariates=covariates, misccovariates=misccovariates,
                          piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates, risk=obs$risk, subset=subset, ci="none")[[2]]
    states <- seq(length.out=x$qmodel$nstates)
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



#' Plot empirical and fitted survival curves
#' 
#' Plot a Kaplan-Meier estimate of the survival probability and compare it with
#' the fitted survival probability from a \code{msm} model.
#' 
#' If the data represent observations of the process at arbitrary times, then
#' the first occurrence of the absorbing state in the data will usually be
#' greater than the actual first transition time to that state.  Therefore the
#' Kaplan-Meier estimate of the survival probability will be an overestimate.
#' 
#' The method of Turnbull (1976) could be used to give a non-parametric
#' estimate of the time to an interval-censored event, and compared to the
#' equivalent estimate from a multi-state model.  This is implemented in the
#' CRAN package \pkg{interval} (Fay and Shaw 2010).
#' 
#' This currently only handles time-homogeneous models.
#' 
#' @param x Output from \code{\link{msm}}, representing a fitted multi-state
#' model object.
#' @param from Non-absorbing state from which to consider survival.  Defaults
#' to state 1.  The fitted probabilities will then be calculated as the
#' transition probabilities from this state to \code{to}.  The empirical
#' survival curve plots survival from the first observation of \code{from}
#' (where this exists) to the first entry time into \code{to}.
#' @param to Absorbing state to consider. Defaults to the highest-labelled
#' absorbing state.
#' @param range Vector of two elements, giving the range of times to plot for.
#' @param covariates Covariate values for which to evaluate the expected
#' probabilities.  This can either be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula, or a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' 
#' but note the empirical curve is plotted for the full population.  To
#' consider subsets for the empirical curve, set \code{survdata=TRUE} to
#' extract the survival data and build a survival plot by hand using
#' \code{\link[survival]{plot.survfit}}.
#' @param ci If \code{"none"} (the default) no confidence intervals are
#' plotted.  If \code{"normal"} or \code{"bootstrap"}, confidence intervals are
#' plotted based on the respective method in \code{\link{pmatrix.msm}}. This is
#' very computationally-intensive, since intervals must be computed at a series
#' of times.
#' @param B Number of bootstrap or normal replicates for the confidence
#' interval.  The default is 100 rather than the usual 1000, since these plots
#' are for rough diagnostic purposes.
#' @param interp If \code{interp="start"} (the default) then the entry time
#' into the absorbing state is assumed to be the time it is first observed in
#' the data.
#' 
#' If \code{interp="midpoint"} then the entry time into the absorbing state is
#' assumed to be halfway between the time it is first observed and the previous
#' observation time. This is generally more reasonable for "progressive" models
#' with observations at arbitrary times.
#' @param legend.pos Vector of the \eqn{x} and \eqn{y} position, respectively,
#' of the legend.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param lty Line type for the fitted curve. See \code{\link{par}}.
#' @param lwd Line width for the fitted curve. See \code{\link{par}}.
#' @param col Colour for the fitted curve. See \code{\link{par}}.
#' @param lty.ci Line type for the fitted curve confidence limits. See
#' \code{\link{par}}.
#' @param lwd.ci Line width for the fitted curve confidence limits. See
#' \code{\link{par}}.
#' @param col.ci Colour for the fitted curve confidence limits. See
#' \code{\link{par}}.
#' @param mark.time Mark the empirical survival curve at each censoring point,
#' see \code{\link[survival]{lines.survfit}}.
#' @param col.surv Colour for the empirical survival curve, passed to
#' \code{\link[survival]{lines.survfit}}. See \code{\link{par}}.
#' @param lty.surv Line type for the empirical survival curve, passed to
#' \code{\link[survival]{lines.survfit}}. See \code{\link{par}}.
#' @param lwd.surv Line width for the empirical survival curve, passed to
#' \code{\link[survival]{lines.survfit}}. See \code{\link{par}}.
#' @param survdata Set to \code{TRUE} to return the survival data frame
#' constructed when plotting the empirical curve.  This can be used for
#' constructing survival plots by hand using
#' \code{\link[survival]{plot.survfit}}.
#' @param ... Other arguments to be passed to the \code{\link{plot}} function
#' which draws the fitted curve, or the \code{\link[survival]{lines.survfit}}
#' function which draws the empirical curve.
#' @seealso \code{\link[survival]{survfit}},
#' \code{\link[survival]{plot.survfit}}, \code{\link{plot.prevalence.msm}}
#' @references Turnbull, B. W. (1976) The empirical distribution function with
#' arbitrarily grouped, censored and truncated data. J. R. Statist. Soc. B 38,
#' 290-295.
#' 
#' Fay, MP and Shaw, PA (2010). Exact and Asymptotic Weighted Logrank Tests for
#' Interval Censored Data: The interval R package. Journal of Statistical
#' Software. http://www.jstatsoft.org/v36/ i02/. 36 (2):1-34.
#' @keywords models
#' @export plot.survfit.msm
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



#' Calculate tables of hazard ratios for covariates on transition intensities
#' 
#' Hazard ratios are computed by exponentiating the estimated covariate effects
#' on the log-transition intensities.  This function is called by
#' \code{\link{summary.msm}}.
#' 
#' 
#' @param x Output from \code{\link{msm}} representing a fitted multi-state
#' model.
#' @param hazard.scale Vector with same elements as number of covariates on
#' transition rates. Corresponds to the increase in each covariate used to
#' calculate its hazard ratio. Defaults to all 1.
#' @param cl Width of the symmetric confidence interval to present.  Defaults
#' to 0.95.
#' @return
#' 
#' A list of tables containing hazard ratio estimates, one table for each
#' covariate.  Each table has three columns, containing the hazard ratio, and
#' an approximate upper and lower confidence limit respectively (assuming
#' normality on the log scale), for each Markov chain transition intensity.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}, \code{\link{summary.msm}},
#' \code{\link{odds.msm}}
#' @keywords models
#' @export hazard.msm
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



#' Calculate tables of odds ratios for covariates on misclassification
#' probabilities
#' 
#' Odds ratios are computed by exponentiating the estimated covariate effects
#' on the logit-misclassification probabilities.
#' 
#' 
#' @param x Output from \code{\link{msm}} representing a fitted multi-state
#' model.
#' @param odds.scale Vector with same elements as number of covariates on
#' misclassification probabilities. Corresponds to the increase in each
#' covariate used to calculate its odds ratio. Defaults to all 1.
#' @param cl Width of the symmetric confidence interval to present.  Defaults
#' to 0.95.
#' @return
#' 
#' A list of tables containing odds ratio estimates, one table for each
#' covariate.  Each table has three columns, containing the odds ratio, and an
#' approximate upper 95\% and lower 95\% confidence limit respectively
#' (assuming normality on the log scale), for each misclassification
#' probability.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}, \code{\link{hazard.msm}}
#' @keywords models
#' @export odds.msm
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
#' the model is not a hidden Markov model and there are no censored state
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



#' Score residuals
#' 
#' Score residuals for detecting outlying subjects.
#' 
#' The score residual for a single subject is
#' 
#' \deqn{U(\theta)^T I(\theta)^{-1} U(\theta)}{U(theta)^T I(theta)^{-1}
#' U(theta)}
#' 
#' where \eqn{U(\theta)}{U(theta)} is the vector of first derivatives of the
#' log-likelihood for that subject at maximum likelihood estimates
#' \eqn{\theta}{theta}, and \eqn{I(\theta)}{theta} is the observed Fisher
#' information matrix, that is, the matrix of second derivatives of minus the
#' log-likelihood for that subject at theta.
#' 
#' Subjects with a higher influence on the maximum likelihood estimates will
#' have higher score residuals.
#' 
#' These are only available for models with analytic derivatives (which
#' includes all non-hidden and most hidden Markov models).
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param plot If \code{TRUE}, display a simple plot of the residuals in
#' subject order, labelled by subject identifiers
#' @return Vector of the residuals, named by subject identifiers.
#' @author Andrew Titman \email{a.titman@@lancaster.ac.uk} (theory), Chris
#' Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk} (code)
#' @keywords models
#' @export scoreresid.msm
scoreresid.msm <- function(x, plot=FALSE){
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    if (!deriv_supported(x$data, x$hmodel, x$cmodel))
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



#' Expected first passage time
#' 
#' Expected time until first reaching a particular state or set of states in a
#' Markov model.
#' 
#' The expected first passage times from each of a set of states
#' \eqn{\mathbf{i}}{i} to to the remaining set of states
#' \eqn{\overline{\mathbf{i}}}{ibar} in the state space, for a model with
#' transition intensity matrix \eqn{Q}, are
#' 
#' \deqn{-Q_{\mathbf{i},\mathbf{i}}^{-1} \mathbf{1}}{-Q_{i,i}^{-1} 1}
#' 
#' where \eqn{\mathbf{1}}{1} is a vector of ones, and
#' \eqn{Q_{\mathbf{i},\mathbf{i}}}{Q_{i,i}} is the square subset of \eqn{Q}
#' pertaining to states \eqn{\mathbf{i}}{i}.
#' 
#' It is equal to the sum of mean sojourn times for all states between the
#' "from" and "to" states in a unidirectional model.  If there is non-zero
#' chance of reaching an absorbing state before reaching \code{tostate}, then
#' it is infinite.  It is trivially zero if the "from" state equals
#' \code{tostate}.
#' 
#' This function currently only handles time-homogeneous Markov models.  For
#' time-inhomogeneous models it will assume that \eqn{Q} equals the average
#' intensity matrix over all times and observed covariates.  Simulation might
#' be used to handle time dependence.
#' 
#' Note this is the \emph{expectation} of first passage time, and the
#' confidence intervals are CIs for this mean, not predictive intervals for the
#' first passage time.  The full distribution of the first passage time to a
#' set of states can be obtained by setting the rows of the intensity matrix
#' \eqn{Q} corresponding to that set of states to zero to make a model where
#' those states are absorbing.  The corresponding transition probability matrix
#' \eqn{Exp(Qt)} then gives the probabilities of having hit or passed that
#' state by a time \eqn{t} (see the example below). This is implemented in
#' \code{\link{ppass.msm}}.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param qmatrix Instead of \code{x}, you can simply supply a transition
#' intensity matrix in \code{qmatrix}.
#' @param tostate State, or set of states supplied as a vector, for which to
#' estimate the first passage time into.  Can be integer, or character matched
#' to the row names of the Q matrix.
#' @param start Starting state (integer).  By default (\code{start="all"}),
#' this will return a vector of expected passage times from each state in turn.
#' 
#' Alternatively, this can be used to obtain the expected first passage time
#' from a \emph{set} of states, rather than single states.  To achieve this,
#' \code{state} is set to a vector of weights, with length equal to the number
#' of states in the model.  These weights should be proportional to the
#' probability of starting in each of the states in the desired set, so that
#' weights of zero are supplied for other states.  The function will calculate
#' the weighted average of the expected passage times from each of the
#' corresponding states.
#' @param covariates Covariate values defining the intensity matrix for the
#' fitted model \code{x}, as supplied to \code{\link{qmatrix.msm}}.
#' @param ci If \code{"normal"}, then calculate a confidence interval by
#' simulating \code{B} random vectors from the asymptotic multivariate normal
#' distribution implied by the maximum likelihood estimates (and covariance
#' matrix) of the log transition intensities and covariate effects.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' 
#' If \code{"none"} (the default) then no confidence interval is calculated.
#' @param cl Width of the symmetric confidence interval, relative to 1.
#' @param B Number of bootstrap replicates.
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @param ... Arguments to pass to \code{\link{MatrixExp}}.
#' @return A vector of expected first passage times, or "hitting times", from
#' each state to the desired state.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{sojourn.msm}}, \code{\link{totlos.msm}},
#' \code{\link{boot.msm}}.
#' @references Norris, J. R. (1997) Markov Chains. Cambridge University Press.
#' @keywords models
#' @examples
#' 
#' twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
#'              c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
#' efpt.msm(qmatrix=twoway4.q, tostate=3)
#' # given in state 1, expected time to reaching state 3 is infinite
#' # since may die (state 4) before entering state 3
#' 
#' # If we remove the death state from the model, EFPTs become finite
#' Q <- twoway4.q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
#' efpt.msm(qmatrix=Q, tostate=3)
#' 
#' # Suppose we cannot die or regress while in state 2, can only go to state 3
#' Q <- twoway4.q; Q[2,4] <- Q[2,1] <- 0; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
#' efpt.msm(qmatrix=Q, tostate=3)
#' # The expected time from 2 to 3 now equals the mean sojourn time in 2.
#' -1/Q[2,2]
#' 
#' # Calculate cumulative distribution of the first passage time
#' # into state 3 for the following three-state model
#' Q <- twoway4.q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
#' # Firstly form a model where the desired hitting state is absorbing
#' Q[3,] <- 0
#' MatrixExp(Q, t=10)[,3]
#' ppass.msm(qmatrix=Q, tot=10)
#' # Given in state 1 at time 0, P(hit 3 by time 10) = 0.479
#' MatrixExp(Q, t=50)[,3]  # P(hit 3 by time 50) = 0.98
#' ppass.msm(qmatrix=Q, tot=50)
#' 
#' 
#' @export efpt.msm
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



#' Passage probabilities
#' 
#' Probabilities of having visited each state by a particular time in a
#' continuous time Markov model.
#' 
#' The passage probabilities to state \eqn{s} are computed by setting the
#' \eqn{s}th row of the transition intensity matrix \eqn{Q} to zero, giving an
#' intensity matrix \eqn{Q*} for a simplified model structure where state
#' \eqn{s} is absorbing.  The probabilities of passage are then equivalent to
#' row \eqn{s} of the transition probability matrix \eqn{Exp(tQ*)} under this
#' simplified model for \eqn{t=}\code{tot}.
#' 
#' Note this is different from the probability of occupying each state at
#' exactly time \eqn{t}, given by \code{\link{pmatrix.msm}}.  The passage
#' probability allows for the possibility of having visited the state before
#' \eqn{t}, but then occupying a different state at \eqn{t}.
#' 
#' The mean of the passage distribution is the expected first passage time,
#' \code{\link{efpt.msm}}.
#' 
#' This function currently only handles time-homogeneous Markov models.  For
#' time-inhomogeneous models the covariates are held constant at the value
#' supplied, by default the column means of the design matrix over all
#' observations.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param qmatrix Instead of \code{x}, you can simply supply a transition
#' intensity matrix in \code{qmatrix}.
#' @param tot Finite time to forecast the passage probabilites for.
#' @param start Starting state (integer).  By default (\code{start="all"}),
#' this will return a matrix one row for each starting state.
#' 
#' Alternatively, this can be used to obtain passage probabilities from a
#' \emph{set} of states, rather than single states.  To achieve this,
#' \code{state} is set to a vector of weights, with length equal to the number
#' of states in the model.  These weights should be proportional to the
#' probability of starting in each of the states in the desired set, so that
#' weights of zero are supplied for other states.  The function will calculate
#' the weighted average of the passage probabilities from each of the
#' corresponding states.
#' @param covariates Covariate values defining the intensity matrix for the
#' fitted model \code{x}, as supplied to \code{\link{qmatrix.msm}}.
#' @param piecewise.times Currently ignored: not implemented for
#' time-inhomogeneous models.
#' @param piecewise.covariates Currently ignored: not implemented for
#' time-inhomogeneous models.
#' @param ci If \code{"normal"}, then calculate a confidence interval by
#' simulating \code{B} random vectors from the asymptotic multivariate normal
#' distribution implied by the maximum likelihood estimates (and covariance
#' matrix) of the log transition intensities and covariate effects.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' 
#' If \code{"none"} (the default) then no confidence interval is calculated.
#' @param cl Width of the symmetric confidence interval, relative to 1.
#' @param B Number of bootstrap replicates.
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @param ... Arguments to pass to \code{\link{MatrixExp}}.
#' @return A matrix whose \eqn{r, s} entry is the probability of having visited
#' state \eqn{s} at least once before time \eqn{t}, given the state at time
#' \eqn{0} is \eqn{r}.  The diagonal entries should all be 1.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{efpt.msm}}, \code{\link{totlos.msm}},
#' \code{\link{boot.msm}}.
#' @references Norris, J. R. (1997) Markov Chains. Cambridge University Press.
#' @keywords models
#' @examples
#' 
#' Q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
#'            c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
#' 
#' ## ppass[1,2](t) converges to 0.5 with t, since given in state 1, the
#' ## probability of going to the absorbing state 4 before visiting state
#' ## 2 is 0.5, and the chance of still being in state 1 at t decreases.
#' 
#' ppass.msm(qmatrix=Q, tot=2)
#' ppass.msm(qmatrix=Q, tot=20)
#' ppass.msm(qmatrix=Q, tot=100)
#' 
#' Q <- Q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
#' 
#' ## Probability of about 1/2 of visiting state 3 by time 10.5, the
#' ## median first passage time
#' 
#' ppass.msm(qmatrix=Q, tot=10.5)
#' 
#' ## Mean first passage time from state 2 to state 3 is 10.02: similar
#' ## to the median
#' 
#' efpt.msm(qmatrix=Q, tostate=3)
#' 
#' @export ppass.msm
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
