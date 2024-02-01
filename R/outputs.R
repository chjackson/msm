### Functions related to extracting intensity and probability matrices from msm models

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
    if (ci=="none"){
      if (sojourn) res <- soj
      else {
        res <- mat
        class(res) <- "msm.est"
      }
    }
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
    }
    class(res) <- "msm.est"
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

#' @export
"[.msm.est" <- function(x, i, j, drop=FALSE){
    Narg <- nargs() - (!missing(drop)) # number of args including x, excluding drop
    if ((missing(i) && missing(j)))
        res <- x
    else if (!is.list(x)){ 
      res <- unclass(x)
      if (missing(j) && (Narg==2))
        res <- res[i]
      else
        res <- res[i,j]
    }
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

#' @noRd
format.ci <- function(x, l, u, noci=NULL, digits=NULL, ...)
{
    if (is.null(noci)) noci <- rep(FALSE, length(x))
    noci[is.na(noci)] <- FALSE
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
#'
#' @param t1 The start of the time interval to estimate the transition
#' probabilities for.
#'
#' @param t2 The end of the time interval to estimate the transition
#' probabilities for.
#'
#' @param times Cut points at which the transition intensity matrix changes.
#' 
#' @param covariates A list with number of components one greater than the
#' length of \code{times}.  Each component of the list is specified in the same
#' way as the \code{covariates} argument to \code{\link{pmatrix.msm}}.  The
#' components correspond to the covariate values in the intervals
#' 
#' \code{(t1, times[1]], (times[1], times[2]], ..., (times[length(times)], t2]}
#' 
#' (assuming that all elements of \code{times} are in the interval \code{(t1,
#' t2)}).
#'
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
#'
#' @param cl Width of the symmetric confidence interval, relative to 1.
#'
#' @param B Number of bootstrap replicates, or number of normal simulations
#' from the distribution of the MLEs
#'
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#'
#' @param qlist A list of transition intensity matrices, of length one greater
#' than the length of \code{times}.  Either this or a fitted model \code{x}
#' must be supplied.  No confidence intervals are available if (just)
#' \code{qlist} is supplied.
#'
#' @param ... Optional arguments to be passed to \code{\link{MatrixExp}} to
#' control the method of computing the matrix exponential.
#'
#' @return The matrix of estimated transition probabilities \eqn{P(t)} for the
#' time interval \code{[t1, tn]}.  That is, the probabilities of occupying
#' state \eqn{s} at time \eqn{t_n}{tn} conditionally on occupying state \eqn{r}
#' at time \eqn{t_1}{t1}.  Rows correspond to "from-state" and columns to
#' "to-state".
#'
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#'
#' @seealso \code{\link{pmatrix.msm}}
#'
#' @keywords models
#'
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
                                  covariates=NULL, # list of lists of covariates, for (, times1], (times1, times2], ...
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
    if (is.null(qlist)) validate_piecewise(times, covariates)
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

validate_piecewise <- function(times, covariates){
  if (!is.numeric(times) || is.unsorted(times)) 
    stop("times should be a vector of numbers in increasing order")
  if ((length(covariates) != length(times) + 1))
    stop("Number of covariate lists must be one greater than the number of cut points")
}


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
    Q <- unclass(qmatrix.msm(x, covariates, ci="none"))
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
#' @export coef.msm
#' @export
coef.msm <- function(object, ...)
{
    if (!inherits(object, "msm")) stop("expected object to be a msm model")
    if (object$emodel$misc)
        object[c("Qmatrices", "Ematrices")]
    else object$Qmatrices
}



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
#' @export logLik.msm
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
        covlabels <- attr(x$Qmatrices, "covlabels") # rewritten as baseline.1 if any covariates called "baseline"
        if (x$foundse) {
            for (i in 1:x$qcmodel$ncovs) {
                cov <- covlabels[i]
                haz.rat <- t(exp(hazard.scale[i]*x$Qmatrices[[cov]]))[keepvec]
                LCL <- t(exp(hazard.scale[i]*(x$Qmatrices[[cov]] - qnorm(1 - 0.5*(1 - cl))*x$QmatricesSE[[cov]]) ))[keepvec]
                UCL <- t(exp(hazard.scale[i]*(x$Qmatrices[[cov]] + qnorm(1 - 0.5*(1 - cl))*x$QmatricesSE[[cov]]) ))[keepvec]
                haz.tab <- cbind(haz.rat, LCL, UCL)
                dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs),
                                          c("HR", "L", "U"))
                haz.list[[x$qcmodel$covlabels[i]]] <- haz.tab
            }
        }
        else {
            for (i in 1:x$qcmodel$ncovs) {
                cov <- x$qcmodel$covlabels[i]
                haz.tab <- as.matrix(t(exp(hazard.scale[i]*x$Qmatrices[[cov]]))[keepvec])
                dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs), "HR")
                haz.list[[x$qcmodel$covlabels[i]]] <- haz.tab
            }
        }
    }
    else haz.list <- "No covariates on transition intensities"
    haz.list
}



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
