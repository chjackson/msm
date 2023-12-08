

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

mattotrans <- function(x, matrix, lower=NULL, upper=NULL, fixed=NULL, keep.diag=FALSE, intmisc="intens"){
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
    if (!is.null(fixed)) res[,4] <- fixed[keep]
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
    covlabels <-  make.unique(c("Baseline", x$qcmodel$covlabels))[-(1)]
    colnames(fres) <- c("Baseline", covlabels)
    rownames(fres) <- rownames(y)
    fres[,1] <- format.ci(y[,1],y[,2],y[,3],y[,4],digits=digits,...)
    im <- t(x$qmodel$imatrix); diag(im) <- -colSums(im); nd <- which(im[im!=0]==1)
    for (i in seq(length.out=x$qcmodel$ncovs)){
      nm <- x$qcmodel$covlabels[[i]]
      nmout <- covlabels[[i]]
      ## FIXME if Baseline is a cov, then Qmatrices are unaffected
      ## but output should be changed. ret sohuld be named by covlabels
      
      hrs <- mattotrans(x, x$Qmatrices[[nm]], x$QmatricesL[[nm]], x$QmatricesU[[nm]],
                          x$QmatricesFixed[[nm]], keep.diag=FALSE)
      hrs[,1:3] <- exp(hrs[,1:3])
      ret[nmout] <- matrix(ncol=3, nrow=nrow(ret), dimnames=list(NULL,colnames(hrs)[1:3]))
      ret[nd,nmout] <- hrs[,1:3,drop=FALSE]
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


