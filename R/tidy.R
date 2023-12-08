##' Tidy the parameter estimates from an msm model
##'
##' @param x Object returned by \code{\link{msm}}, representing a fitted
##'   multi-state model.
##'
##' @param ... Other arguments (currently unused).
##'
##' @return A "tibble", with one row for each parameter and the following
##'   columns describing the parameter.
##'
##' * `parclass`: Class of parameters: `intens` (transition intensities), `hr`
##'   (hazard ratios representing effects of covariates on intensities), and
##'   their transformed versions `logintens` (log intensities) and `loghr` (log
##'   hazard ratios).
##'
##'   For "misclassification" models fitted with the `ematrix` argument to `msm`,
##'   other classes of parameters include `misc` (misclassification
##'   probabilities), `logitmisc` (misclassification log odds), `or_misc` and
##'   `logor_misc` (effects of covariates on misclassification probabilities, as
##'   odds ratios or log odds ratios, with the first state as the reference
##'   category).
##'
##'   For hidden Markov models fitted with the `hmodel` argument to `msm`, the
##'   parameter class called `hmm` comprises the parameters of the distributions
##'   of the outcome conditionally on the hidden state.  Covariates on the
##'   location parameter of these distributions are included in class `hmmcov`.
##'   If initial state occupancy probabilities are estimated, these are included
##'   in class `initp` (or `initlogodds` for the log odds transforms of these),
##'   and any covariates on these probabilities are included in class `initpcov`.
##'
##' * `state`: Starting state of the transition for transition intensities, and
##'    true state for misclassification probabilities or hidden Markov model parameters.
##'
##' * `tostate`: Ending state of the transition for transition intensities, and 
##'    observed state for misclassification probabilities 
##'
##' * `term`: Name of the covariate for covariate effects, or "baseline" for the 
##'   baseline intensity or analogous parameter value. 
##'   Note that the "baseline" parameters are the parameters with covariates
##'   set to their mean values in the data (stored in e.g. `x$qcmodel$covmeans`),
##'   unless `msm` was called with `center=FALSE`.
##'
##' * `estimate`, `std.error`, `conf.low`, `conf.high`:  Parameter estimate, 
##'   standard error, and lower and upper confidence limits. 
##'
##' * `statistic`, `p.value`: For covariate effects, the Z-test statistic and p-value
##'   for a test of the null hypothesis that the covariate effect is zero, based
##'   on the estimate and standard error. 
##' 
##' @md
##' @importFrom generics tidy
##' @importFrom tibble tibble
##'
##' @export
tidy.msm <- function(x, ...){
  tidynames <- c("parclass","state","tostate","term","estimate")
  tidycinames <- c("std.error","conf.low","conf.high")
  
  qkeep <- which(x$qmodel$imatrix==1, arr.ind=TRUE)
  
  res <- tidy_mats(x$Qmatrices, qkeep)
  res$parclass <- ifelse(res$term=="baseline", "intens",
                         ifelse(res$term=="logbaseline",
                                "logintens", "loghr"))
  ## workaround if any covariates named "baseline"
  res$term <- map_covnames(res$term, from=attr(x$Qmatrices,"covlabels"), to=attr(x$Qmatrices,"covlabels.orig")) 
  if (x$foundse){
    xnames <- c("QmatricesSE", "QmatricesL", "QmatricesU")
    for (i in seq_along(xnames)){
      res[[tidycinames[i]]] <- tidy_mats(x[[xnames[i]]], qkeep)$estimate
    }
    tidynames <- c(tidynames, tidycinames)
  }
  ## convert log hazard ratios to hazard ratios
  coveffs <- res[res$parclass == "loghr",]
  if (nrow(coveffs) > 0){
    coveffs$parclass <- "hr"
    coveffs$std.error <- exp(coveffs$estimate) * coveffs$std.error
    coveffs$estimate <- exp(coveffs$estimate)
    if (x$foundse){
      for (i in c("conf.low","conf.high"))
        coveffs[[i]] <- exp(coveffs[[i]])
    }
    res <- rbind(res, coveffs)
  }
  
  if (x$emodel$misc){
    qkeep <- which(x$emodel$imatrix==1, arr.ind=TRUE)
    rese <- tidy_mats(x$Ematrices, qkeep)
    rese$parclass <- ifelse(rese$term=="baseline", "misc",
                           ifelse(rese$term=="logbaseline",
                                  "logitmisc", "logor_misc"))
    rese$term <- map_covnames(rese$term, from=attr(x$Ematrices,"covlabels"), to=attr(x$Ematrices,"covlabels.orig")) 
    if (x$foundse){
      xnames <- c("EmatricesSE", "EmatricesL", "EmatricesU")
      for (i in seq_along(xnames)){
        rese[[tidycinames[i]]] <- tidy_mats(x[[xnames[i]]], qkeep)$estimate
      }
    }
    ## convert log odds ratios to odds ratios
    coveffs <- rese[rese$parclass == "logor_misc",]
    if (nrow(coveffs) > 0){
      coveffs$parclass <- "or_misc"
      coveffs$estimate <- exp(coveffs$estimate)
      if (x$foundse){
        coveffs$std.error <- coveffs$estimate * coveffs$std.error
        for (i in c("conf.low","conf.high"))
          coveffs[[i]] <- exp(coveffs[[i]])
      }
      rese <- rbind(rese, coveffs)
    }
    res <- rbind(res, rese)
  }

  res$term[res$term=="logbaseline"] <- "baseline"

  if (x$hmodel$hidden && !x$emodel$misc) { 
    resh <- tidy.hmodel(x)
    resh$tostate <- NA
    res <- rbind(res, resh[,colnames(res)])
  }
  res <- tibble::tibble(res)[tidynames]
  
  ## test statistics and p-values
  covs <- which(res$parclass %in% c("loghr","logor_misc","hmmcov","initpcov"))
  if (length(covs) > 0 && x$foundse){
    res$statistic <- res$p.value <- NA
    res$statistic[covs] <- res$estimate[covs] / res$std.error[covs]
    res$p.value[covs] <- 2 * pnorm(-abs(res$statistic[covs]))
    res$statistic[res$parclass %in% c("hr","or_misc")] <-
      res$statistic[res$parclass %in% c("loghr","logor_misc")]
    res$p.value[res$parclass %in% c("hr","or_misc")] <-
      res$p.value[res$parclass %in% c("loghr","logor_misc")]
  }

  res
# perhaps this could be an argument  
#  statenames <- rownames(x$qmodel$imatrix)
#  if (!is.null(statenames)){
#    res$fromname <- statenames[res$from]
#    res$toname <- statenames[res$to]
# res$state <- statenames[res$state]
#  }
  
}

map_covnames <- function(x, from, to){
  if (!is.null(from)){
    for (i in seq_along(from)){
      x[x==from[i]] <- to[i]
    }
  }
  x
}

## Tidier for a list of matrices with one component for the baseline intensity
## matrix and further components for covariate effects on this 

tidy_mats <- function(x, qkeep=NULL){
  if (is.null(qkeep)){
    qkeep <- which(x$estimates > 0, arr.ind=TRUE)
  }
  colnames(qkeep) <- c("state","tostate")
  resq <- lapply(x, function(y){
    cbind(data.frame(res=unclass(y)[qkeep]), qkeep)
  })
  for (i in seq_along(resq)){ 
    resq[[i]]$term <- names(resq)[i]
  }
  resq <- do.call("rbind", resq)
  rownames(resq) <- NULL
  names(resq)[names(resq)=="res"] <- "estimate"
  resq
}

##' Tidy the output of pmatrix.msm and similar functions
##' 
##' This is the method for the generic `tidy` function that is 
##' used for tidying the output of \code{\link{qmatrix.msm}}, \code{\link{pmatrix.msm}}, 
##' \code{\link{ematrix.msm}}, \code{\link{pnext.msm}} or \code{\link{ppass.msm}}.
##'  This should be called as 
##' \code{tidy()}, not \code{tidy.msm.est()} or \code{tidy.qmatrix()} or anything else.
##' 
##' @param x Output of \code{\link{qmatrix.msm}}, \code{\link{pmatrix.msm}}, 
##' \code{\link{ematrix.msm}}, \code{\link{pnext.msm}} or \code{\link{ppass.msm}},
##'  which all return objects of class \code{"msm.est"}.
##' 
##' @param ... Further arguments (unused).
##' 
##' @export
tidy.msm.est <- function(x, ...){
  if (is.matrix(x)) x <- list(estimates=x) # no CIs available
  tm <- tidy_mats(x)
  oldnames <- c("estimates","SE","L","U")
  tidynames <- c("estimate","std.error","conf.low","conf.high")
  for (i in seq_along(oldnames)){
    tm$term[tm$term==oldnames[i]] <- tidynames[i]
  }
  res <- reshape(tm, direction="wide", idvar=c("state","tostate"), timevar="term")
  names(res) <- gsub("estimate\\.","",names(res))
  statenames <- rownames(x$estimates)
  if (!is.null(statenames)){
    res$statename <- statenames[res$state]
    res$tostatename <- statenames[res$tostate]
  }
  tibble::tibble(res)
}

## Tidier for a "hmodel" object, which is one of the components of a "msm" object
## for a hidden Markov model. 

tidy.hmodel <- function(x){
  xh <- x$hmodel
  p <- x$paramdata
  res <- data.frame(parclass = "hmm",
                    state = xh$parstate,
                    term = xh$plabs,
                    estimate = xh$pars)
  if (x$foundse){
    hbasepars <- which(!p$plabs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov"))
    hse <- sqrt(diag(x$covmat[hbasepars,hbasepars])) 
    np <- length(hbasepars)
    hse <- dgexpit(p$params[hbasepars],  # delta method
                   xh$ranges[1:np,"lower"], xh$ranges[1:np,"upper"]) * hse
    res$std.error <- hse
    cis <- setNames(as.data.frame(xh$ci), c("conf.low","conf.high"))
    res <- cbind(res, cis)
  }
  
  hcovpars <- which(p$plabs == "hcov")
  if (length(hcovpars) > 0){
    ce <- data.frame(parclass = "hmmcov", 
                     state = xh$coveffstate,
                     term = xh$covlabels,
                     estimate = xh$coveffect)
    if (x$foundse){
      hcovse <- sqrt(diag(x$covmat[hcovpars,hcovpars]))
      ce$std.error <- hcovse
      cis <- setNames(as.data.frame(xh$covci), c("conf.low","conf.high"))
      ce <- cbind(ce, cis)
    }
    res <- rbind(res, ce)
  }
  
  if (xh$est.initprobs){
    est <- if (x$foundse) xh$initprobs[,"Estimate"] else xh$initprobs
    initp <- data.frame(parclass = "initp",
                        state = 1:x$hmodel$nstates,
                        term = "baseline",
                        estimate = est)
    iloest <- p$params[p$plabs=="initp"]
    ilose <- sqrt(diag(p$covmat))[p$plabs=="initp"]
    initlo <- data.frame(parclass = "initlogodds",
                         state = 2:x$hmodel$nstates,
                         term = "baseline",
                         estimate = iloest)
    if (x$foundse){
      initp <- cbind(initp, data.frame(
                        std.error = NA,
                        conf.low = xh$initprobs[,"LCL"],
                        conf.high = xh$initprobs[,"UCL"]))
      initlo <- cbind(initlo, data.frame(
        std.error = ilose,
        conf.low = iloest - qnorm(0.975) * ilose,
        conf.high = iloest + qnorm(0.975) * ilose))
    }
    res <- rbind(res, initp, initlo)

    if (xh$nicoveffs > 0){
      icest <- if (x$foundse) xh$icoveffect[,"Estimate"] else xh$icoveffect
      icovlabels <- names(icest)
      ic <- data.frame(parclass = "initpcov",
                       state = 2:xh$nstates,
                       term = icovlabels,
                       estimate = icest)
      if (x$foundse){
        cis <- data.frame(
          std.error = (xh$icoveffect[,"UCL"] - xh$icoveffect[,"LCL"])/(2*qnorm(0.975)),
          conf.low = xh$icoveffect[,"LCL"],
          conf.high = xh$icoveffect[,"UCL"])
        ic <- cbind(ic, cis)
      }
      res <- rbind(res, ic)
    }
  }

  if (x$foundse)
    res$std.error[is.na(res$conf.low)] <- NA
    
  tibble::tibble(res)

}

#' Tidy the output of prevalence.msm
#' 
#' Note this should be called as \code{tidy()} not \code{tidy.msm.prevalence()} or anything else, as this is 
#' a method for the generic \code{tidy()} function.
#' 
#' @param x Output of \code{\link{prevalence.msm}}.
#' 
#' @return A tibble with one row per combination of output type (count or percentage)
#' and state, and columns for observed value, expected value and confidence 
#' limits for the expected value (if available).
#' 
#' @param ... Further arguments (unused).
#' 
#' @export
tidy.msm.prevalence <- function(x, ...){
  if (is.list(x$Expected))
    tidy_msm_prevalence_ci(x,...)
  else
    tidy_msm_prevalence_noci(x,...)
}

tidy_msm_prevalence_noci <- function(x,...){
  for (i in c("Observed percentages","Expected percentages")){
    x[[i]] <- cbind(x[[i]], "Total" = 100)
  }
  for (i in seq_along(x)) { 
    x[[i]] <- cbind(output = names(x[i]), 
                    time = as.numeric(rownames(x[[i]])),
                    as.data.frame(x[[i]]))
  }
  x <- do.call(rbind, x)
  x$Total <- NULL
  x <- reshape(x, direction="long", varying=3:ncol(x), 
               v.names="prevalence", timevar="state", idvar=c("output","time"))
  tx <- tibble::tibble(x)
  txobs <- tx[tx$output %in% c("Observed","Observed percentages"),]
  txexp <- tx[tx$output %in% c("Expected","Expected percentages"),]
  tx <- txobs
  names(tx)[names(tx)=="prevalence"] <- "observed"
  tx$expected <- txexp$prevalence
  tx$output[tx$output=="Observed"] <- "count"
  tx$output[tx$output=="Observed percentages"] <- "percentage"
  tx
}

tidy_msm_prevalence_ci <- function(x,...){
  nt <- dim(x$Expected$ci)[2]
  nst <- ncol(x$Expected$estimates)-1
  ecis <- rbind(cbind(summary="observed", 
                      unname(as.data.frame(x$Observed[,-nt]))),
                cbind(summary="expected", 
                      unname(as.data.frame(x$Expected$estimates[,-nt]))),
                cbind(summary="conf.low", 
                      unname(as.data.frame(x$Expected$ci[,-nt,1]))),
                cbind(summary="conf.high", 
                      unname(as.data.frame(x$Expected$ci[,-nt,2]))))
  ecis$output <- "count"
  epcis <- rbind(cbind(summary="observed", 
                       unname(as.data.frame(x$Observed[,-nt]))),
                 cbind(summary="expected", 
                       unname(as.data.frame(x$`Expected percentages`$estimates[,-nt]))),
                 cbind(summary="conf.lower", 
                       unname(as.data.frame(x$`Expected percentages`$ci[,,1]))),
                 cbind(summary="conf.high", 
                       unname(as.data.frame(x$`Expected percentages`$ci[,,2]))))
  epcis$output <- "percentage"
  ecis <- rbind(ecis, epcis)
  ecis$time <- as.numeric(rownames(x$Expected$estimates))
  cilong <- reshape(ecis, direction="long", varying=1 + 1:nst, 
                    v.names="est", 
                    timevar="state", 
                    idvar=c("output","time","summary"))
  ciwide <- reshape(cilong, direction="wide", idvar=c("output","time","state"), 
                    timevar="summary", v.names="est")
  rownames(ciwide) <- NULL
  names(ciwide)[4:7] <- gsub("est.","", names(ciwide)[4:7])
  tibble::tibble(ciwide)
}


#' Tidy the output of totlos.msm and similar functions
#' 
#' Note this should be called as \code{tidy()} not \code{tidy.msm.totlos()} or anything else, as this is 
#' a method for the generic \code{tidy()} function.
#' 
#' @param x Output of \code{\link{totlos.msm}}, \code{\link{envisits.msm}}
#' or \code{\link{efpt.msm}}, which return objects of class \code{"msm.estbystate"}.
#' 
#' @return A tibble with one row per state, and columns for the estimate, and
#'  confidence intervals if available.
#' 
#' @param ... Further arguments (unused).
#' 
#' @export
tidy.msm.estbystate <- function(x, ...){
  if (!is.matrix(x)) x <- matrix(x, nrow=1)
  x <- as.data.frame(t(x))
  names(x) <- if(ncol(x)==1) "estimate" else c("estimate","conf.low","conf.high")
  x <- cbind(state = 1:nrow(x), x)
  tibble::tibble(x)
}
