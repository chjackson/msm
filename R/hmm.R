### FUNCTIONS FOR HIDDEN MARKOV MODELS IN CONTINUOUS TIME
### WITH ARBITRARY RESPONSE DISTRIBUTION

print.hmmMVdist <- function(x, ...)
{
    cat(sprintf("Multivariate hidden Markov model with %d outcomes:\n", length(x)))
    for (i in x) print(i)
}

print.hmmdist <- function(x, ...)
{
    cat("Hidden Markov model", x$label, "distribution\n\n")
    pnames <- if(x$label=="categorical") paste("P(",seq(x$pars[1]),")",sep="") else names(x$pars)
    pars <- if(x$label=="categorical") x$pars[3:(2+x$pars[1])] else x$pars
    cat("Parameters: ", paste(paste(pnames, pars, sep=" = "), collapse=", "))
    cat("\n")
}

msm.check.hmodel <- function(hmodel, nstates)
  {
      if (is.null(hmodel)) stop("Hidden model not specified")
      if (!is.list(hmodel)) stop("Hidden model should be a list")
      if (length(hmodel) != nstates) stop("hmodel of length ", length(hmodel), ", expected ", nstates)
      for (i in hmodel) {
          if (!inherits(i, "hmmdist")) stop("hmodel should be a list of HMM distribution objects")
      }
  }

msm.check.hcovariates <- function(hcovariates, qmodel)
  {
      if (!is.list(hcovariates)) stop("hcovariates should be a list")
      if (length(hcovariates) != qmodel$nstates)
        stop("hcovariates of length ", length(hcovariates), ", expected ", qmodel$nstates)
      for (i in hcovariates) {
          if (!is.null(i) )
            if ( class(i) != "formula") stop("hcovariates should be a list of formulae or NULLs")
      }
  }

msm.form.hmodel <- function(hmodel, hconstraint=NULL, initprobs=NULL, est.initprobs)
  {
      nst <- length(hmodel)
      if (is.null(initprobs))
          initprobs <- if (est.initprobs) rep(1/nst, nst) else c(1, rep(0, nst-1))
      else {
          if (!is.numeric(initprobs)) stop("initprobs should be numeric")
          if (is.matrix(initprobs)) {
              if (ncol(initprobs) != nst) stop("initprobs matrix has ", ncol(initprobs), " columns, should be number of states = ", nst)
              if (est.initprobs) { warning("Not estimating initial state occupancy probabilities since supplied as a matrix") }
              initprobs <- initprobs / rowSums(initprobs)
              est.initprobs <- FALSE
          }
          else  {
              if (length(initprobs) != nst) stop("initprobs of length ", length(initprobs), ", should be ", nst)
              initprobs <- initprobs / sum(initprobs)
              if (est.initprobs && any(initprobs==1)) {
                  est.initprobs <- FALSE
                  warning("Not estimating initial state occupancy probabilities, since some are fixed to 1")
              }
          }
      }
      nipars <- if (est.initprobs) nst else 0
      if (any(sapply(hmodel, inherits, "hmmMVdist"))){
          hmod <- msm.form.mvhmodel(hmodel)
      } else hmod <- msm.form.univhmodel(hmodel)
      hmod <- c(list(hidden=TRUE, nstates=nst, fitted=FALSE, nipars=nipars,
                     initprobs=initprobs, est.initprobs=est.initprobs, ematrix=FALSE),
                hmod)
      class(hmod) <- "hmodel"
      hmod
  }

msm.form.univhmodel <- function(hmodel){
    nst <- length(hmodel)    
    labels <- sapply(hmodel, function(x) x$label)
    models <- match(labels, .msm.HMODELS)
    pars <- lapply(hmodel, function(x) x$pars)
    plabs <- lapply(hmodel, function(x) names(x$pars))
    ## where non-misclassified outcome for hmmIdent distribution is not specified, this is
    ## just the state
    pars[labels=="identity"][sapply(pars[labels=="identity"], length) == 0] <- which(labels=="identity")
    plabs[labels=="identity"] <- "which"
    names(plabs) <- paste("state", 1:nst, sep=".")
    npars <- sapply(pars, length)
    parstate <- rep(1:nst, npars)
    firstpar <- c(0, cumsum(npars)[-nst])
    pars <- as.numeric(unlist(pars))
    plabs <- unlist(plabs)
    locpars <- which(plabs == rep(.msm.LOCPARS[labels], npars))
    names(pars) <- plabs
    list(models=models, labels=labels, npars=npars, nout=rep(1, nst), mv=FALSE,
         totpars=sum(npars), pars=pars, plabs=plabs, parstate=parstate,
         firstpar=firstpar, locpars=locpars)
}

### CHANGES FROM UNIVARIATE TO MULTIVARIATE OUTCOME HMMS
### models should be matrix instead of vector, new dim given by nout   (same for labels)
### npars was vector with one for each state.  now matrix
### pars was ragged array, mapped to states with parstate.  should still be.  (same for plabs)
### need new vector parout to map pars to outcomes. 
### firstpar, used in C to point to first parameters for a univariate model.  was one for each state.  now needs to be a matrix

msm.form.mvhmodel <- function(hmodel){
    nst <- length(hmodel)    
    for (i in seq_along(hmodel))
        if (!inherits(hmodel[[i]], "hmmMVdist")) hmodel[[i]] <- list(hmodel[[i]])

    nout <- sapply(hmodel, length)
    models <- labels <- firstpar <- matrix(nrow=max(nout), ncol=nst)
    npars <- matrix(0, nrow=max(nout), ncol=nst)
    pars <- plabs <- parstate <- parout <- vector(nst, mode="list")
    for (i in 1:nst){
        labels[1:nout[i],i] <- sapply(hmodel[[i]], function(x) x$label)
        models[1:nout[i],i] <- match(labels[1:nout[i],i], .msm.HMODELS)
        pars[[i]] <- lapply(hmodel[[i]], function(x)x$pars)
        npars[1:nout[i],i] <- sapply(pars[[i]], length)
        pars[[i]][labels[,i]=="identity" & npars[,i]==0] <- i # TESTME - hmmIdent with no arg: par is the state
        plabs[[i]][labels[,i]=="identity" & npars[,i]==0] <- "which"
        plabs[[i]] <- lapply(hmodel[[i]], function(x)names(x$pars))
        parstate[[i]] <- rep(i, sum(npars[,i]))
        parout[[i]] <- rep(1:nout[i], npars[1:nout[i],i])
    }
    firstpar <- matrix(c(0, cumsum(npars)[-length(npars)]), nrow=max(nout), ncol=nst)
    firstpar[npars==0] <- NA
    pars <- unlist(pars); plabs <- unlist(plabs); parout <- unlist(parout); parstate <- unlist(parstate)
    locpars <- which(plabs == rep(.msm.LOCPARS[labels], npars))
    names(pars) <- plabs
    list(models=models, labels=labels, npars=npars, nout=nout, mv=TRUE,
         totpars=sum(npars), pars=pars, plabs=plabs, parstate=parstate, parout=parout,
         firstpar=firstpar, locpars=locpars)  
}

## NOTE removed whichcovh, whichcovh.orig

msm.form.hcmodel <- function(hmodel, mm, hcovinits, hconstraint)
  {
      nst <- hmodel$nstates
      ncovs <- if (is.null(mm)) rep(0, nst) else sapply(mm, function(x) {ncol(x)-1})
      ncovs2 <- rep(rep(0, nst), hmodel$npars)
      ncovs2[hmodel$locpars] <- ncovs[hmodel$parstate[hmodel$locpars]]
      coveffstate <- rep(1:nst, tapply(ncovs2, hmodel$parstate, sum))
      if (is.null(hcovinits)){
          coveffect <- rep(0, sum(ncovs2))
      }
      else {
          if (!(sum(ncovs2) == length(unlist(hcovinits)))) {
              warning("Initial values for hidden covariate effects do not match numbers of covariates, ignoring")
              coveffect <- rep(0, sum(ncovs2))
          } else coveffect <- unlist(hcovinits)
          if (!is.numeric(coveffect)) {
              warning("hcovinits should be numeric")
              coveffect <- rep(0, sum(ncovs2))
          }
      }
      covlabels <- lapply(mm, function(x) colnames(x)[-1])
      covlabels <- unlist(covlabels[hmodel$parstate[hmodel$locpars]])
      names(coveffect) <- covlabels
      hcmod <- list(ncovs=ncovs2, coveffect=coveffect, covlabels=covlabels,
                    coveffstate=coveffstate, ncoveffs=length(coveffect))
      hmodel <- c(hmodel, hcmod)
      hmodel$plabs[hmodel$plabs=="hcov"] <- paste("hcov.",covlabels,sep="")
      class(hmodel) <- "hmodel"
      hmodel
  }

### NOTE whichcovi removed

msm.form.icmodel <- function(hmodel, mm, icovinits) {
  nst <- hmodel$nstates
  nicovs <- ncol(mm) - 1
  nicovs <- rep(nicovs, nst-1)
  if (!is.matrix(hmodel$initprobs))
      nicovs[hmodel$initprobs[-1] == 0] <- 0 # don't estimate cov effects on probs which are fixed to zero
  if (is.null(icovinits))
    icoveffect <- rep(0, sum(nicovs))
  else {
      icoveffect <- unlist(icovinits)
      if (!(length(icoveffect) == sum(nicovs))) {
          warning("Initial values for initial state covariate effects do not match numbers of covariates, ignoring")
          icoveffect <- rep(0, sum(nicovs))
      }
      else if (!is.numeric(icoveffect)) {
          warning("icovinits should be numeric")
          icoveffect <- rep(0, sum(nicovs))
      }
  }
  names(icoveffect) <- rep(colnames(mm)[-1], each=sum(nicovs>0))
  icmod <- list(nicovs=nicovs, icoveffect=icoveffect, nicoveffs=length(icoveffect))
  hmodel <- c(hmodel, icmod)
  class(hmodel) <- "hmodel"
  hmodel
}

## Convert old-style misclassification model specification to a new-style HMM with categorical response

msm.emodel2hmodel <- function(emodel, qmodel)
  {
      nst <- qmodel$nstates
      if (emodel$misc) {
          hidden <- TRUE
          nepars <- rowSums(emodel$imatrix)
          models <- ifelse(apply(emodel$ematrix,1,function(x)any(x==1)), 2, 1)
          npars <- ifelse(models==1, 2 + nst, 1)
          pars <- plabs <- vector(nst, mode="list")
          parstate <- rep(1:nst, npars)
          names(pars) <- names(plabs) <- paste("state", 1:nst, sep=".")
          for (i in seq(nst)) {
              if (models[i]==1) {
                  ppars <- emodel$ematrix[i,]
                  plab <- rep("p", nst)
                  plab[ppars==0] <- "p0"
                  plab[i] <- "pbase"
                  plabs[[i]] <- c("ncats", "basecat", plab)
                  ## Baseline category is the probability of no misclassification (diagonal of ematrix)
                  pars[[i]] <- c(nst, i, ppars)
              }
              else {
                  plabs[[i]] <- "which"
                  pars[[i]] <- which(emodel$ematrix[i,]==1)
              }
          }
          firstpar <- c(0, cumsum(npars)[-qmodel$nstates])
          pars <- unlist(pars)
          plabs <- unlist(plabs)
          names(pars) <- plabs
          labels <- .msm.HMODELS[models]
          locpars <- which(plabs == rep(.msm.LOCPARS[labels], npars))
          hmod <- list(hidden=TRUE, fitted=FALSE, nstates=nst, models=models, labels=labels,
                       ematrix=TRUE, ## remember if obtained from ematrix, since could change meaning of obstrue
                       nout=rep(1, nst), mv=FALSE,
                       npars=npars, totpars=sum(npars), locpars=locpars,
                       pars=pars, plabs=plabs, parstate=parstate, firstpar=firstpar, nipars=emodel$nipars, initprobs=emodel$initprobs, est.initprobs=emodel$est.initprobs)
          hmod$constr <- msm.econstr2hconstr(emodel$constr, hmod)
      }
      else {
          hmod <- list(hidden=FALSE, fitted=FALSE, ematrix=FALSE, models=rep(0, qmodel$nstates), npars=0, ndpars=0)
      }
      class(hmod) <- "hmodel"
      hmod
  }

msm.misccov2hcov <- function(misccovariates, emodel)
  {
      nst <- nrow(emodel$imatrix)
      whichst <- rep(1:nst, rowSums(emodel$imatrix))
      hcov <- vector(nst, mode="list")
      for (i in 1:nst) {
          if (!any(whichst==i)) hcov[[i]] <-  ~ 1
          else hcov[[i]] <- misccovariates
      }
      hcov
  }

msm.misccovinits2hcovinits <- function(misccovinits, hcovariates, emodel, ecmodel)
  {
      nst <- nrow(emodel$imatrix)
      whichst <- rep(1:nst, rowSums(emodel$imatrix))
      hcovinits <- vector(nst, mode="list")
      for (i in 1:nst) {
          if (is.null(misccovinits)) hcovinits[[i]] <- rep(0, ecmodel$ncovs * rowSums(emodel$imatrix)[i])
          else if (!any(whichst==i)) hcovinits[[i]] <- numeric(0)
          else hcovinits[[i]] <- as.vector(t(sapply(misccovinits, function(x)x[whichst==i])))
      }
      hcovinits
  }

msm.econstr2hconstr <- function(econstr, hmodel)
  {
      constr <- seq(length=hmodel$totpars)
      for (i in unique(econstr)) {
          constr[hmodel$plabs == "p"][econstr == i] <-
            min(constr[hmodel$plabs == "p"][econstr == i])
      }
      match(constr, unique(constr))
  }

print.hmodel <- function(x, ...)
  {
      ci <- (x$fitted && x$foundse)
      cols <- if (ci) c("Estimate","LCL","UCL") else ("")
      if (!x$hidden)
        cat("Non-hidden Markov model\n")
      else {
          cat("Hidden Markov model, ")
          nst <- x$nstates
          cat(nst, "states\n")
          if (x$est.initprobs){
              cat("Initial state occupancy probabilities: ")
              if (x$nipars > 0) {
                  cat("\n")
                  print(x$initprobs)
                  cat("\n")
                  if (any(x$nicovs > 0)) {
                      cat("Covariates on log odds of initial states relative to state 1\n")
                      print(x$icoveffect)
                  }
                  cat("\n")
              }
              else cat(paste(x$initprobs, collapse=","), "\n\n")
          }
          for (i in 1:nst) {
              if (x$mv){
                  cat("State", i, "\n")
                  for (j in 1:x$nout[i]){
                      cat("Outcome", j, "-", x$labels[j,i], "distribution\n")
                      if (x$labels[j,i]=="categorical")
                          pars <- print.hmmcat(x, i, j, mv=TRUE)   ## TESTME
                      else {
                          inds <- x$parstate==i & x$parout==j
                          pars <- as.matrix(x$pars[inds])
                          if (ci) pars <- cbind(pars, matrix(x$ci[inds,], ncol=2))
                          dimnames(pars) <- list(x$plabs[inds], cols)
                      }
                      ## covs not supported for the moment
                      print(pars)
                      cat("\n")
                  }
              } else { 
                  cat("State", i, "-", x$labels[i], "distribution\n")
                  cat("Parameters: \n")
                  if (x$label[i]=="categorical")
                      pars <- print.hmmcat(x, i)
                  else {
                      pars <- as.matrix(x$pars[x$parstate==i])
                      if (ci) pars <- cbind(pars, matrix(x$ci[x$parstate==i ,], ncol=2))
                      dimnames(pars) <- list(x$plabs[x$parstate==i], cols)
                  }
                  if (any(x$ncovs[x$parstate==i] > 0)){
                      coveffs <- as.matrix(x$coveffect[x$coveffstate==i])
                      if (ci) coveffs <- cbind(coveffs, matrix(x$covci[x$coveffstate==i,], ncol=2))
                      rownames(coveffs) <- x$covlabels[x$coveffstate==i]
                      pars <- rbind(pars, coveffs)
                  }
                  print(pars)
                  cat("\n")
              }
          }
      }
  }

print.hmmcat <- function(x, i, j=NULL, mv=FALSE)
{
    inds <- if (mv) x$parstate==i & x$parout==j else x$parstate==i
    pars <- x$pars[inds]
    res <- matrix(pars[3:(2+pars[1])], ncol=1)
    rownames(res) <- paste("P(",seq(length=pars[1]),")",sep="")
    if (x$fitted && x$foundse) {
        ci <- matrix(x$ci[inds,], ncol=2)
        res <- cbind(res,  ci[3:(2+pars[1]),])
        colnames(res) <- c("Estimate","LCL","UCL")
    }
    else colnames(res) <- "prob"
    res
}

msm.form.hconstraint <- function(constraint, hmodel)
  {
      constr <- seq(length=hmodel$totpars)
      for (con in names(constraint)) {
          if ( ! (con %in% c(hmodel$plabs, hmodel$covlabels)))
            stop("parameter \"", con, "\" in hconstraint unknown")
          if (con %in% hmodel$plabs) {
              tc <- constraint[[con]]
              np <- length(tc)
              if (np != sum(hmodel$plabs==con))
                stop("constraint for \"", con, "\" of length ",
                     np, ", should be ", sum(hmodel$plabs==con))
              for (i in unique(tc))
                constr[hmodel$plabs == con][tc == i] <-
                  min(constr[hmodel$plabs == con][tc == i])
          }
      }
      match(constr, unique(constr))
  }

msm.form.hcovconstraint <- function(constraint, hmodel)
  {
      constr <- seq(length=hmodel$ncoveffs)
      for (con in names(constraint)) {
          if ( ! (con %in% c(hmodel$plabs, hmodel$covlabels)))
            stop("parameter \"", con, "\" in hconstraint unknown")
          if (con %in% hmodel$covlabels) {
              tc <- constraint[[con]]
              np <- length(tc)
              if (np != sum(hmodel$covlabels==con))
                stop("constraint for \"", con, "\" of length ",
                     np, ", should be ", sum(hmodel$covlabels==con))
              for (i in unique(tc))
                constr[hmodel$covlabels == con][tc == i] <-
                  min(constr[hmodel$covlabels == con][tc == i])
          }
      }
      match(constr, unique(constr))
  }

## User-supplied range constraints on HMM parameters.

msm.form.hranges <- function(ranges, hmodel)
{
    hranges <- .msm.PARRANGES[hmodel$plabs,]
    if (hmodel$ncoveffs>0) {
        hcovranges <- matrix(rep(c(-Inf, Inf), hmodel$ncoveffs), nrow=hmodel$ncoveffs, byrow=TRUE)
        rownames(hcovranges) <- hmodel$covlabels
        hranges <- rbind(hranges, hcovranges)
    }
    if (!is.null(ranges)) {
        if (!is.list(ranges)) stop("expected \"hranges\" to be a list")
        for (i in names(ranges)) {
            if ( ! (i %in% c(hmodel$plabs, hmodel$covlabels)))
                stop("parameter \"", i, "\" in \"hranges\" unknown")
            ran.default <- hranges[rownames(hranges)==i,,drop=FALSE]
            ran.user <- do.call("cbind", ranges[[i]])
            for (j in seq_len(nrow(ran.user))) {
                if (ran.user[j,1] < ran.default[j,1]) {
                    warning("User-supplied lower bound of ",ran.user[j,1]," for ", i,
                            " less than theoretical minimum of ", ran.default[j,1], ", ignoring")
                    ran.user[j,1] <- ran.default[j,1]
                }
                if (ran.user[j,2] > ran.default[j,2]) {
                    warning("User-supplied upper bound of ",ran.user[j,2]," for ", i,
                            " less than theoretical maximum of ", ran.default[j,2], ", ignoring")
                    ran.user[j,2] <- ran.default[j,2]
                }
                hranges[rownames(hranges)==i,][j,] <- ran.user[j,]
            }
        }
    }
    ## ideally should be strict here if estimated, but allow inits on Inf boundary if not estimated.
    for (i in seq(along=hmodel$pars)){
        if (!in.range(hmodel$pars[i], hranges[i,], strict=FALSE))
            stop("Initial value ", hmodel$pars[i], " of parameter \"", hmodel$plabs[i], "\" outside allowed range ",
                 "[", paste(hranges[i,], collapse=","), "]")
    }
    hranges
}

in.range <- function(x, interval, strict=FALSE) {
    if (!is.numeric(interval) || length(interval)!=2) stop("interval should be a numeric vector of length 2")
    if (!is.numeric(x)) stop("x should be numeric")
    if (strict)
        ( (x > interval[1]) & (x < interval[2]) )
    else
        ( (x >= interval[1]) & (x <= interval[2]) )
}

msm.form.initprobs <- function(hmodel, initprobs, mf){
    npts <- attr(mf,"npts")
    if (!is.null(hmodel$phase.states) && is.null(initprobs)) {
        hmodel$initprobs <- matrix(0, nrow=npts, ncol=hmodel$nstates)
        initstate <- mf$"(state)"[!duplicated(mf$"(subject)")]
        hmodel$initprobs[cbind(1:npts, match(initstate, hmodel$pars))] <- 1
        if (hmodel$est.initprobs)
            warning("Not estimating initial state occupancy probabilities: assuming everyone starts at first phase")
        hmodel$est.initprobs <- FALSE
    }
    if (!hmodel$est.initprobs) {
        if (is.matrix(hmodel$initprobs)) {
            if (nrow(hmodel$initprobs) != npts)
                stop("initial state occupancy probability should have ", npts, " (number of subjects) rows if supplied as a matrix, found ",nrow(hmodel$initprobs))
        }
    }
    hmodel
}
