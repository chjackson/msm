
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
    class(res) <- "msm.prevalence"
    if (plot) plot.prevalence.msm(x, mintime=min(times), maxtime=max(times), timezero=timezero, initstates=initstates,
                                  interp=interp, censtime=censtime, covariates=covariates, misccovariates=misccovariates,
                                  piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates, ...)
    res
}

#' @export
print.msm.prevalence <- function(x,...){
  print(unclass(x))
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
#' @export
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
