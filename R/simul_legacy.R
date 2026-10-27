## Original pure-R implementations of sim.msm and simmulti.msm (msm <= 1.8.2),
## kept to check that the C engines simulate from the same distribution.

simmulti_legacy <- function(data,           # data frame with subject, times, covariates...
                         qmatrix,        # intensity matrix
                         covariates=NULL,  # initial values
                         death = FALSE,  # vector of indicators for "death" states, ie absorbing states whose entry time is known exactly,
                                        # but with unknown transient state at previous instant
                         start,         # starting states of the process, defaults to all 1.
                         ematrix = NULL,# misclassification matrix
                         misccovariates = NULL, # covariates on misclassification probabilities
                         hmodel = NULL,  # hidden Markov model formula
                         hcovariates = NULL,   # covariate effects on hidden Markov model response distribution
                         censor.states = NULL,
                         drop.absorb = TRUE
                         )
{

### Check consistency of qmatrix and covariate inits
    nstates <- nrow(qmatrix)
    msm.check.qmatrix(qmatrix)
    qmatrix <- msm.fixdiag.qmatrix(qmatrix)

### Subject, time and state
    if (!("subject" %in% names(data)))
        data$subject <- rep(1, nrow(data))
    if (!("time" %in% names(data)))
        stop("\"time\" column missing from data")
    data <- as.data.frame(data)
    subject <- data[,"subject"]
    time <- data[,"time"]
    if (is.unsorted(subject)){
        warning("Data are not ordered by subject ID and time - output will be ordered.")
        data <- data[order(subject, time),]
    }
    if (any(duplicated(data[,c("subject", "time")]))){
        warning("Data contain duplicated observation times for a subject - removing duplicates.")
        data <- data[!duplicated(data[,c("subject", "time")]), ]
    }
    subject <- data[,"subject"]; time <- data[,"time"];
    cens <- if (any(colnames(data)=="cens")) data[,"cens"] else rep(0, length(subject))
    msm.check.times(time, subject)
    times <- split(time, subject)
    cens <- split(cens, subject)
    n <- length(unique(subject))

### Covariates on intensities
    covnames <- names(covariates)
    ncovs <- length(covnames)
    misscovs <- setdiff(covnames, names(data))
    if (length(misscovs) > 0)
        stop("Covariates ", paste(misscovs, collapse=", "), " not found in data")
    covs <- if (ncovs > 0) lapply(split(data[,covnames], subject), as.matrix) else NULL
    allcovs <- covnames

### Covariates on misclassification
    misccovnames <- names(misccovariates)
    nmisccovs <- length(misccovnames)
    misscovs <- setdiff(misccovnames, names(data))
    if (length(misscovs) > 0)
        stop("Misclassification covariates ", paste(misscovs, collapse=", "), " not found in data")
    misccovs <- if (nmisccovs > 0) lapply(split(data[,setdiff(misccovnames,covnames)], subject), as.matrix) else NULL
    allmisccovs <- misccovnames

### Covariates on HMM
    if (!is.null(hcovariates)) {
        if (is.null(hmodel)) stop("hcovariates specified, but no hmodel")
        hcovnames <- unique(names(unlist(hcovariates)))
        if (length(hcovariates) != nstates)
            stop("hcovariates of length ", length(hcovariates), ", expected ", nstates)
        msm.check.hmodel(hmodel, nstates)
        misscovs <- setdiff(hcovnames, names(data))
        if (length(misscovs) > 0)
            stop("Covariates ", paste(misscovs, collapse=", "), " not found in data")
        hcovs <- lapply(split(data[,setdiff(hcovnames, setdiff(misccovnames,covnames))], subject), as.matrix)
    }
    else hcovs <- hcovnames <- NULL

### Extra variables to return
    extravars <- setdiff(names(data), unique(c("subject","time", "cens", covnames, misccovnames, hcovnames)))
    extradat <- if(length(extravars)>0) split(data[,extravars,drop=FALSE], subject) else NULL

### Starting states
    if (missing(start)) start <- rep(1, n)
    else if (length(start) == 1) start <- rep(start, n)
    else if (length(start) != n)
        stop("Supplied ", length(start), " starting states, expected 1 or ", n)

    nq <- length(qmatrix[qmatrix > 0])
    misspeccovs <- covnames[sapply(covariates, length) != nq]
    if (length(misspeccovs) > 0)
        stop("Initial values for covariates ", paste(misspeccovs, collapse=", "), " should be of length ", nq)
    beta <- do.call("rbind", as.list(covariates))
    ne <- length(ematrix[ematrix > 0])
    misspeccovs <- covnames[sapply(misccovariates, length) != ne]
    if (length(misspeccovs) > 0)
        stop("Initial values for misclassification covariates ", paste(misspeccovs, collapse=", "), " should be of length ", ne)
    beta.misc <- do.call("rbind", as.list(misccovariates))

### Check death argument. Logical values allowed for backwards compatibility
### (TRUE means final state is death, FALSE means no death state)
    statelist <- if (nstates==2) "1, 2" else if (nstates==3) "1, 2, 3" else paste("1, 2, ... ,",nstates)
    if (is.logical(death) && death==TRUE)  {death <- nstates}
    else if (is.logical(death) && death==FALSE) {death <- 0}
    else if (length(setdiff(unique(death), 1:nstates)) > 0)
        stop(paste("Death states indicator contains states not in",statelist))

### Simulate a realisation for each person
    state <- numeric()
    keep.data <- numeric()
    subj <- split(subject, subject)
    subj.num <- match(subject, unique(subject))
    for (pt in 1:n)
    {
        sim.mod <- sim_msm_legacy(qmatrix, max(times[[pt]]), covs[[pt]], beta, times[[pt]], start[pt], min(times[[pt]]))
        obsd <- getobs_legacy(sim.mod, times[[pt]], death, drop.absorb)
        pt.data <- cbind(subj[[pt]][obsd$keep], obsd$time, obsd$state, cens[[pt]][obsd$keep])
        if (!is.null(covnames)) pt.data <- cbind(pt.data, covs[[pt]][obsd$keep,,drop=FALSE])
        if (!is.null(misccovnames)) pt.data <- cbind(pt.data, misccovs[[pt]][obsd$keep,,drop=FALSE])
        if (!is.null(hcovnames)) pt.data <- cbind(pt.data, hcovs[[pt]][obsd$keep,,drop=FALSE])
        if (!is.null(extravars)) pt.data <- cbind(pt.data, extradat[[pt]][obsd$keep,,drop=FALSE])
        inds <- which(subj.num==pt)
        pt.data <- cbind(pt.data, keep=inds[obsd$keep])
        keep.data <- rbind(keep.data, pt.data)
      }
    colnames(keep.data) <- c("subject","time","state","cens",union(union(covnames,misccovnames),hcovnames),extravars,"keep")
    keep.data <- as.data.frame(keep.data)

### Simulate some misclassification or a HMM conditionally on the underlying state
    if (!is.null(ematrix)) {
        if (!all(dim(ematrix) == dim(qmatrix)))
            stop("Dimensions of qmatrix and ematrix should be equal")
        keep.data <- cbind(keep.data, obs=simmisc.msm(keep.data$state, ematrix, beta.misc, keep.data[,misccovnames,drop=FALSE]))
    }
    else if (!is.null(hmodel))
        keep.data <- cbind(keep.data, obs=simhidden.msm(keep.data$state, hmodel, nstates, hcovariates, keep.data[,hcovnames,drop=FALSE]))

### Replace state at censor times by censoring indicators
    censor <- unique(keep.data$cens[keep.data$cens != 0])
    if (is.null(censor.states)) censor.states <- 1:(nstates-1)
    if (!is.null(keep.data$obs))
      keep.data$obs <- ifelse(keep.data$cens > 0 & keep.data$state %in% censor.states, keep.data$cens, keep.data$obs)
    else
      keep.data$state <- ifelse(keep.data$cens > 0 & keep.data$state %in% censor.states, keep.data$cens, keep.data$state)
    keep.data$cens <- NULL

#    attr(keep.data, "keep") <- obsd$keep
    keep.data
}

sim_msm_legacy <- function(qmatrix,   # intensity matrix
                    maxtime,   # maximum time for realisations
                    covs=NULL,     # covariate matrix, nobs rows, ncovs cols
                    beta=NULL,     # matrix of cov effects on qmatrix. ncovs rows, nintens cols.
                    obstimes=0, # times at which time-dependent covariates change
                    start = 1,     # starting state
                    mintime = 0    # time to start from
                    )
{
    ## Keep only times where time-dependent covariates change
    if (!is.null(covs)) {
        covs2 <- collapse_covs_legacy(covs)
        covs <- covs2$covs
        obstimes <- obstimes[covs2$ind]
    }
    else {obstimes <- mintime; covs <- beta <- 0}
    if (is.vector(beta)) beta <- matrix(beta, ncol=length(beta))
    nct <- length(obstimes)
    nstates <- nrow(qmatrix)
    ## Form an array of qmatrices, one for each covariate change-time
    qmatrices <- array(rep(t(qmatrix), nct), dim=c(dim(qmatrix), nct))
    qmatrices[rep(t(qmatrix)>0, nct)] <- qmatrices[rep(t(qmatrix)>0, nct)] * as.numeric(exp(t(beta)%*%t(covs))) # nintens*nobs
    for (i in 1:nct)
        qmatrices[,,i] <- msm.fixdiag.qmatrix(t(qmatrices[,,i]))
    cur.t <- mintime; cur.st <- next.st <- start; rem.times <- obstimes; t.ind <- 1
    nsim <- 0; max.nsim <- 10
    simstates <- simtimes <- numeric(max.nsim) ## allocate memory up-front for simulated outcome
    absorb <- absorbing.msm(qmatrix=qmatrix)
    ## Simulate up to maxtime or absorption
    while (cur.t < maxtime) {
        nsim <- nsim + 1
        cur.st <- next.st
        simstates[nsim] <- cur.st; simtimes[nsim] <- cur.t
        if (cur.st %in% absorb) break;
        rate <- -qmatrices[cur.st,cur.st, t.ind:length(obstimes)]
        nextlag <- rpexp(1, rate, rem.times-rem.times[1])
        cur.t <- cur.t + nextlag
        t.ind <- which.min((cur.t - obstimes)[cur.t - obstimes > 0])
        rem.times <- cur.t
        if (any(obstimes > cur.t))
            rem.times <- c(rem.times, obstimes[(t.ind+1): length(obstimes)])
        cur.q <- qmatrices[,, t.ind]
        next.st <- resample((1:nstates)[-cur.st], size=1, prob = cur.q[cur.st, -cur.st])
        if (nsim > max.nsim) { ## need more memory for simulated outcome, allocate twice as much
            simstates <- c(simstates, numeric(max.nsim))
            simtimes <- c(simtimes, numeric(max.nsim))
            max.nsim <- max.nsim*2
        }
    }
    ## If process hasn't absorbed by the end, then include a censoring time
    if (cur.t >= maxtime) {
        nsim <- nsim+1
        simstates[nsim] <- cur.st
        simtimes[nsim] <- maxtime
    }
    list(states = simstates[1:nsim], times = simtimes[1:nsim], qmatrix = qmatrix)
}


## Drop rows of a covariate matrix which are identical to previous row
## Similar method to R's unique.data.frame

collapse_covs_legacy <- function(covs)
{
    if (nrow(covs)==1) list(covs=covs, ind=1)
    else {
        pcovs <- apply(covs, 1, function(x) paste(x, collapse="\r"))
        lpcovs <- c("\r", pcovs[1:(length(pcovs)-1)])
        ind <- pcovs!=lpcovs
        list(covs=covs[ind,,drop=FALSE], ind=which(ind))
    }
}

### Given a simulated Markov model, get the current state at various observation times
### By default, only keep one observation in the absorbing state

getobs_legacy <- function(sim, obstimes, death=FALSE, drop.absorb=TRUE)
{
    absorb <- absorbing.msm(qmatrix=sim$qmatrix)
                                        # Only keep one observation in the absorbing state
    if (drop.absorb && any(sim$states %in% absorb)) {
        if (any(sim$states %in% death))
            keep <- which(obstimes < max(sim$times))
        else {
            lo <- c(-Inf, obstimes[1:(length(obstimes)-1)])
            keep <- which(lo <= max(sim$times))
        }
    }
    else keep <- 1 : length(obstimes)
    obstimes <- obstimes[keep]
    state <- sim$states[rowSums(outer(obstimes, sim$times, ">="))]
    time <- obstimes
    if (any(sim$states %in% death)) { # Keep the exact death time if required
        state <- c(state, sim$states[sim$states %in% death])
        time <- c(time, sim$times[sim$states %in% death])
        state <- state[order(time)]
        time <- time[order(time)]
        keep <- c(keep, max(keep)+1)
    }
    list(state = state, time = time, keep=keep)
}
