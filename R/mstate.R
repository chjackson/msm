## Function to convert data for a msm model fit to data for a coxph model fit

msm2Surv <- function(data, # data frame
                     subject, time, state, # names of subject, time and state variables (character)
                     covs=NULL, # names of covariates (character vector)
                     Q # transition intensity matrix.  should be zero where transitions are disallowed.
                     ) {
    if (missing(subject)) stop("subject variable not given")
    if (missing(time)) stop("time variable not given")
    if (missing(state)) stop("state variable not given")
    fpt <- !duplicated(data[,subject]) # indicator for patient's first observation
    lpt <- !duplicated(data[,subject], fromLast=TRUE) # ... last observation ...
    nev <- nrow(data[!lpt,])
    ## Data frame of observed events
    ev <- data.frame(id=data[!lpt,subject], from=data[!lpt,state],to=data[!fpt,state],
                     Tstart=data[!lpt,time], Tstop=data[!fpt,time], time=data[!fpt,time]-data[!lpt,time],
                     status=rep(1,nev))
    if (is.null(covs)) covs <- setdiff(colnames(data), c(subject, time, state))
    ## rename any covariates which clash with standard names in the returned data
    for (i in c("id", "from", "to","Tstart","Tstop","time","status","trans"))
        covs[covs==i] <- colnames(data)[colnames(data)==i] <- paste(i, ".2", sep="")
    for (i in covs)
        ev[,i] <- data[!lpt, i]
    neq <- sum(ev$Tstart == ev$Tstop)
    if (neq > 0) {
        warning("Omitting ",neq, " rows with two observations at the same time")
        ev <- ev[ev$Tstart < ev$Tstop,]
    }
    diag(Q) <- 0; Q[Q>0] <- 1
    if (is.null(rownames(Q))) rownames(Q) <- 1:nrow(Q)
    if (is.null(colnames(Q))) colnames(Q) <- 1:ncol(Q)
    Qf <- Q[ev$from,]
    Qf[cbind(1:nrow(Qf), ev$to)] <- 0
    nto <- rowSums(Qf)
    ncens <- sum(nto)
    cto <- which(t(Qf)==1,arr.ind=TRUE)[,1]
    ## Data frame of censored events
    cens <- data.frame(id=rep(ev$id, nto), from=rep(ev$from, nto), to=cto,
                       Tstart=rep(ev$Tstart, nto), Tstop=rep(ev$Tstop, nto),
                       time=rep(ev$Tstop, nto) - rep(ev$Tstart, nto),
                       status=rep(0, ncens))
    for (i in covs)
        cens[,i] <- rep(ev[,i], nto)
    surv <- rbind(ev, cens)
    surv <- surv[order(surv$id, surv$Tstart, surv$to),]
    surv <- surv[!(surv$from==surv$to),]
    Qi <- t(Q); Qi[Qi==1] <- seq_along(which(t(Q)==1)); Qi <- t(Qi)
    surv$trans <- Qi[cbind(surv$from,surv$to)]
    rownames(surv) <- NULL
    tmat <- t(Q)
    tmat[t(Q)==1] <- seq(along=tmat[t(Q)==1])
    tmat <- t(tmat); tmat[Q==0] <- NA
    names(dimnames(tmat)) <- c("from","to")
    attr(surv, "trans") <-  tmat
    class(surv) <- c("msdata","data.frame")
    surv
}
