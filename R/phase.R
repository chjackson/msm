msm.phase2qmodel <- function(qmodel, phase.states, inits, qconstraint, analyticp, use.expm){
    if (any(!(phase.states %in% 1:qmodel$nstates))) stop("phase.states should be in 1,...,",qmodel$nstates)
    markov.states <- setdiff(1:qmodel$nstates, phase.states)
    nst <- 2*length(phase.states) + length(markov.states)
    reps <- rep(1, qmodel$nstates)
    reps[phase.states] <- 2
    pars <- rep(1:qmodel$nstates, reps) # index of each new state in old states

    if (!is.null(inits)){
        if (!is.list(inits)) stop("phase.inits should be a list")
        if (!length(inits)==length(phase.states))
            stop(sprintf("phase.inits of length %d, but there are %d phased states", length(inits), length(phase.states)))
        for (i in seq(along=inits)){
            if (!length(inits[[i]])==2) stop(sprintf("phase.inits[[%d]] list of length %d, should be 2", i, length(inits[[i]])))
            if (is.null(names(inits[[i]]))) names(inits[[i]]) <- c("trans","exit")
            if (is.vector(inits[[i]]$trans)) inits[[i]]$trans <- matrix(inits[[i]]$trans, nrow=1)
            if (is.vector(inits[[i]]$exit)) inits[[i]]$exit <- matrix(inits[[i]]$exit, nrow=1)
            if (length(inits[[i]]$trans)!=1) stop(sprintf("phase.inits[[%d]]$trans of length %d, should be 1", i, length(inits[[i]]$trans)))
            if (ncol(inits[[i]]$exit)!=2) stop(sprintf("phase.inits[[%d]]$exit has %d columns, should be 2", i, length(inits[[i]]$exit)))
        }
    }
    ## TODO warn for identical exit rates

    ## Form new transition intensity matrix
    qmatrix.new <- matrix(0, nrow=nst, ncol=nst)
    rn <- rownames(qmodel$qmatrix)[pars]; cn <- colnames(qmodel$qmatrix)[pars]
    pp <- pars %in% phase.states
    rn[pp] <- paste(rn[pp], " [P", rep(1:2, length(phase.states)), "]", sep="")
    cn[pp] <- paste(cn[pp], " [P", rep(1:2, length(phase.states)), "]", sep="")
    ## also form short state names here (used in viterbi)
    sn <- pars; sn[pp] <- paste0(sn[pp], rep(c("P1","P2"), length(phase.states)))
    dimnames(qmatrix.new) <- list(rn, cn)

    ## Transition rates between non-phased states
    mpars <- which(pars %in% markov.states)
    qmatrix.new[mpars,mpars] <- qmodel$qmatrix[markov.states,markov.states,drop=FALSE]
    ## Phase entry transition rates
    phase1 <- which(pars%in%phase.states & !duplicated(pars))
    qmatrix.new[mpars,phase1] <- qmodel$qmatrix[markov.states,phase.states,drop=FALSE]

    ## Unsure if we need to keep all of these.  TODO document them 
    qaux <- list(phase.states=phase.states, markov.states=markov.states,
                 phase.reps=reps, phase.pars=pars,
                 oldstates=pars, phase.labs=sn,
                 imatrix.orig=qmodel$imatrix, qmatrix.orig=qmodel$qmatrix,
                 pdests=list(), pdests.orig=list(),
                 phase1.ind=numeric(), phase2.ind=numeric()) 
    for (i in seq(along=phase.states)){
        ## possible destination states from current phase state in
        ## original and expanded model respectively
        dests <- which(qmodel$imatrix[phase.states[i],]==1)
        dests.new <- which(pars %in% dests & !duplicated(pars))
        erates <- qmodel$qmatrix[phase.states[i],dests]
        if (is.null(inits)) {
            ## Phase exit transition rates
            erates1 <- erates*0.8
            erates2 <- erates*1.2
        }
        else {
            if (nrow(inits[[i]]$exit) != length(dests)){
                plural <- if (nrow(inits[[i]]$exit) > 1) "s" else ""
                stop(sprintf("phase.inits[[%d]]$exit has %d row%s, but there are %d exit states from this state", i, nrow(inits[[i]]$exit), plural, length(dests)))
            }
            erates1 <- inits[[i]]$exit[,1]; erates2 <- inits[[i]]$exit[,2]
        }
        phase1 <- which(pars==phase.states[i])[1]
        phase2 <- which(pars==phase.states[i])[2]
        qmatrix.new[phase1,dests.new] <- erates1
        qmatrix.new[phase2,dests.new] <- erates2
        ## default to 0.5 for p2
        if (is.null(inits)) {
            trans <- 0.5/sum(erates1)
        }
        else trans <- inits[[i]]$trans
        qmatrix.new[phase1,phase2] <- trans
        qaux$pdests.orig[[i]] <- dests; qaux$pdests[[i]] <- dests.new
        qaux$phase1.ind[i] <- phase1; qaux$phase2.ind[i] <- phase2
    }
    q.new <- msm.form.qmodel(qmatrix=qmatrix.new, qconstraint=qconstraint,
                             analyticp=analyticp, use.expm=use.expm, phase.states=NULL)
    q.new <- c(q.new, qaux)
    q.new
}

msm.phase2hmodel <- function(qmodel, hmodel){   
    if (is.null(hmodel)) {
        hmodel <- vector(qmodel$nstates, mode="list")
        for (i in 1:qmodel$nstates) hmodel[[i]] <- hmmIdent(qmodel$phase.pars[i])
    }
    ## dummy initprobs here: actually defined later in msm.form.initprobs, since needs knowledge of the data
    hmodel <- c(msm.form.hmodel(hmodel=hmodel,
                                initprobs=c(1,rep(0,qmodel$nstates-1)), 
                                est.initprobs=FALSE),
                list(phase.states=qmodel$phase.states))
    hmodel
}

## Convert parameters of fitted phase-type model to mixture representation

phasemeans.msm <- function(x, covariates="mean", ci=c("none","normal","bootstrap"), cl=0.95, B=1000, cores=NULL){
    ps <- x$qmodel$phase.states
    Q <- qmatrix.msm(x, ci="none", covariates=covariates)
    res <- matrix(nrow=length(ps), ncol=3)
    pnames <- c("Short stay mean", "Long stay mean", "Long stay probability")
    rownames(res) <- rownames(x$qmodel$imatrix.orig)[x$qmodel$phase.states]
    colnames(res) <- pnames
    for (i in seq_along(ps)){
        p1 <- x$qmodel$phase1.ind[i]; p2 <- x$qmodel$phase2.ind[i]
        lam1 <- Q[p1, p2]
        mu1 <- sum(Q[p1, -c(p1, p2)])
        mu2 <- sum(Q[p2, -p2])
        mean1 <- 1/(lam1 + mu1)
        mean2 <- mean1 + 1/mu2
        prob2 <- lam1*mean1
        res[i,] <- c(mean1, mean2, prob2)
    }
    ci <- match.arg(ci)
    p.ci <- switch(ci,
                   normal = phasemeans.normci.msm(x=x, covariates=covariates, cl=cl, B=B),
                   bootstrap = phasemeans.ci.msm(x=x, covariates=covariates, cl=cl, B=B, cores=cores),
                   none = NULL)
    res <- if (ci=="none") res else list(estimates = res, L=p.ci[,,1], U=p.ci[,,2])
    class(res) <- "msm.est"
    res
}

d2phase <- function(x, l1, mu1, mu2, log=FALSE){
    t <- x
    ret <- numeric(length(t))
    ret[t<0] <- 0
    ret[l1<0 | mu1<0 | mu2<0] <- NaN
    ind <- (l1>=0 & mu1>=0 & mu2>=0 & t>=0)
    l1 <- rep(l1, length=length(t)); mu1 <- rep(mu1, length=length(t)); mu2 <- rep(mu2, length=length(t)); 
    if(any(ind)) { 
        ## same as MatrixExp(Q*t)[1,] %*% c(mu1, mu2, 0)
        if (log) {
            tmp <- ifelse(l1+mu1==mu2,
                          -(l1+mu1)*t + log(mu1 + (l1+mu1)*l1*t),
                          log((-(l1+mu1)*exp(-(l1+mu1)*t)*(-mu1+mu2) + mu2*l1*exp(-mu2*t))/(l1+mu1-mu2)))
        } else {
            tmp <- ifelse(l1+mu1==mu2,
                          exp(-(l1+mu1)*t)*(mu1 + (l1+mu1)*l1*t),
                          (-(l1+mu1)*exp(-(l1+mu1)*t)*(-mu1+mu2) + mu2*l1*exp(-mu2*t))/(l1+mu1-mu2))
        }
        ret[ind] <- tmp[ind]
    }
    ret
}

p2phase <- function(q, l1, mu1, mu2, lower.tail=TRUE, log.p=FALSE){
    t <- q
    ret <- numeric(length(t))
    ret[t<0] <- 0
    ret[l1<0 | mu1<0 | mu2<0] <- NaN
    ind <- (l1>=0 & mu1>=0 & mu2>=0 & t>=0)
    if(any(ind)) { 
        if (l1+mu1==mu2) {
            if (!lower.tail && log)
                tmp <- -(l1+mu1)*t + log(1 + l1*t)    
            else {
                tmp <- 1 - exp(-(l1+mu1)*t)*(1 + l1*t)
                if (!lower.tail) tmp <- 1 - tmp
                if (log.p) tmp <- log(tmp)
            }
        }
        else { 
            ## same as MatrixExp(Q*t)[1,3]
            tmp <- 1 + exp(-(l1+mu1)*t)*(-mu1+mu2)/(l1+mu1-mu2) - l1*exp(-mu2*t)/(l1+mu1-mu2)
            if (!lower.tail) tmp <- 1 - tmp
            if (log.p) tmp <- log(tmp)
        }
        ret[ind] <- tmp[ind]        
    }
    ret
}

q2phase <- function(p, l1, mu1, mu2, lower.tail=TRUE, log.p=FALSE){
    qgeneric(p2phase, p=p, l1=l1, mu1=mu1, mu2=mu2, lower.tail=lower.tail, log.p=log.p)
}

r2phase <- function(n, l1, mu1, mu2){
    if (length(n) > 1) n <- length(n)
    ret <- numeric(n)
    ret[l1<0 | mu1<0 | mu2<0] <- NaN
    ind <- (l1>=0 & mu1>=0 & mu2>=0)
    if(any(ind)) {
        l1 <- rep(l1, length=n); mu1 <- rep(mu1, length=n); mu2 <- rep(mu2, length=n)
        l1 <- l1[ind]; mu1 <- mu1[ind]; mu2 <- mu2[ind]; n <- sum(ind)
        ptrans <- l1/(l1+mu1)
        ret[ind] <- rexp(n,l1+mu1) + rbinom(n, 1, ptrans)*rexp(n, mu2)
    }
    ret
}

h2phase <- function(x, l1, mu1, mu2, log=FALSE) {
    t <- x
    ret <- numeric(length(x))
    ret[t<0] <- 0
    ret[l1<0 | mu1<0 | mu2<0] <- NaN
    ind <- (l1>=0 & mu1>=0 & mu2>=0 & t>=0)
    if(any(ind)) {
        if (l1+mu1==mu2) {
            tmp <- (mu1 + (l1+mu1)*l1*t) / (1 + l1*t) # clearer formula
            if (log) tmp <- log(tmp)
        }
        else {
            if (log)
                tmp <- d2phase(t,l1,mu1,mu2,log=TRUE) - p2phase(t,l1,mu1,mu2,lower.tail=FALSE,log.p=TRUE)
            else
                tmp <- d2phase(t,l1,mu1,mu2) / p2phase(t,l1,mu1,mu2,lower.tail=FALSE)
        }
        ret[ind] <- tmp[ind]
    }
    ret
}
