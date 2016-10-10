## forget about data.orig, this was just without pci imputation
## not supported whichcovh.orig
## could use unname to remove names on stuff

recreate.olddata <- function(x) {

    x$data <- expand.data(x)
    mf <- x$data$mf; mf.agg <- x$data$mf.agg; mm.cov <- x$data$mm.cov; mm.cov.agg <- x$data$mm.cov.agg;
    mm.mcov <- x$data$mm.mcov;    mm.hcov <- x$data$mm.hcov; mm.icov <- x$data$mm.icov

    get.covdata <- function(mm, mf) {
        list(covlabels=colnames(mm)[-1],
             ncovs=ncol(mm)-1,
             covfactor=sapply(mf, is.factor),
             covfactorlevels=lapply(mf, levels),
             covmat=mm,
             covmat.orig=mf,
             covlabels.orig=colnames(mf) ## forget about covrows.kept
             )
    }

    covdata <- get.covdata(mm.cov, mf)
    misccovdata <- if (x$emodel$misc) get.covdata(mm.mcov, mf) else NULL
    if (x$hmodel$hidden) {
        hcovdata <- list()
        for (i in seq(x$qmodel$nstates))
            hcovdata[[i]] <- get.covdata(mm.hcov[[i]], mf)
        icovdata <- get.covdata(mm.icov, mf)
        whichcovh <- lapply(mm.hcov, function(x)match(colnames(x)[-1], colnames(mm)[-1]))
        for (i in seq(x$qmodel$nstates)) hcovdata[[i]]$whichcov <- whichcovh[[i]]
    } else hcovdata <- icovdata <- NULL

    ## TESTME - watch for factors, interactions and intercept, reduce if no covs.
    mmh <- if (!is.null(mm.hcov)) do.call("cbind",mm.hcov) else NULL
    mm <- cbind(mm.cov, mm.icov, mmh)
    mm <- mm[,unique(colnames(mm))]
    covdata$whichcov <- match(colnames(mm.cov)[-1], colnames(mm)[-1])
    covdata$whichcov.orig <- match(attr(mf, "covnames.q"), attr(mf, "covnames"))
    misccovdata$whichcov <- match(colnames(mm.mcov)[-1], colnames(mm)[-1]) # names should be in mm.hcov
    icovdata$whichcov <- match(colnames(mm.icov)[-1], colnames(mm)[-1])

    ret <- list(
        fromstate = model.extract(mf.agg, "fromstate"),
        tostate = model.extract(mf.agg, "tostate"),
        timelag = model.extract(mf.agg, "timelag"),
        nocc = model.extract(mf.agg, "nocc"),
        whicha = model.extract(mf.agg, "whicha"),
        noccsum = model.extract(mf.agg, "noccsum"),
        obstype = if (!is.null(mf.agg)) model.extract(mf.agg, "obstype") else model.extract(mf, "obstype"),
        covmat = mm.cov.agg[,-1,drop=FALSE],
        covdata = covdata, misccovdata = misccovdata, hcovdata = hcovdata, icovdata = icovdata,
        npts = length(unique(model.extract(mf, "subject"))),
        covlabels = names(attr(mf, "covmeans"))[-1],
        covmeans = attr(mf, "covmeans")[-1],
        nobs = if (!is.null(mf.agg)) nrow(mf.agg) else nrow(mf),
        ntrans = sum(duplicated(model.extract(mf, "subject"))),
        time = mf[,2],
        state = mf[,1],
        subject = model.extract(mf, "subject"),
        n = nrow(mf),
        obstype.obs = model.extract(mf, "obstype"),
        firstobs = c(which(!duplicated(model.extract(mf, "subject"))), nrow(mf)+1),
        obstrue = model.extract(mf, "obstrue"),
        pci.imp = model.extract(mf, "pci.imp"),
        cov = mm[,-1,drop=FALSE],
        cov.orig = mf[,attr(mf,"covnames")],
        covlabels.orig = attr(mf,"covnames")
        )
    ret
}
