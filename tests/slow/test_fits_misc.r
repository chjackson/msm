context("msm misclassification model fits")

misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=ematrix, death = 4)
miscnew.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q,
                   death = 4, 
                   hmodel=list(
                   hmmCat(prob=c(0.9, 0.1, 0, 0)),
                   hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                   hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()))

test_that("cav misclassification model fit",{
    expect_equal(3951.82919869367, misc.msm$minus2loglik, tol=1e-06)
    if(interactive()) save(misc.msm, file="~/msm/devel/models/1.4/misc.msm.rda")
    stopifnot(isTRUE(all.equal(3951.82919869367, miscnew.msm$minus2loglik, tol=1e-06)))
    expect_equal(miscnew.msm$minus2loglik, misc.msm$minus2loglik)
})

misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,  qmatrix = oneway4.q, ematrix=ematrix, death = 4, misccovariates = ~dage + sex, fixedpars=c(11, 17))
misccovnew.msm <- msm(state ~ years, subject = PTNUM, data = cav,  qmatrix = oneway4.q, death = 4,
                      fixedpars=c(11,17), # sex on state 2|1, 2|3
                      hmodel=list(
                      hmmCat(prob=c(0.9, 0.1, 0, 0)),
                      hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                      hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()),
                      hcovariates=list(~dage + sex, ~dage + sex, ~dage + sex, ~1)
                      )
test_that("cav misclassification model fit with misclassification covariates",{
    expect_equal(3934.69056624596, misccov.msm$minus2loglik, tol=1e-06)
    expect_equal(misccovnew.msm$minus2loglik, misccov.msm$minus2loglik, tol=1e-06)
    if(interactive()) save(misccov.msm, file="~/work/msm/devel/models/1.4/misccov.msm.rda")
})

misccovboth.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=ematrix, death = 4, covariates = ~ I(pdiag=="IHD"),  misccovariates = ~dage, fixedpars=9)
test_that("cav misclassification model fit with both sorts of covariate",{
    expect_equal(3894.49350137067, misccovboth.msm$minus2loglik, tol=1e-06)
    if(interactive()) save(misccovboth.msm, file="~/work/msm/devel/models/1.4/misccovboth.msm.rda")
})

test_that("baseline intensity constraints in misclassification models",{
    miscqc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=4:7,
                    qconstraint = c(1, 2, 1, 2, 3))
    expect_equal(4209.65938095232, miscqc.msm$minus2loglik, tol=1e-06)
})

test_that("baseline misclassification constraints in misclassification models",{
    ematrix2 <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.11, 0),c(0, 0.11, 0, 0),c(0, 0, 0, 0))
    miscec.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                    qmatrix = oneway4.q, ematrix=ematrix2, death = 4, fixedpars=1:5,
                    econstraint = c(1, 1, 2, 2))
    stopifnot(isTRUE(all.equal(4163.71595374043, miscec.msm$minus2loglik, tol=1e-06)))
    e <- ematrix.msm(miscec.msm)
    stopifnot(isTRUE(all.equal(0.969811772846615, e$estimates[1,1], tol=1e-06)))
})

test_that("intensity covariate constraints in misclassification models",{
    misccc.msm <- msm(state ~ years, subject = PTNUM, data = cav, fixedpars=c(1:5, 9:12),
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4,
                    covariates = ~ sex, covinits=list(sex=c(0, 0, 0.1, 0, 0)),
                    constraint = list(sex = c(1, 2, 1, 2, 3))  )
    stopifnot(isTRUE(all.equal(4276.38456321703, misccc.msm$minus2loglik, tol=1e-06)))
    q <- qmatrix.msm(misccc.msm, covariates=0)
    stopifnot(isTRUE(all.equal(-0.176204765732341, q$estimates[1,1])))
})

test_that("misclassification covariate constraints in misclassification models",{
    miscecc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                       qmatrix = oneway4.q, ematrix=ematrix, death = 4,
                       misccovariates = ~dage + sex, fixedpars=c(1:5, 7,  11,12, 15),
                       misccovinits = list(dage=c(0.01,0.01,0.001,0.001), sex=c(0.0131,0.0132,0.0133,0.0134)),
                       miscconstraint = list(dage = c(1, 1, 2, 2)), control=list(maxit=10000))
    stopifnot(isTRUE(all.equal(4094.07318840828, miscecc.msm$minus2loglik, tol=1e-06)))
    e <- ematrix.msm(miscecc.msm)
    stopifnot(isTRUE(all.equal(0.988468090015414, e$estimates[1,1], tol=1e-04)))
    e <- ematrix.msm(miscecc.msm, covariates=0)
    stopifnot(isTRUE(all.equal(0.986856424900654, e$estimates[1,1], tol=1e-06)))
})

test_that("estimating initprobs",{
    nsubj <- 50; nobspt <- 6
    sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length.out=nobspt),
                         x = rnorm(nsubj*nobspt), y = rnorm(nsubj*nobspt)* 5 + 2 )
    (three.q <- rbind(c(0, exp(-3), exp(-6)), c(0, 0, exp(-3)), c(0, 0, 0)))
    ematrix3 <- rbind(c(0, 0.1, 0), c(0.1, 0, 0), c(0,0,0))
    initprobs <- c(0.5, 0.5, 0)
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, ematrix=ematrix3, start=sample(1:3, 50, prob=initprobs, replace=TRUE))
    miscip.msm <- msm(obs ~ time, subject = subject, data = sim2.df, qmatrix = three.q, ematrix=ematrix3, initprobs=c(0.4, 0.5, 0.05), fixedpars=7, est.initprobs=TRUE, control=list(maxit=10000))
    stopifnot(miscip.msm$hmodel$initprobs["State 2","LCL"] < 0.5 && 0.5 < miscip.msm$hmodel$initprobs["State 2","UCL"])
})

test_that("estimating covariate effects on initprobs",{
    nsubj <- 500; nobspt <- 6
    set.seed(22061976)
    sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt),
                         time = seq(0, 20, length.out=nobspt),
                         x = rnorm(nsubj*nobspt))
    (three.q <- (rbind(c(0, exp(-3), exp(-6)),
                       c(0, 0, exp(-3)),
                       c(0, 0, 0))))
    ematrix3 <- rbind(c(0, 0.1, 0), c(0.1, 0, 0), c(0,0,0))
    ip.base <- c(0.2, 0.3, 0.5)
    beta <- 2 # covariate effect on log(ip2/ip1) and log(ip3/ip1)
    ipl2 <- exp(log(ip.base[2]/ip.base[1]) + beta*sim.df$x[sim.df$time==0])
    ipl3 <- exp(log(ip.base[3]/ip.base[1]) + beta*sim.df$x[sim.df$time==0])
    initprobs <- cbind(1, ipl2, ipl3)/(1 + ipl2 + ipl3)
    start <- numeric(nsubj)
    for (i in 1:nsubj)
      start[i] <- sample(1:3, 1, prob=initprobs[i,], replace=TRUE)
    sim2.df <- simmulti.msm(sim.df, qmatrix=three.q, ematrix=ematrix3,
                            start=start, covariates=list(x=c(0,0,0)))
    misc.msm <- msm(obs ~ time, subject = subject, data = sim2.df,
                    qmatrix = three.q, ematrix=ematrix3, center=FALSE,
                    initcovariates = ~ x,
                    initcovinits = list(x=c(0.1,0.01)),
                    est.initprobs=TRUE, fixedpars=1:4)
    stopifnot(misc.msm$hmodel$initprobs["State 2","LCL"] < ip.base[2] &&
              ip.base[2] < misc.msm$hmodel$initprobs["State 2","UCL"])
    stopifnot(misc.msm$hmodel$initprobs["State 3","LCL"] < ip.base[3] &&
              ip.base[3] < misc.msm$hmodel$initprobs["State 3","UCL"])
    stopifnot(misc.msm$hmodel$icoveffect["x, State 2","LCL"] < beta
              && beta < misc.msm$hmodel$icoveffect["x, State 2","UCL"])
    stopifnot(misc.msm$hmodel$icoveffect["x, State 3","LCL"] < beta
              && beta < misc.msm$hmodel$icoveffect["x, State 3","UCL"])

    ## with structural zeros (Jeffrey Eaton & Tara Mangal's bug fix)
    set.seed(1)
    ip.base <- c(0.2, 0, 0.8) # zero in middle 
    beta <- 2 # covariate effect on log(ip2/ip1) and log(ip3/ip1)
    ipl2 <- 0
    ipl3 <- exp(log(ip.base[3]/ip.base[1]) + beta*sim.df$x[sim.df$time==0])
    initprobs <- cbind(1, ipl2, ipl3)/(1 + ipl2 + ipl3)
    start <- numeric(nsubj)
    for (i in 1:nsubj)
      start[i] <- sample(1:3, 1, prob=initprobs[i,], replace=TRUE)
    sim2.df <- simmulti.msm(sim.df, qmatrix=three.q, ematrix=ematrix3,
                            start=start, covariates=list(x=c(0,0,0)))
    misc.msm <- msm(obs ~ time, subject = subject, data = sim2.df,
                    qmatrix = three.q, ematrix=ematrix3, center=FALSE,
                    initprobs=c(0.5, 0, 0.5),
                    initcovariates = ~ x,  initcovinits = list(x=0),
                    est.initprobs=TRUE, fixedpars=1:4)
    stopifnot(misc.msm$hmodel$initprobs["State 3","LCL"] < ip.base[3] &&
              ip.base[3] < misc.msm$hmodel$initprobs["State 3","UCL"])
    stopifnot(misc.msm$hmodel$icoveffect["x, State 3","LCL"] < beta
              && beta < misc.msm$hmodel$icoveffect["x, State 3","UCL"])

    ## structural zeros: zero at end
    set.seed(1)
    ip.base <- c(0.2, 0.8, 0)
    beta <- 2 # covariate effect on log(ip2/ip1) and log(ip3/ip1)
    ipl2 <- exp(log(ip.base[2]/ip.base[1]) + beta*sim.df$x[sim.df$time==0])
    ipl3 <- 0
    initprobs <- cbind(1, ipl2, ipl3)/(1 + ipl2 + ipl3)
    start <- numeric(nsubj)
    for (i in 1:nsubj)
      start[i] <- sample(1:3, 1, prob=initprobs[i,], replace=TRUE)
    sim2.df <- simmulti.msm(sim.df, qmatrix=three.q, ematrix=ematrix3,
                            start=start, covariates=list(x=c(0,0,0)))
    misc.msm <- msm(obs ~ time, subject = subject, data = sim2.df,
                    qmatrix = three.q, ematrix=ematrix3, center=FALSE,
                    initprobs=c(0.5, 0.5, 0),
                    initcovariates = ~ x,  initcovinits = list(x=0),
                    est.initprobs=TRUE, fixedpars=1:4)
    stopifnot(misc.msm$hmodel$initprobs["State 2","LCL"] < ip.base[2] &&
              ip.base[2] < misc.msm$hmodel$initprobs["State 2","UCL"])
    stopifnot(misc.msm$hmodel$icoveffect["x, State 2","LCL"] < beta
              && beta < misc.msm$hmodel$icoveffect["x, State 2","UCL"])
})

test_that("obstrue",{
    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, obstrue=firstobs,
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE)
    stopifnot(all.equal(misc.msm$minus2loglik, 4165.84711809003))
### test against dummy covariate hack
    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE, center=FALSE,
                    misccovariates=~firstobs, misccovinits = list(firstobs=rep(-1e+06,4)))
    stopifnot(all.equal(misc.msm$minus2loglik, 4165.84711809003, tol=1e-06))
})

context("outputs of misclassification models")

test_that("ematrix.msm",{
    e <- ematrix.msm(misc.msm)
    expect_equal(0.00766164690017842, e$estimates[1,2], tol=1e-06)
    expect_equal(0.00334014173130528, e$SE[1,2], tol=1e-06)
    expect_equal(0.00325308687951399, e$L[1,2], tol=1e-06)
    expect_equal(0.0179371415700995, e$U[1,2], tol=1e-06)
    e <- ematrix.msm(misccov.msm)
    expect_equal(0.005888168, e$estimates[1,2], tol=1e-06)
    expect_equal(0.004122881, e$SE[1,2], tol=1e-06)
    e <- ematrix.msm(misccov.msm, covariates=0)
    expect_equal(0.00304673983030711, e$estimates[1,2], tol=1e-06)
    expect_equal(0.00480361545354035, e$SE[1,2], tol=1e-06)
    expect_equal(0.000266862963906964, e$L[1,2], tol=1e-06)
    expect_equal(0.116160648623402, e$U[1,2], tol=1e-06)
    e <- ematrix.msm(misccov.msm, covariates=list(dage=50, sex=0))
    expect_equal(0.00975889769191679, e$estimates[1,2], tol=1e-06)
    expect_error(ematrix.msm("Foo"), "expected .+ msm model")
    expect_warning(ematrix.msm(misc.msm, covariates=list(foo=1)), "no covariates")
    expect_warning(ematrix.msm(misccov.msm, covariates=list(foo=1)), "Covariate .+ unknown")
})

test_that("viterbi.msm",{
    vit <- viterbi.msm(misc.msm)[viterbi.msm(misc.msm)$subject==100063,]
    expect_equal(c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2), vit$fitted, tol=1e-06)
    ## posterior probabilities
    expect_true(all(vit$pstate >= 0-0.000001 & vit$pstate <= 1.000001))
    expect_equal(rowSums(vit$pstate), rep(1, nrow(vit)), tol=1e-06)
    expect_equal(vit$fitted, apply(vit$pstate, 1, which.max))
    expect_error(viterbi.msm("foo"), "expected .+ msm model")
})

test_that("viterbi.msm and ematrix.msm with non-HMM",{
    cav.msm <- msm(state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, fixedpars=TRUE)
    vit <- viterbi.msm(cav.msm)
    stopifnot(all.equal(vit$observed, vit$fitted))
    expect_null(ematrix.msm(cav.msm)) # not a misc model
})

test_that("odds.msm",{
    odds <- odds.msm(misccov.msm)
    expect_equal(0.925025938925706, odds$dage[2,2], tol=1e-06)
    expect_equal(0.939115450192996, odds$dage[1,2], tol=1e-06)
    expect_equal(5.30030000947669, odds$sex[2,3], tol=1e-06)
    expect_equal(1, odds$sex[4,3], tol=1e-06)
    expect_error(odds.msm("foo"), "expected .+ msm model")
    expect_error(odds.msm(misccov.msm, odds.scale="foo"), "non-numeric")
    expect_error(odds.msm(misccov.msm, odds.scale=c(1,2,3,4)), "scale of length 4, expected 2")
    expect_error(odds.msm(misccov.msm, odds.scale=c(1,2,3)), "scale of length 3, expected 2")
    expect_equal(odds.msm(misccov.msm, odds.scale=2)$dage[1,1], odds.msm(misccov.msm)$dage[1,1]^2)
    expect_equal(odds.msm(misccov.msm, odds.scale=c(2,3))$sex[1,1], odds.msm(misccov.msm)$sex[1,1]^3)
})

test_that("non misclassification-specific output functions on misclassification models",{
    q <- qmatrix.msm(misccov.msm)
    expect_equal(0.227327094280783, q$estimates[2,3], tol=1e-06)
    soj <- sojourn.msm(misccov.msm)
    expect_equal(6.81099600426226, soj[1,1])
    p <- pmatrix.msm(misccov.msm, 10)
    expect_equal(0.119622293176754, p[1,3], tol=1e-03)
    p <- prevalence.msm(misccov.msm)
    expect_equal(158, p$Observed[5,4], tol=1e-03)
    expect_equal(134.412855042637, p$Expected[5,4], tol=1e-03)
    tot <- totlos.msm(misccov.msm)
    expect_equal(c(6.81099600481148, 2.84328019082745, 2.08683929661159), as.numeric(tot[1:3]), tol=1e-06)
})

test_that("outputs of hmmCat misclassification models",{
    expect_null(ematrix.msm(miscnew.msm))
    pars <- miscnew.msm$hmodel$pars
    expect_equivalent(pars[names(pars) %in% c("pbase","p","p0")], as.numeric(t(misc.msm$Ematrices$baseline[1:3,])))
    vit <- viterbi.msm(misc.msm)
    vitnew <- viterbi.msm(miscnew.msm)
    expect_equal(vit, vitnew)
    expect_error(odds.msm(misccovnew.msm),"Requires a misclassification model specified with ematrix")
    expect_equivalent(odds.msm(misccov.msm)$dage[,1], exp(misccovnew.msm$hmodel$coveffect[misccovnew.msm$hmodel$covlabels=="dage"]))
    expect_equivalent(odds.msm(misccov.msm)$sex[,1], exp(misccovnew.msm$hmodel$coveffect[misccovnew.msm$hmodel$covlabels=="sex"]))
    expect_equal(qmatrix.msm(misccovnew.msm,covariates=0), qmatrix.msm(misccov.msm,covariates=0))
})

test_that("viterbi.msm bug with fixedpars",{
    v <- viterbi.msm(misc.msm)
    for (pt in unique(cav$PTNUM)){
        subs <- cav[cav$PTNUM==pt,]
        x <- v[v$subject==pt,]
        miscfix.msm <- msm(state ~ years, subject = PTNUM, data = subs, qmatrix = qmatrix.msm(misc.msm)$estimates, ematrix=ematrix.msm(misc.msm)$estimates, death = 4, fixedpars=TRUE)
        y <- viterbi.msm(miscfix.msm)
        stopifnot(all(x$fitted==y$fitted))
    }
})
