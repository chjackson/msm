## depends on psor.msm

context("analytic derivatives of likelihood")

test_that("derivatives by subject: sum to overall derivative",{
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), fixedpars=FALSE)
    q.mle <- psor.msm$paramdata$opt$par
    deriv.overall <- deriv.msm(q.mle, expand.data(psor.msm), psor.msm$qmodel, psor.msm$qcmodel, psor.msm$cmodel, psor.msm$hmodel, psor.msm$paramdata)
    deriv.subj <- Ccall.msm(q.mle, do.what="deriv.subj", expand.data(psor.msm), psor.msm$qmodel, psor.msm$qcmodel, psor.msm$cmodel, psor.msm$hmodel, psor.msm$paramdata)
    expect_equal(deriv.overall, colSums(deriv.subj))
})

options(msm.test.analytic.derivatives=TRUE)
err <- 1e-04

test_that("analytic derivatives match numeric",{
    cav.msm <- msm(state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE)
    expect_lt(deriv_error(cav.msm), err)
    cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = FALSE, fixedpars=TRUE)
    expect_lt(deriv_error(cav.msm), err)
    cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qconstraint = c(1,1,2,2,2,3,3),  qmatrix = twoway4.q, death = FALSE, fixedpars=TRUE)
    expect_lt(deriv_error(cav.msm), err)
    psor.0.q <- rbind(c(0,0.1,0,0),c(0,0,0.2,0),c(0,0,0,0.3),c(0,0,0,0))
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.0.q, fixedpars=TRUE)
    expect_lt(deriv_error(psor.msm), err)
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.0.q, covariates = ~ollwsdrt+hieffusn,
                    constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), fixedpars=TRUE)
    expect_lt(deriv_error(psor.msm), err)
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.0.q, covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), death=TRUE, fixedpars=TRUE)
    expect_lt(deriv_error(psor.msm), err)
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.0.q, covariates = ~ollwsdrt+hieffusn, death=TRUE, fixedpars=TRUE)
    expect_lt(deriv_error(psor.msm), err)

    msmtest5 <- msm(state ~ time, qmatrix = fiveq, subject = ptnum, data = bos, exacttimes=TRUE, fixedpars=TRUE)
    expect_lt(deriv_error(msmtest5), err)
    msmtest5 <- msm(state ~ time, qmatrix = fiveq, subject = ptnum, data = bos, exacttimes=TRUE, qconstraint=c(1,2,1,2,1,2,1), fixedpars=TRUE)
    expect_lt(deriv_error(msmtest5), err)
    msmtest5 <- msm(state ~ time, qmatrix = fiveq, covariates = ~time, subject = ptnum, data = bos, exacttimes=TRUE, fixedpars=TRUE)
    expect_lt(deriv_error(msmtest5), err)
    msmtest5 <- msm(state ~ time, qmatrix = fiveq, covariates = ~time, constraint=list(time=c(1,2,1,2,1,2,2)),
                   subject = ptnum, data = bos, exacttimes=TRUE, fixedpars=TRUE)
    expect_lt(deriv_error(msmtest5), err)
})

if (0) {
### NOTE: NUMERIC DERIVS BREAK WITH THESE MATRICES when analyticp=TRUE: CLOSE TO REPEATED EIGENVALUES.
    psor.1.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
    psor.1.q <- rbind(c(0,0.1,0,0),c(0,0,0.10001,0),c(0,0,0,0.1001),c(0,0,0,0))
    diag(psor.1.q) <- -rowSums(psor.1.q)
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.1.q, fixedpars=TRUE)
    psor.msm$paramdata$deriv.test
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.1.q, fixedpars=TRUE, analyticp=FALSE)
    psor.msm$paramdata$deriv.test
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.1.q, qconstraint=c(1,1,2), fixedpars=TRUE)
    psor.msm$paramdata$deriv.test
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.1.q, qconstraint=c(1,1,2), fixedpars=TRUE, analyticp=FALSE)
    psor.msm$paramdata$deriv.test
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.1.q, covariates = ~ollwsdrt+hieffusn,
                    constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), fixedpars=FALSE, analyticp=TRUE)
    psor.msm$paramdata$deriv.test
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.1.q, covariates = ~ollwsdrt+hieffusn,
                    constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), fixedpars=FALSE, analyticp=FALSE)
    psor.msm$paramdata$deriv.test
}

context("analytic derivatives of likelihood in HMMs")

test.df <- data.frame(time=1:2, obs=c(1,1), x=c(1,2), y=c(3,4))
test_that("Categorical, 2 obs",{
    tm <- msm(obs ~ time, qmatrix=rbind(c(0,1,0),c(0,0,0),c(0,0,0)), ematrix=rbind(c(0.8,0.1,0.1),c(0.1,0.9,0),c(0,0,0)), data=test.df, fixedpars=TRUE)
    expect_lt(deriv_error(tm), err)
})

test_that("Categorical, lots of obs ",{
    nobs <- 100
    test.df <- data.frame(time=1:nobs, obs=sample(c(1,2),size=nobs,replace=TRUE), x=c(1,2), y=c(3,4))
    tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), ematrix=rbind(c(0.8,0.2),c(0.9,0.1)), data=test.df, fixedpars=TRUE)
    expect_lt(deriv_error(tm), err)
})

test_that("Categorical, a covariate",{
    tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), hmodel=list(hmmCat(c(0.8,0.2)),hmmCat(c(0.9,0.1))), hcovariates=list(~x,~1),  data=test.df, fixedpars=TRUE)
    expect_lt(deriv_error(tm), err)
})

test_that("Categorical, a covariate on more than one state",{
    tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), hmodel=list(hmmCat(c(0.8,0.2)),hmmCat(c(0.9,0.1))), hcovariates=list(~x,~x),  data=test.df, fixedpars=TRUE)
    expect_lt(deriv_error(tm), err)
})

test_that("Categorical, 4 potential obs",{
    tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), hmodel=list(hmmCat(c(0.8,0.1,0.05,0.05)),hmmCat(c(0.05,0.9,0.02,0.03))),  data=test.df, fixedpars=TRUE)
    expect_lt(deriv_error(tm), err)
})

test_that("Derivatives not supported with misclassification constraints",{
    expect_warning(tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), ematrix=rbind(c(0.8,0.2),c(0.9,0.1)), econstraint=c(1,1), data=test.df, fixedpars=TRUE), "Analytic derivatives not available")
    expect_warning(tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), hmodel=list(hmmCat(c(0.8,0.1,0.05,0.05)),hmmCat(c(0.05,0.9,0.02,0.03))),  data=test.df, hconstraint=list(p=c(1,1,2,3,4,5)), fixedpars=TRUE), "Analytic derivatives not available")
    expect_warning(tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), hmodel=list(hmmCat(c(0.8,0.1,0.05,0.05)),hmmCat(c(0.05,0.9,0.02,0.03))),  data=test.df, hcovariates=list(~x+y,~x+y), hconstraint=list(p=c(1,1,2,3,4,5),x=c(1,2,2,3,4,5),y=c(1,2,3,3,3,3)), fixedpars=TRUE), "Analytic derivatives not available")
    expect_warning(tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), hmodel=list(hmmCat(c(0.8,0.2)),hmmCat(c(0.9,0.1))), hcovariates=list(~x,~x),  hconstraint=list(x=c(1,1)), data=test.df, fixedpars=TRUE), "Analytic derivatives not available")
})

test_that("Derivatives with CAV misclassification model",{
    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav[1:200,], qmatrix = oneway4.q, ematrix=ematrix, misccovariates = ~dage + sex, covariates = ~ dage, covinits = list(dage=c(0.1,0.2,0.3,0.4,0.5)), misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016)), fixedpars=TRUE)
    expect_lt(deriv_error(misc.msm), err)

    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav[1:2,], qmatrix = oneway4.q, ematrix=ematrix, fixedpars=TRUE)
    expect_lt(deriv_error(misc.msm), err)

    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav[1:2,],
                    qmatrix = rbind(c(0,0.5,0),c(0,0,0.5),c(0,0,0)),
                    hmodel=list(hmmCat(c(0.9,0.1,0)), hmmCat(c(0.1,0.8,0.1)), hmmCat(c(0,0.1,0.9))),
                    fixedpars=TRUE)
    expect_lt(deriv_error(misc.msm), err)
})
          
## others in slow/test_fits_hmm.r

test_that("simple exponential",{
    nobs <- 3
    test.df <- data.frame(time=1:nobs, obs=c(rexp(nobs,c(sample(c(1,2),size=nobs,replace=TRUE)))))
    tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), hmodel=list(hmmExp(1.5),hmmExp(2)), data=test.df, fixedpars=TRUE)
    expect_lt(deriv_error(tm), err)
})

test_that("Information matrix",{
    nobs <- 1000
    set.seed(1)
    test.df <- data.frame(time=1:nobs, obs=sample(c(1,2),size=nobs,replace=TRUE))
    p1 <- 0.2; p2 <- 0.2;  pr1 <- c(1-p1, p1) # P obs(1,2) | true 1
    pr2 <- c(p2, 1-p2) # P obs(1,2) | true 2
    (tm <- msm(obs ~ time, qmatrix=rbind(c(0,1),c(0,0)), hmodel=list(hmmCat(pr1),hmmCat(pr2)), data=test.df, fixedpars=TRUE, hessian=TRUE))
    expect_equal(c(0.88475550512612, 0.202501688704573, -0.474183198550202), tm$paramdata$info[1:3], tol=1e-05)
    tm$paramdata$opt$hessian
})

set.seed(22061976)
nsubj <- 100; nobspt <- 6
sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length=nobspt),
                     x = rnorm(nsubj*nobspt), y = rnorm(nsubj*nobspt)* 5 + 20)
three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))

set.seed(22061976)
nsubj <- 100; nobspt <- 6
sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length=nobspt),
                     x = rnorm(nsubj*nobspt), y = rnorm(nsubj*nobspt)* 5 + 20)

test_that("poisson",{
    hmodel3 <- list(hmmPois(6), hmmPois(12), hmmIdent(999))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, hmodel = hmodel3)
    sim.hid <- msm(obs ~ time, subject=subject, data=sim2.df, qmatrix=three.q, hmodel=hmodel3, fixedpars=TRUE)
    expect_lt(deriv_error(sim.hid), err)
})

test_that("binomial",{
    hmodel3 <- list(hmmBinom(10, 0.1), hmmBinom(20, 0.3), hmmIdent(999))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, hmodel = hmodel3)
    sim.hid <- msm(obs ~ time, subject=subject, data=sim2.df, qmatrix=three.q, hmodel=hmodel3, fixedpars=TRUE)
    expect_lt(deriv_error(sim.hid), err)
})

test_that("negative binomial",{
    hmodel3 <- list(hmmNBinom(10, 0.1), hmmNBinom(20, 0.3), hmmIdent(999))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, hmodel = hmodel3)
    sim.hid <- msm(obs ~ time, subject=subject, data=sim2.df, qmatrix=three.q, hmodel=hmodel3, fixedpars=TRUE)
    expect_lt(deriv_error(sim.hid), err)
})

test_that("beta",{
### some kind of underflow with about 200 obs or more.   big derivs, prob poorly identified model
    hmodel3 <- list(hmmBeta(0.5,0.5), hmmBeta(2, 2), hmmIdent(999))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, hmodel = hmodel3)
    sim.hid <- msm(obs ~ time, subject=subject, data=sim2.df[1:100,], qmatrix=three.q, hmodel=hmodel3, fixedpars=TRUE)
    expect_lt(deriv_error(sim.hid), err)
})

test_that("t",{
    hmodel3 <- list(hmmT(1, 2, 2), hmmT(4, 2, 3), hmmIdent(999))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, hmodel = hmodel3)
    sim.hid <- msm(obs ~ time, subject=subject, data=sim2.df[1:100,], qmatrix=three.q, hmodel=hmodel3, fixedpars=TRUE)
    expect_lt(deriv_error(sim.hid), err)
})

options(msm.test.analytic.derivatives=NULL)
