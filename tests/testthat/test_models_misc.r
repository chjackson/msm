context("msm misclassification model likelihoods")

test_that("cav misclassification model with no covariates",{
    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=ematrix, deathexact = 4, fixedpars=TRUE)
    expect_equal(4296.9155995778, misc.msm$minus2loglik, tol=1e-06)
    miscnew.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                       qmatrix = oneway4.q, deathexact = 4, fixedpars=TRUE,
                       hmodel=list(
                       hmmCat(prob=c(0.9, 0.1, 0, 0)),
                       hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                       hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent())
                       )
    expect_equal(miscnew.msm$minus2loglik, misc.msm$minus2loglik)  
})

test_that("cav misclassification model with covariates on transition rates",{
    misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=ematrix, deathexact = 4, fixedpars = TRUE, covariates = ~ sex, covinits=list(sex=rep(0.1, 5)))
    expect_equal(4299.38058878142, misccov.msm$minus2loglik, tol=1e-06)
})

test_that("cav misclassification model with covariates on misclassification probabilities",{
    misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                       qmatrix = oneway4.q, ematrix=ematrix, deathexact = 4, fixedpars=TRUE,
                       misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016)))
    expect_equal(4306.3077053482, misccov.msm$minus2loglik, tol=1e-06)
    misccovnew.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                          qmatrix = oneway4.q, deathexact = 4, fixedpars=TRUE, center=TRUE,
                          hmodel=list(
                          hmmCat(prob=c(0.9, 0.1, 0, 0)),
                          hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                          hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()),
                          hcovariates=list(~dage + sex, ~dage + sex, ~dage + sex, ~1),
                          hcovinits = list(c(0.01,-0.013), c(0.02,-0.014,0.03,-0.015), c(0.04,-0.016), NULL)
                          )
    expect_equal(misccov.msm$minus2loglik, misccovnew.msm$minus2loglik, tol=1e-06)
})

test_that("misclassification model with no misclassification reduces to simple",{
    nomisc.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = twoway4.q, ematrix=matrix(0, nrow=4, ncol=4), deathexact = 4, fixedpars=TRUE)
    simple.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = twoway4.q, deathexact = 4, fixedpars=TRUE)
    expect_equal(nomisc.msm$minus2loglik, simple.msm$minus2loglik)
})

test_that("misclassification model with obstrue",{
    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=ematrix, deathexact = 4, fixedpars=TRUE, obstrue=firstobs)
    expect_equal(4165.84711809003, misc.msm$minus2loglik, tol=1e-06)
})

test_that("misclassification model with exact times",{
    miscexact.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=ematrix, exacttimes=TRUE, fixedpars=TRUE)
    expect_equal(4864.14764195147, miscexact.msm$minus2loglik, tol=1e-06)
})

test_that("misclassification model with initprobs",{
    miscinits.msm <- msm(state ~ years, subject = PTNUM, data = cav,  qmatrix = oneway4.q, ematrix=ematrix, deathexact = 4, initprobs=c(0.7, 0.1, 0.1, 0.1), fixedpars=TRUE)
    expect_equal(4725.9078185031, miscinits.msm$minus2loglik, tol=1e-06)
})

test_that("misclassification model with censoring",{
    misccens.msm <- msm(state ~ years, subject = PTNUM, data = cav.cens,
                    qmatrix = oneway4.q, ematrix=ematrix, deathexact=TRUE, censor=99, fixedpars=TRUE)
    expect_equal(4025.42265024404, misccens.msm$minus2loglik, tol=1e-06)
})

test_that("misclassification model with two types of censoring",{
    misccens.msm <- msm(state ~ years, subject=PTNUM, data=cav.cens2, qmatrix=oneway4.q, ematrix=ematrix, censor=c(99, 999), deathexact=4, censor.states=list(c(1,2,3), c(2,3)), fixedpars=TRUE)
    expect_equal(3811.69640533587, misccens.msm$minus2loglik, tol=1e-06)
    cav.cens2$obstrue <- as.numeric(cav.cens2$state %in% c(999))
    misccens.msm <- msm(state ~ years, subject=PTNUM, data=cav.cens2, qmatrix=oneway4.q, ematrix=ematrix, censor=c(99, 999), deathexact=4, obstrue=obstrue, censor.states=list(c(1,2,3), c(2,3)), fixedpars=TRUE)
    expect_equal(3822.04540210944, misccens.msm$minus2loglik, tol=1e-06)
})

test_that("misclassification model with no misclassification reduces to simple, with censoring",{
    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav.cens,
                    qmatrix = twoway4.q, ematrix=matrix(0, nrow=4, ncol=4), censor=99, deathexact=TRUE, fixedpars=TRUE)
    simple.msm <- msm(state ~ years, subject = PTNUM, data = cav.cens, qmatrix = twoway4.q, deathexact=TRUE, censor=99, fixedpars=TRUE)
    expect_equal(misc.msm$minus2loglik, simple.msm$minus2loglik, tol=1e-06)
    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav.cens,
                    qmatrix = twoway4.q, ematrix=matrix(0, nrow=4, ncol=4), censor=99, fixedpars=TRUE)
    simple.msm <- msm(state ~ years, subject = PTNUM, data = cav.cens, qmatrix = twoway4.q, censor=99, fixedpars=TRUE)
    expect_equal(misc.msm$minus2loglik, simple.msm$minus2loglik, tol=1e-06)    
})

test_that("misclassification model with no misclassification reduces to simple, using hmmCat",{
    miscnew.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = twoway4.q,
                       hmodel=list(hmmCat(prob=c(1, 0, 0, 0)), hmmCat(prob=c(0, 1, 0, 0)), hmmCat(prob=c(0, 0, 1, 0)), hmmIdent()), fixedpars=TRUE)
    simple.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = twoway4.q, fixedpars=TRUE)
    expect_equal(miscnew.msm$minus2loglik, simple.msm$minus2loglik, tol=1e-06)    
})

test_that("can't mix ematrix and hcovariates",{
    expect_error( misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                                     qmatrix = oneway4.q, ematrix=ematrix, deathexact = 4, fixedpars=1:17,
                                     hcovariates=list(~dage + sex, ~dage + sex, ~dage + sex, ~1),
                                     hcovinits = list(c(0.01,0.013), c(0.01,0.013,0.01,0.013), c(0.01,0.013), NULL) ), "hcovariates have been specified, but no hmodel")
})

test_that("data inconsistent with initprobs/ematrix",{
    cav2 <- cav
    cav2$state[c(1,8)] <- 3
    expect_warning(msm(state ~ years, subject = PTNUM, data = cav2, qmatrix = oneway4.q, ematrix=ematrix, deathexact = 4, fixedpars=TRUE), "First observation .+ is impossible")
})
          
test_that("various errors",{
    wrong.e <- "foo"
    expect_error(misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=wrong.e, deathexact = 4, fixedpars=TRUE),"ematrix should be a numeric matrix")
    wrong.e <- 1
    expect_error(misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=wrong.e, deathexact = 4, fixedpars=TRUE),"ematrix should be a numeric matrix")
    wrong.e <- cbind(c(0,1,2), c(0,1,2))
    expect_error(misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=wrong.e, deathexact = 4, fixedpars=TRUE),"Number of rows and columns of ematrix should be equal")
    wrong.e <- cbind(c(0,1), c(0,2))
    expect_error(misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=wrong.e, deathexact = 4, fixedpars=TRUE),"Dimensions of qmatrix and ematrix should be the same")
    expect_warning(msm(state ~ years, subject = PTNUM, data = cav, qmatrix = oneway4.q, ematrix=ematrix, deathexact = 4, fixedpars=TRUE, gen.inits=TRUE), "gen.inits not supported for hidden Markov models, ignoring")
})
