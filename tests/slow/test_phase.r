context("phase-type models")

test_that("one phased state", {
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, fixedpars=TRUE,
                    phase.states=c(1))
    expect_equal(psor.msm$minus2loglik, 1308.35590680014, tol=1e-06)
}
)
test_that("one phased state, model fit", {
    ## num overflow in lik with bfgs, even with fnscale. NM method does better. TODO does bad inits cause that?
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, 
                    phase.states=c(1), method="Nelder-Mead")
    expect_equal(psor.msm$minus2loglik, 1234.568505251202, tol=1e-06)
    psor0.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q)
    ## better fit compared to unphased
    AIC(psor.msm, psor0.msm)

    ## viterbi
    vit <- viterbi.msm(psor.msm)
    ms <- vit$observed %in% psor.msm$qmodel$markov.states
    expect_equal(vit$observed[ms], as.numeric(as.character(vit$fitted[ms])))
    ps <- vit$observed %in% psor.msm$qmodel$phase.states
    expect_equal(grep("P", vit$fitted[ps]), seq_along(vit$fitted[ps]))

    ## outputs are for aggregate state space 
    pmatrix.msm(psor.msm, ci="normal", B=100)
    prevalence.msm(psor.msm)

    ## phase means
    phasemeans.msm(psor.msm, ci="normal")
#    phasemeans.msm(psor.msm, ci="boot", B=4, cores=4)
})

test_that("two phased states", {
    ## note phase transition within state 2 appears hard to estimate.
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, fixedpars=TRUE,
                    phase.states=c(1,2))
    expect_equal(psor.msm$minus2loglik, 1303.768173132893, tol=1e-06)
})

test_that("supplying initial values",{
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1,2),
                     phase.inits = list(
                     list(trans=c(0.42032), exit=c(0.17404,0.04809)),
                     list(trans=c(0.42032), exit=c(0.17404,0.04809))),                 
                     fixedpars=TRUE)
    expect_equal(psor.msm$minus2loglik, 1280.788860486683, tol=1e-06)
})

test_that("errors in initial values",{
    expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1,2), phase.inits = "foo"), "phase.inits should be a list")
    expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1,2), phase.inits = list(1,2,3)), "phase.inits of length 3, but there are 2 phased states")
    expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1,2), phase.inits = list(1,2)), "phase.inits.+ list of length 1, should be 2")
    expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1), phase.inits = list(list(trans=c(0.42032, 1), exit=c(0.17404,0.04809)))), "phase.inits.+trans of length 2, should be 1")
    expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1), phase.inits = list(list(trans=c(0.42032), exit=c(0.17404,0.04809,1,1)))), "phase.inits.+exit has 4 columns, should be 2")
    expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1), phase.inits = list(list(trans=c(0.42032), exit=matrix(c(0.17404,0.04809,1,1),ncol=2,byrow=TRUE)))), "phase.inits.+exit has 2 rows, but there are 1 exit states from this state")
})

test_that("with covariates", {
## Names for transition-specific covs are interpreted as state numbers in new state space
    psorc.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, # fixedpars=TRUE,
                     covariates = list("2-3" = ~ ollwsdrt),
                     phase.states=c(1), method="Nelder-Mead", control=list(fnscale=1000,maxit=10000))
    psorc.msm
})

test_that("with qconstraint", {
    psorc.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, # fixedpars=TRUE,
                     qconstraint = c(1,1,2,3,4),
                     phase.states=c(1), method="Nelder-Mead", control=list(fnscale=1000,maxit=10000))
    psorc.msm
})

test_that("with censoring", {
    psor2 <- psor
    psor2$state[c(16,18)] <- 99
    psorc.msm <- msm(state ~ months, subject=ptnum, data=psor2, qmatrix = psor.q, fixedpars=TRUE,
                     censor = 99, censor.states=c(3,4))
    expect_equal(psorc.msm$minus2loglik, 1290.189703452163, tol=1e-06)
})

test_that("HMM on top", {
### Note 1.6 clarifies obstrue not supported for hmmCat: so test result in 1.5 wrong
    miscnew.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                       qmatrix = oneway4.q, death = 5, #obstrue=firstobs,
                       fixedpars=TRUE,
                       # could try to fit, but weakly identifiable.  fit seems improved though (-2LL 3864 vs 3951)
                       # control=list(trace=1,REPORT=1,fnscale=4000, maxit=10000),
                       phase.states = 1, 
                       hmodel=list(
                       hmmCat(prob=c(0.9, 0.1, 0, 0)),
                       hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                       hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent(4))
                       )
    expect_equal(miscnew.msm$minus2loglik, 4357.70908245511, tol=1e-06)  
})
