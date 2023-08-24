context("phase-type models")

test_that("one phased state, model fit", {
    ## num overflow in lik with bfgs, even with fnscale. NM method does better. TODO does bad inits cause that?
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, 
                    phase.states=c(1), pci=5, method="Nelder-Mead", control=list(fnscale=1000))
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
