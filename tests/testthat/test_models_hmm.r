context("Hidden Markov model likelihoods")

three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))

four.q <-  rbind(c(0, exp(-6), 0, exp(-9)), c(0, 0, exp(-6.01), exp(-9)), c(0, 0, 0, exp(-6.02)), c(0, 0, 0, 0))

five.q <-  rbind(c(0, exp(-6), 0, 0, exp(-9)),
                 c(0, 0, exp(-6.01), 0, exp(-9)),
                 c(0, 0, 0, exp(-6.02), exp(-6.03)),
                 c(0, 0, 0, 0, exp(-6.04)),
                 c(0, 0, 0, 0, 0))

hmodel3 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))
hmodel4 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=72.5, sd=10), hmmNorm(mean=42.5, sd=18), hmmIdent(999))
hmodel5 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=72.5, sd=10), hmmNorm(mean=62.5, sd=10), hmmNorm(mean=42.5, sd=18), hmmIdent(999))


test_that("HMM normal likelihoods: FEV data",{
    (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3, fixedpars=TRUE))
    expect_equal(52388.7381942858, fev3.hid$minus2loglik, tol=1e-06)
    (fev4.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=four.q, deathexact=4, hmodel=hmodel4, fixedpars=TRUE))
    expect_equal(50223.9497625937, fev4.hid$minus2loglik, tol=1e-06)
    (fev5.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=five.q, deathexact=5, hmodel=hmodel5, fixedpars=TRUE))
    expect_equal(49937.9840668066, fev5.hid$minus2loglik, tol=1e-06)
})

test_that("HMM with obstrue",{
    fev$obstrue <- as.numeric(fev$fev==999)
    fev$fev2 <- fev$fev; fev$fev2[fev$obstrue==1] <- 3
    hmodel32 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(3))
    (fev3.hid <- msm(fev2 ~ days, subject=ptnum, obstrue=obstrue, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel32, fixedpars=TRUE))
    expect_equal(52388.7381942858, fev3.hid$minus2loglik, tol=1e-06)

    fev$obstrue <- "foo"
    expect_error(fev3.hid <- msm(fev2 ~ days, subject=ptnum, obstrue=obstrue, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel32, fixedpars=TRUE), "obstrue should be logical or numeric")
    fev$obstrue <- as.numeric(fev$fev==999); fev$obstrue[1] <- 10
    expect_error(fev3.hid <- msm(fev2 ~ days, subject=ptnum, obstrue=obstrue, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel32, fixedpars=TRUE), "Interpreting \"obstrue\" as containing true states, but it contains values not in 0,1,...,3")

    ## can supply true state in obstrue
    fev$obstrue <- as.numeric(fev$fev==999); fev$obstrue[fev$obstrue==1] <- 3
    fev3.hid <- msm(fev2 ~ days, subject=ptnum, obstrue=obstrue, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel32, fixedpars=TRUE)
    expect_equal(52388.7381942858, fev3.hid$minus2loglik, tol=1e-06)
})

test_that("HMM normal likelihoods: FEV data: covariate on outcome",{
    (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3,
                     hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                     hconstraint = list(acute = c(1,1)),
                     fixedpars=TRUE, center=FALSE))
    expect_equal(52134.2372359988, fev3.hid$minus2loglik, tol=1e-06)
    (fev4.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=four.q, deathexact=4, hmodel=hmodel4,
                     hcovariates=list(~acute, ~acute, ~acute, NULL), hcovinits = list(-8, -8, -8, NULL),
                     hconstraint = list(acute = c(1,1,1)),
                     fixedpars=TRUE, center=FALSE))
    expect_equal(50095.8606697016, fev4.hid$minus2loglik, tol=1e-06)
    (fev5.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=five.q, deathexact=5, hmodel=hmodel5,
                     hcovariates=list(~acute, ~acute, ~acute, ~acute, NULL), hcovinits = list(-8, -8, -8, -8, NULL),
                     hconstraint = list(acute = c(1,1,1,1)),
                     fixedpars=TRUE, center=FALSE))
    expect_equal(49839.1627881087, fev5.hid$minus2loglik, tol=1e-06)
})

context("Hidden Markov model error handling")

test_that("errors in hmodel",{
    ## hmodel with wrong named parameters
    expect_error(
        hmodel3 <- list(hmmMETNorm(mean=100, sd=16, splat=8, lower=80, upper=Inf, meanerr=0),
                        hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                        hmmIdent(999)), "unused argument \\(splat"
        )
    ## hmodel with extra unnamed parameter
    expect_error(
        hmodel3 <- list(hmmMETNorm(mean=100, sd=16, sderr=8, lower=80, upper=Inf, meanerr=0),
                        hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                        hmmIdent(999, 3)), "unused argument \\(3"
        )
    ## hmodel with parameters unnamed, but correct number - OK.
    hmodel.OK <- list(hmmMETNorm(100, 16, 8, 80, Inf, 0),
                      hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                      hmmIdent(999))
    ## hmodel with too few parameters (named or unnamed)
    expect_error(
        hmodel3 <- list(hmmMETNorm(100, 16, 80, Inf),
                        hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                        hmmIdent(999)), "Parameter sderr for hmmMETNorm not supplied"
        )
### Initial values for certain parameters in wrong ranges.
    hmodel3 <- list(hmmMETNorm(mean=100, sd=-16, sderr=8, lower=80, upper=Inf, meanerr=0),
                    hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                    hmmIdent(999))
    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3, fixedpars=1:9), "Initial value -16 of parameter \"sd\" outside allowed range")

    expect_error(
        hmodel3 <- list(hmmMETNorm(mean=100, sd=16, sderr="splat", lower=80, upper=Inf, meanerr=0),
                        hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                        hmmIdent(999)), "Expected numeric values"
        )
    expect_error(
        hmodel3 <- list(
            hmmBinom(size=0.1, prob=0.5),
            hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
            hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()),
        "Value of size should be integer")
})


test_that("error handling in HMM fits",{
    hmodel3 <- list(hmmMETNorm(mean=100, sd=16, sderr=8, lower=80, upper=Inf, meanerr=0),
                    hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                    hmmIdent(999))
    hmodel4 <- list(hmmMETNorm(mean=100, sd=16, sderr=8, lower=80, upper=Inf, meanerr=0),
                    hmmMEUnif(sderr=8, lower=65, upper=80, meanerr=0),
                    hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=65, meanerr=0),
                    hmmIdent(999))
    hmodel5 <- list(hmmMETNorm(mean=100, sd=16, sderr=8, lower=80, upper=Inf, meanerr=0),
                    hmmMEUnif(sderr=8, lower=65, upper=80, meanerr=0),
                    hmmMEUnif(sderr=8, lower=50, upper=65, meanerr=0),
                    hmmMETNorm(mean=42, sd=18, sderr=8, lower=0, upper=50, meanerr=0),
                    hmmIdent(999))

### Wrong number of states in the hmodel, versus the qmatrix.
    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel4, fixedpars=1:9), "hmodel of length 4")
### Rubbish in hmodel list
    hmodel.rubbish <- list("should be a ", " list of ", "hmodel objects")
    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel.rubbish, fixedpars=1:9), "hmodel should be a list of HMM distribution objects")
    ## Death state wrong (not HMM-specific error, but no harm putting it in this file anyway)
    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=10, hmodel=hmodel3, fixedpars=1:9), "Death states indicator contains states not in 1, 2, 3")

    ## Covariate list of wrong length
    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3, hcovariates=list(~acute, ~acute), hcovinits = list(-8, -8, NULL),), "hcovariates of length 2, expected 3")

    ## Covariate list has rubbish in it
    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3, hcovariates=list("should be formulae", "rubbish", NULL), hcovinits = list(-8, -8, NULL)), "hcovariates should be a list of formulae or NULLs")

    ## Covariates not in the data
    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3, hcovariates=list(~acute, ~nonexistent, NULL), hcovinits = list(-8, -8, NULL)), "object .+ not found")

    ## Covariate inits of wrong length (just warning, ignores)
    expect_warning(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3, hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, -8), fixedpars=TRUE), "Initial values for hidden covariate effects do not match numbers of covariates")

    ## Rubbish in hcovinits (just warning, ignores)
    expect_warning(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3,hcovariates=list(~acute, ~acute, NULL), hcovinits = list("fooo", -8, NULL), fixedpars=TRUE), "hcovinits should be numeric")

    ## hconstraint has unknown parameters
    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3, hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL), hconstraint = list(nonexistent = c(1,1), acute = c(1,1))), "parameter .+ in hconstraint unknown")

    ## hconstraint has rubbish in (note character vectors are allowed)
    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3, hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL), hconstraint = list("rubbish", acute = c(1,1))), "parameter .+ in hconstraint unknown")

    expect_error(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, hmodel=hmodel3, hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL), hconstraint = list(sderr=c("wrong length"), acute = c(1,1))), "constraint for \"sderr\" of length 1, should be 2")

})
