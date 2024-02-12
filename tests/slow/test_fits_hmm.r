context("Hidden Markov models: model fits and slow tests")

three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))
four.q <-  rbind(c(0, exp(-6), 0, exp(-9)), c(0, 0, exp(-6.01), exp(-9)), c(0, 0, 0, exp(-6.02)), c(0, 0, 0, 0))
five.q <-  rbind(c(0, exp(-6), 0, 0, exp(-9)),
                 c(0, 0, exp(-6.01), 0, exp(-9)),
                 c(0, 0, 0, exp(-6.02), exp(-6.03)),
                 c(0, 0, 0, 0, exp(-6.04)),
                 c(0, 0, 0, 0, 0))
hmodel1 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))
fev1.msm <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel1,
                hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                hconstraint = list(acute = c(1,1)), center=FALSE)

test_that("HMM normal model fit",{
    stopifnot(isTRUE(all.equal(51597.8909140275, fev1.msm$minus2loglik, tol=1e-06)))
    q <- qmatrix.msm(fev1.msm)
    stopifnot(isTRUE(all.equal(0.000801187055541291, q$estimates[2,3], tol=1e-05)))
})

test_that("Viterbi with normal HMMs",{
    keep <- fev$ptnum==1 & fev$fev<999
    vit <- viterbi.msm(fev1.msm)[keep,]
    max1 <- max(vit$time[vit$fitted==1])
    min2 <- min(vit$time[vit$fitted==2])
    stopifnot(isTRUE(all.equal(2377, max1)))
    stopifnot(isTRUE(all.equal(2406, min2))) # one obs different in 1.3, assume because of memory bug fixed
    if (interactive())  {
        plot(fev$days[keep], fev$fev[keep], type="l", ylab=expression(paste("% baseline ", FEV[1])), xlab="Days after transplant")
        abline(v = mean(max1,min2), lty=2)
        text(max1 - 500, 50, "STATE 1")
        text(min2 + 500, 50, "STATE 2")
    }
})

test_that("HMM measurement error outcomes",{ 
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
    (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3, fixedpars=TRUE))
    expect_equal(52567.734400762, fev3.hid$minus2loglik, tol=1e-06)
    (fev4.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=four.q, death=4, hmodel=hmodel4, fixedpars=TRUE))
    expect_equal(50202.2171575953, fev4.hid$minus2loglik, tol=1e-06)
    (fev5.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=five.q, death=5, hmodel=hmodel5, fixedpars=TRUE))
    expect_equal(49308.4942559547, fev5.hid$minus2loglik, tol=1e-06)

    (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3,
                     hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                     fixedpars=TRUE, center=FALSE))
    expect_equal(52219.7372318393, fev3.hid$minus2loglik, tol=1e-06)
    (fev4.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=four.q, death=4, hmodel=hmodel4,
                     hcovariates=list(~acute, ~acute, ~acute, NULL), hcovinits = list(-8, -8, -8, NULL),
                     fixedpars=TRUE, center=FALSE))
    expect_equal(49933.2168265887, fev4.hid$minus2loglik, tol=1e-06)
    (fev5.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=five.q, death=5, hmodel=hmodel5,
                     hcovariates=list(~acute, ~acute, ~acute, ~acute, NULL), hcovinits = list(-8, -8, -8, -8, NULL), fixedpars=TRUE, center=FALSE))
    expect_equal(49167.8668910928, fev5.hid$minus2loglik, tol=1e-06)
})



context("analytic derivatives in continuous HMMs")

test_that("normal",{
    hmodel3 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))
    (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3, fixedpars=TRUE))
})

test_that("log-normal",{
    hmodel3 <- list(hmmLNorm(mean=log(100), sd=1), hmmLNorm(mean=log(54), sd=1), hmmIdent(999))
    (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3, fixedpars=TRUE))
})

test_that("gamma",{
    hmodel3 <- list(hmmGamma(shape=100, rate=1), hmmGamma(shape=54, rate=1), hmmIdent(999))
    (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3, fixedpars=TRUE))
})

test_that("weibull",{
    hmodel3 <- list(hmmWeibull(shape=1, scale=100), hmmWeibull(shape=1, scale=54), hmmIdent(999))
    (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3, fixedpars=TRUE))
})

test_that("Constraints",{
    hmodel3 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))
    fev$x <- fev$acute+2
    fev$y <- fev$acute+10
    options(msm.test.analytic.derivatives=TRUE)
    (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev[1:100,], qmatrix=three.q, death=3, hmodel=hmodel3, hcovariates=list(~x+y, ~x+y, NULL), hconstraint=list(sd=c(1,1),x=c(1,1)), fixedpars=TRUE))
    expect_lt(deriv_error(fev3.hid), 1e-04)
    options(msm.test.analytic.derivatives=NULL)
})

test_that("HMMs with one observation for a subject", {
  fevno1 <- fev
  fevno1 <- fevno1[fevno1$ptnum != 1,]
  fev1 <- fev
  fev1 <- fev1[!(fev1$ptnum==1 & fev1$days>191), ]
  head(fev1)
  head(fevno1)
  fev1.msm <- msm(fev ~ days, subject=ptnum, data=fev1, qmatrix=three.q, 
                  death=3, hmodel=hmodel1, fixedpars=TRUE)
  fev1.msm
  
  fevno1.msm <- msm(fev ~ days, subject=ptnum, data=fevno1, qmatrix=three.q, 
                    death=3, hmodel=hmodel1, fixedpars=TRUE)
  
  ## Check gives same result as inserting a uninformative censored state 
  ## at the second obs 
  censrow <- data.frame(ptnum=1, days=200, fev=-99, acute=0)
  fevcens <- rbind(fev1[1,], censrow, fev1[-1,])
  fevcens$obstrue <- NA
  fevcens$obstrue[2] <- 1 # watch this syntax for censor+hmm! See help(msm)
  fevcens.msm <- msm(fev ~ days, subject=ptnum, data=fevcens, qmatrix=three.q, 
                     death=3, hmodel=hmodel1, censor=-99, censor.states = 1:3,
                     obstrue = obstrue, fixedpars=TRUE)
  fevcens.msm

  expect_equal(logLik.msm(fev1.msm), logLik.msm(fevcens.msm))
  
  expect_true(logLik.msm(fevno1.msm) != logLik.msm(fevcens.msm))
})
