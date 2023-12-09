three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))

test_that("hmodel2list",{
  hmodel3 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))
  (fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, 
                   hmodel=hmodel3, fixedpars=TRUE))

  (fev3.hid2 <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, deathexact=3, 
                    hmodel=hmodel2list(fev3.hid$hmodel), 
                    fixedpars=TRUE))
  expect_equal(fev3.hid2$minus2loglik, fev3.hid$minus2loglik)
  
  hm <- hmodel2list(fev3.hid$hmodel)
  hml <- hmodel2list(fev3.hid$hmodel, hmmdist = FALSE)
  expect_equal(hm[[1]], do.call(hmmNorm, hml[[1]]))
})

test_that("hmodel2list with hmmCat",{
  misccovnew.msm <- msm(state ~ years, subject = PTNUM, data = cav[1:1000,],  
                        qmatrix = oneway4.q, death = 4,
                        fixedpars=TRUE,
                        hmodel=list(
                          hmmCat(prob=c(0.9, 0.1, 0, 0)),
                          hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                          hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()),
                        hcovariates=list(~dage + sex, ~dage + sex, ~dage + sex, ~1))
  
  misccovnew2.msm <- msm(state ~ years, subject = PTNUM, data = cav[1:1000,],  
                         qmatrix = oneway4.q, death = 4,
                         fixedpars=TRUE,
                         hmodel = hmodel2list(misccovnew.msm$hmodel),
                         hcovariates=list(~dage + sex, ~dage + sex, ~dage + sex, ~1))
  expect_equal(misccovnew2.msm$minus2loglik, misccovnew.msm$minus2loglik)
  hm <- hmodel2list(misccovnew.msm$hmodel)
  hml <- hmodel2list(misccovnew.msm$hmodel, hmmdist = FALSE)
  expect_equal(hm[[1]], do.call(hmmCat, hml[[1]]))
  expect_equal(hm[[2]], do.call(hmmCat, hml[[2]]))
})

## For univariate, each element has one element per par 
## For multivariate, each element has one element per outcome dist, which in turn has
## one element per par

test_that("hmodel2list with hmmMV",{
  fev$fevmat <- cbind(fev$fev, fev$fev*runif(nrow(fev), 0.95, 1.05))
  hmodel3mv <- list(hmmMV(hmmNorm(mean=100, sd=16), 
                          hmmNorm(mean=95, sd=18)), 
                    hmmMV(hmmNorm(mean=54, sd=16), 
                          hmmNorm(mean=50, sd=18)), 
                    hmmIdent(999))
  (fev3.hid <- msm(fevmat ~ days, subject=ptnum, 
                   data=fev, qmatrix=three.q, deathexact=3, 
                   hmodel=hmodel3mv, fixedpars=TRUE))
  
  (fev3.hid2 <- msm(fevmat ~ days, subject=ptnum, 
                   data=fev, qmatrix=three.q, deathexact=3, 
                   hmodel=hmodel2list(fev3.hid$hmodel), fixedpars=TRUE))
  expect_equal(fev3.hid$minus2loglik, fev3.hid2$minus2loglik)
  
  hm <- hmodel2list(fev3.hid$hmodel)
  hml <- hmodel2list(fev3.hid$hmodel, hmmdist=FALSE)
  expect_equal(hm[[1]], hmmMV(do.call(hmmNorm, hml[[1]][[1]]),
                              do.call(hmmNorm, hml[[1]][[2]])))
  expect_equal(hm[[2]], hmmMV(do.call(hmmNorm, hml[[2]][[1]]),
                              do.call(hmmNorm, hml[[2]][[2]])))
})
