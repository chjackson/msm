psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  
                covariates = ~ollwsdrt)

test_that("ppass.msm",{
  pp <- ppass.msm(psor.msm, tot=10)
  pm <- pmatrix.msm(psor.msm, t=10)
  expect_equal(pp[,4], pm[,4]) # state 4 is absorbing
  pp <- ppass.msm(qmatrix=twoway4.q, tot=1000)
  expect_equal(pp[1,2], 0.5)
  expect_warning(ppass.msm(qmatrix=twoway4.q, tot=100, ci="normal"), 
                 "No fitted model supplied: not calculating confidence intervals")
})

psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  
                covariates = ~ollwsdrt, pci=5)

test_that("ppass.msm with pci",{
  ## if tot is before 5, then do by hand
  qq <- qmatrix.msm(psor.msm, covariates=0, ci="none")
  qq[3,] <- 0
  expect_equal(ppass.msm(psor.msm, tot=3, covariates=0)[,3],
               MatrixExp(qq, t=3)[,3])
               
  ## if tot is 7, then do one cycle of piecewise p 
  qq0 <- qmatrix.msm(psor.msm, covariates=0, ci="none")
  qq1 <- qmatrix.msm(psor.msm, covariates=list(ollwsdrt=0,timeperiod="[5,Inf)"), ci="none")
  qq0[3,] <- qq1[3,] <- 0
  expect_equal(
    pmatrix.piecewise.msm(qlist = list(qq0,qq1), t1=0, t2=7, times=5)[,3],
    ppass.msm(psor.msm, tot=7, covariates=0, ci="none")[,3])

  pp <- ppass.msm(psor.msm, tot=7, ci="normal", B=10)
  expect_gt(pp[["U"]][1,3], pp[["estimates"]][1,3] )
})

psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  
                covariates = ~ollwsdrt + months)

test_that("ppass.msm with non-pci time-dependent models",{
  qq0 <- qmatrix.msm(psor.msm, covariates=0, ci="none")
  qq1 <- qmatrix.msm(psor.msm, covariates=list(ollwsdrt=0,months=3), ci="none")
  qq0[3,] <- qq1[3,] <- 0
  pp <- ppass.msm(psor.msm, tot=7, piecewise.times = 3,
                  piecewise.covariates = list(list(ollwsdrt=0, months=0),
                                              list(ollwsdrt=0, months=3)))
  expect_equal(pp[,3], 
               pmatrix.piecewise.msm(qlist = list(qq0,qq1), t1=0, t2=7, times=3)[,3],
  )

  expect_error(ppass.msm(psor.msm, tot=7, piecewise.times = 3,
                         piecewise.covariates = list(list(ollwsdrt=0, months=3))),
               "Number of covariate lists must be one greater")
  expect_error(ppass.msm(psor.msm, tot=7, piecewise.times = c(3,2),
                         piecewise.covariates = list(list(ollwsdrt=0, months=3))),
               "times should be a vector of numbers in increasing order")
})
