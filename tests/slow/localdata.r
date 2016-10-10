context("Tests with local datasets")

c2.df <- read.table("~/work/msm/tests/jeanluc/donneesaveccancerPT.txt", header=TRUE)
qx <- rbind( c(0, 0.005, 0, 0, 0), c(0, 0, 0.01, 0.02,0), c(0, 0, 0, 0.04, 0.03), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))

test_that("multiple death states",{
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df,
                  qmatrix=qx, death=c(4, 5), method="BFGS", fixedpars = TRUE)
    expect_equal(70646.727505836, c2.msm$minus2loglik, tol=1e-06)
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df,
                  qmatrix=qx, method="BFGS", fixedpars = TRUE)
    expect_equal(62915.1638036017, c2.msm$minus2loglik, tol=1e-06)
})

test_that("obstype vector with multiple deaths",{
    obstype <- ifelse(c2.df$state %in% c(4,5), 3, 1)
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, qmatrix=qx,
                  obstype=obstype, method="BFGS", fixedpars = TRUE)
    expect_equal(70646.727505836, c2.msm$minus2loglik, tol=1e-06)
    obstype <- rep(1, length(c2.df$state))
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, obstype=obstype,
                  qmatrix=qx, method="BFGS", fixedpars = TRUE)
    expect_equal(62915.1638036017, c2.msm$minus2loglik, tol=1e-06)
})

test_that("multiple death states with misclassification",{
  qx <- rbind( c(0, 0.005, 0, 0, 0), c(0, 0, 0.01, 0.02,0), c(0, 0, 0, 0.04, 0.03), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
  ex <- rbind( c(0, 0, 0, 0, 0), c(0, 0, 0.1, 0, 0), c(0, 0.1, 0, 0, 0), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0) )
  c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, qmatrix=qx, ematrix=ex, death=c(4, 5), method="BFGS", fixedpars = TRUE)
  expect_equal(70084.3665626129, c2.msm$minus2loglik, tol=1e-06)
  d45 <- rep(1, nrow(c2.df)); d45[c2.df$state %in% c(4,5)] <- 3
  c22.msm <- msm(state~years, subject=PTNUM, data=c2.df, qmatrix=qx, ematrix=ex, obstype=d45, method="BFGS", fixedpars = TRUE)
  expect_equal(c22.msm$minus2loglik, c2.msm$minus2loglik, tol=1e-06)
})

test_that("Diabetic retinopathy data from MARKOV",{
    marsh.df <- read.table("~/work/msm/tests/markov/test.dat", col.names=c("subject","eyes","time","duration","hba1"))
    marsh.df$hba1 <- marsh.df$hba1 - mean(marsh.df$hba1)
    marsh.msm <- msm(eyes ~ time, subject=subject, qmatrix = rbind(c(0,0.02039,0,0), c(0.007874,0,0.01012,0), c(0,0.01393,0,0.01045), c(0,0,0,0)),  covariates = ~ hba1, data = marsh.df, fixedpars=TRUE)
    expect_equal(335.897217906310, marsh.msm$minus2loglik, tol=1e-06)
    marsh.msm <-  msm(eyes ~ time, subject=subject, qmatrix = rbind(c(0,0.02039,0,0), c(0.007874,0,0.01012,0), c(0,0.01393,0,0.01045), c(0,0,0,0)), covariates = ~ hba1, data = marsh.df) # Nelder-Mead runs out of iterations, test fixed in 1.4
    expect_equal(310.100684750764, marsh.msm$minus2loglik, tol=1e-06)
    expect_equal(-0.0635250146505878, qmatrix.msm(marsh.msm, covariates=0)$estimates[1,1], tol=1e-06)
})
