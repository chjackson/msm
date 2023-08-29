test_that("bootstrap iterations that returned an error are dropped",{
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor[1:100,], qmatrix = psor.q, 
                  covariates = ~ollwsdrt+hieffusn,fixedpars=FALSE, control=list(maxit=10))
  random_error <- function(x){if (rbinom(1,1,0.5)) stop("Error") else 1}
  set.seed(1)
  q.list <- boot.msm(psor.msm, random_error, B=10)
  expect_equal(length(q.list), 6)
})

