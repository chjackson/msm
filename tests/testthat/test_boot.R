test_that("bootstrap iterations that returned an error are dropped",{
  suppressWarnings(
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor[1:100,], 
                    qmatrix = rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0)), 
                    covariates = ~ollwsdrt+hieffusn,fixedpars=FALSE, control=list(maxit=10))
  )
  random_error <- function(x){if (rbinom(1,1,0.5)) stop("Error") else 1}
  set.seed(1)
  suppressWarnings(
    q.list <- boot.msm(psor.msm, random_error, B=10)
  )
  expect_lt(length(q.list), 10)
})

