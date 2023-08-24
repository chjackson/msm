test_that("pearson.msm help example",{
  psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor,
                  qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                  constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
  p <- pearson.msm(psor.msm, timegroups=2, intervalgroups=2, covgroups=2)
  expect_equal(p$test[["p"]],0)
})

#test_that("pearson with bootstrap",{
#  psor.msm <- msm(state ~ months, subject=ptnum, data=psor[1:29,],
#                  qmatrix = psor.q)
#  p <- pearson.msm(psor.msm, timegroups=2, intervalgroups=2, covgroups=1, boot=TRUE, B=10)
#})