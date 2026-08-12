test_that("pearson.msm help example",{
  psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor,
                  qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                  constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
  p <- pearson.msm(psor.msm, timegroups=2, intervalgroups=2, covgroups=2)
  expect_equal(p$test[["p"]],0)

  # with exact death times
  set.seed(1)
  cav.msm <- msm( state ~ years, subject=PTNUM, data = cav[1:1000,],
                  qmatrix = twoway4.q, deathexact = TRUE, fixedpars=FALSE)
  p <- pearson.msm(cav.msm)
  expect_equal(p$test[["stat"]],67,tol=1)
})

#test_that("pearson with bootstrap",{
#  psor.msm <- msm(state ~ months, subject=ptnum, data=psor[1:29,],
#                  qmatrix = psor.q)
#  p <- pearson.msm(psor.msm, timegroups=2, intervalgroups=2, covgroups=1, boot=TRUE, B=10)
                                        #})

test_that("pearson.msm, models with interactions",{
  skip_on_cran()
  expect_no_error({
    psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor,
                    qmatrix = psor.q, covariates = ~ollwsdrt:hieffusn)
    pearson.msm(psor.msm, timegroups=2, intervalgroups=2, covgroups=2)
  })
})

test_that("pearson.msm with next.obstime and exact death times",{

  cav.msm <- msm(state ~ years, subject=PTNUM, data=cav, qmatrix=twoway4.q, death=4)

  cav$lastobs <- rev(!duplicated(rev(cav$PTNUM)))
  cav$nexttime <- ifelse(cav$lastobs, NA, c(cav$years[2:nrow(cav)], NA))
  cav$next.scheduled.time <- ifelse(cav$state==4 | cav$lastobs, 
                                    cav$years + 1, cav$nexttime)
  expect_no_warning({
    gof <- pearson.msm(cav.msm, next.obstime = cav$next.scheduled.time, 
                     maxtimes = max(cav$years) + 1)
  })
})
