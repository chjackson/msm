test_that("simfitted.msm",{ 
  skip_on_cran()
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt)
  expect_equal(simfitted.msm(psor.msm)$time[1], psor$months[1])
  
  psorc <- psor
  psorc$state[2] <- 99 
  psorc.msm <- msm(state ~ months, subject=ptnum, data=psorc, qmatrix = psor.q, 
                   covariates = ~ollwsdrt, censor=99)
  set.seed(1)
  simdat <- simfitted.msm(psorc.msm)$state[2]
  expect_equal(simdat, psorc$state[2])
})
