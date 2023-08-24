psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  
                covariates = ~ollwsdrt+hieffusn, 
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), 
                control=list(fnscale=1))

## Merge states 2 and 3
psor$state3 <- ifelse(psor$state==3, 2,
                      ifelse(psor$state==4, 3, psor$state))
psor3.q <- psor.q[1:3,1:3]
psor3.msm <- msm(state3 ~ months, subject=ptnum, data=psor, qmatrix = psor3.q,  
                 covariates = ~ollwsdrt+hieffusn, 
                 constraint = list(hieffusn=c(1,1),ollwsdrt=c(1,1)), 
                 control=list(fnscale=1))

test_that("DRAIC",{
  d <- draic.msm(psor.msm, psor3.msm)
  expect_true(is.numeric(d$draic))
  dl <- draic.msm(psor.msm, psor3.msm, likelihood.only=TRUE)
  expect_equal(d$lik.restricted["complex","-LL"], dl[["complex"]])
  expect_equal(as.numeric(dl[["complex"]]), 415.032735368145, tolerance=1e-06)
  expect_true(is.numeric(draic.msm(psor.msm, psor3.msm, tl=0.99)$ti["Lower"]))
})

test_that("DRAIC with observed information", {
  skip_on_cran()
  expect_error(draic.msm(psor.msm, psor3.msm, information="observed"), NA)
})

test_that("DRLCV",{
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor[1:29,], qmatrix = psor.q)
    psor3.msm <- msm(state3 ~ months, subject=ptnum, data=psor[1:29,], qmatrix = psor3.q)
    expect_error(drlcv.msm(psor.msm, psor3.msm, verbose=FALSE), NA)
})

