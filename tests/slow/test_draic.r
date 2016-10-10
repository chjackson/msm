psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), control=list(fnscale=1))

## Merge states 2 and 3
psor$state3 <- ifelse(psor$state==3, 2,
               ifelse(psor$state==4, 3, psor$state))
psor3.q <- psor.q[1:3,1:3]
psor3.msm <- msm(state3 ~ months, subject=ptnum, data=psor, qmatrix = psor3.q,  covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1),ollwsdrt=c(1,1)), control=list(fnscale=1))

test_that("DRAIC",{
    expect_equal(draic.msm(psor.msm, psor3.msm)$draic, 0.02885907, tol=1e-06)
    expect_equal(as.numeric(draic.msm(psor.msm, psor3.msm, likelihood.only=TRUE)["complex"]), 415.032735368145, tolerance=1e-06)
    expect_equal(as.numeric(draic.msm(psor.msm, psor3.msm, tl=0.99)$ti["Lower"]), -0.0281113467856842, tol=1e-06)
})

if (0){
test_that("DRLCV",{
if (interactive()){
    draic.msm(psor.msm, psor3.msm, information="observed")
    expect_equal(drlcv.msm(psor.msm, psor3.msm)$drlcv, 0.03716199, tolerance=1e-06)
}
})
}
