## Make sure everything runs at least
## hard to check e.g. coverage
## e.g. width should be greater than asymp SE from Hessian. 

set.seed(1)
psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, 
                covariates = ~ollwsdrt+hieffusn,
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), fixedpars=FALSE)

test_that("bootstrap CIs run for all functions that support them",{
    skip_on_cran()
    q.list <- boot.msm(psor.msm, function(x)x$Qmatrices$baseline, B=2)
    expect_error(qmatrix.msm(psor.msm, ci="boot", B=2), NA)
    expect_error(pmatrix.msm(psor.msm, ci="boot", B=2), NA)
    totlos.msm(psor.msm, ci="boot", B=2, tot=1000)
    qratio.msm(psor.msm, c(1,2), c(2,3), ci="boot", B=2)
    pnext.msm(psor.msm, ci="boot", B=2)
    efpt.msm(psor.msm, tostate=3, ci="boot", B=2) 
    ppass.msm(psor.msm, tot=10, ci="boot", B=2)
    
    cavcens.msm <- msm(state ~ years, subject=PTNUM, data=cav.cens[500:1000,], 
                       qmatrix=twoway4.q, censor=99, fixedpars=c(2,5,7))
    pmatrix.msm(cavcens.msm, ci="boot", B=2)
    pci.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, 
                   covariates = ~ollwsdrt+hieffusn, 
                   constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), pci=5, fixedpars=FALSE)
    expect_error(pmatrix.msm(pci.msm, ci="boot", B=2), NA) 
})

test_that("normal CIs run for all functions that support them",{
    qmatrix.msm(psor.msm, ci="normal", B=10)
    pmatrix.msm(psor.msm, ci="normal", B=10)
    totlos.msm(psor.msm, ci="normal", B=10, tot=1000)
    qratio.msm(psor.msm, c(1,2), c(2,3), ci="normal", B=10)
    pnext.msm(psor.msm, ci="normal", B=10)
    efpt.msm(psor.msm, tostate=4, ci="normal", B=10)
    expect_error(ppass.msm(psor.msm, tot=10, ci="normal", B=10) , NA)
})


psor2 <- psor
psor2$ollwsdrt <- sample(c(0,1,2), size=nrow(psor), replace=TRUE)
psor2$ollwsdrt <- factor(ifelse(psor2$ollwsdrt==0, "foo", ifelse(psor2$ollwsdrt==1, "bar", "boing")))

test_that("bootstrap CIs with factor covariates",{
    skip_on_cran()
    psor2.msm <- msm(state ~ months, subject=ptnum, data=psor2, 
                     qmatrix = psor.q, covariates = ~ollwsdrt + hieffusn,
                     constraint = list(hieffusn=c(1,1,1)))
    expect_true(is.list(pmatrix.msm(psor2.msm, ci="boot", B=3)))
    expect_true(is.list(pmatrix.msm(psor2.msm, ci="boot", B=2,
                                    covariates=list(hieffusn=0, ollwsdrt="foo"))))
})

psor2$ptnum <- factor(psor2$ptnum)

test_that("bootstrap CIs with factor subject IDs and factor covariates",{
    skip_on_cran()
    psor2.msm <- msm(state ~ months, subject=ptnum, data=psor2, 
                     
                     qmatrix = psor.q, covariates = ~ollwsdrt + hieffusn,
                     constraint = list(hieffusn=c(1,1,1)))
    expect_true(is.list(pmatrix.msm(psor2.msm, ci="boot", B=3)))
    expect_true(is.list(pmatrix.msm(psor2.msm, ci="boot", B=2, 
                                    
                                    covariates=list(hieffusn=0, ollwsdrt="foo"))))
})

if (!covr::in_covr()){
  test_that("bootstrap CIs with parallel processing",{
    skip_on_cran()
    expect_error({
      t.par <- system.time(pmatrix.msm(psor.msm, ci="boot", cores=4, B=8))
                                        #    t.ser <- system.time(pmatrix.msm(psor.msm, ci="boot", cores=1, B=8))
                                        #    expect_true(t.par["elapsed"] < t.ser["elapsed"]) # if B is bigger 
    }, NA)
  })
}

