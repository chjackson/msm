context("slow model outputs")

psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))

test_that("totlos.msm",{
    tot <- totlos.msm(psor.msm, tot=Inf)
    expect_equal(c(10.4237238934583, 6.08618568645837, 3.93108395480513, Inf), as.numeric(tot), tol=1e-04)
    tot <- totlos.msm(psor.msm, fromt=1, tot=30)
    expect_equal(c(8.88391824336528, 5.28131329931934, 3.18628092862045), as.numeric(tot[1:3]), tol=1e-04)
    tot <- totlos.msm(psor.msm, start=2, fromt=10, tot=30, end=1:3)
    expect_equal(c(0, 1.1329680757691, 1.50674417311231), as.numeric(tot), tol=1e-04)
    expect_equal(-557.449730608585, as.numeric(logLik.msm(psor.msm)), tol=1e-04)
    expect_error(totlos.msm("foo"))
    expect_error(totlos.msm(psor.msm, start="foo"))
    expect_error(totlos.msm(psor.msm, start=-1))
    expect_error(totlos.msm(psor.msm, start=1, fromt=1, tot=0))
    expect_error(totlos.msm(psor.msm, start=1, fromt=c(1,2), tot=c(3,4)))
    expect_error(totlos.msm(psor.msm, start=1, fromt="foo", tot=2))
    expect_error(totlos.msm(psor.msm, start=1, fromt=-3, tot=-2))
})

test_that("totlos.msm analytic matches numeric integration",{
    f <- function(time) {
        y <- numeric(length(time))
        for (i in seq(along=y))
            y[i] <- pmatrix.msm(psor.msm, time[i], ci="none")[1,2]
        y
    }
    expect_equivalent(integrate(f, 0, 10)$value, totlos.msm(psor.msm, tot=10)[2])

    expect_equal(totlos.msm(psor.msm, tot=10), totlos.msm(psor.msm, tot=10, num.integ=TRUE))
    expect_equal(totlos.msm(psor.msm, fromt=1, tot=10), totlos.msm(psor.msm, fromt=1, tot=10, num.integ=TRUE))
    expect_equal(totlos.msm(psor.msm, tot=100), totlos.msm(psor.msm, tot=100, num.integ=TRUE))
    expect_equal(totlos.msm(psor.msm, tot=1000), totlos.msm(psor.msm, tot=1000, num.integ=TRUE))
    expect_equal(totlos.msm(psor.msm, tot=10000), totlos.msm(psor.msm, tot=10000, num.integ=TRUE))
})

test_that("factor covariates", {
    cav$pdiag2 <- as.character(cav$pdiag); cav$pdiag2[cav$pdiag %in% c("CVCM","Hyper","Other","Restr")] <- "Other"; cav$pdiag2 <- factor(cav$pdiag2, levels=c("IDC","IHD","Other"))
    cavfaccov.msm <- msm(state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, covariates = ~ pdiag2)
    cavfaccov2.msm <- msm(state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q,  covariates = ~ sex + pdiag2)

    expect_equal(msm:::factorcov2numeric.msm(list(pdiag2 = "IHD", sex=1), cavfaccov2.msm), list(sex=1,pdiag2IHD=1,pdiag2Other=0))
    expect_equal(msm:::factorcov2numeric.msm(list(sex=1, pdiag2 = "IHD"), cavfaccov2.msm), list(sex=1,pdiag2IHD=1,pdiag2Other=0))
    expect_warning(msm:::factorcov2numeric.msm(list(1, pdiag2 = "IHD"), cavfaccov2.msm), "Covariate .+ unknown")
    expect_warning(msm:::factorcov2numeric.msm(list(pdiag2 = "IHD", 1), cavfaccov2.msm), "Covariate .+ unknown")
    expect_equal(msm:::factorcov2numeric.msm(list(1, "IHD"), cavfaccov2.msm), list(sex=1,pdiag2IHD=1,pdiag2Other=0))
    expect_error(msm:::factorcov2numeric.msm(list("IHD", 1), cavfaccov2.msm), "Level \"1\" of covariate pdiag2 unknown")

    expect_equal(qmatrix.msm(cavfaccov.msm, covariates=0), qmatrix.msm(cavfaccov.msm, covariates=list(pdiag2="IDC")))
    expect_equal(qmatrix.msm(cavfaccov2.msm, covariates=0), qmatrix.msm(cavfaccov2.msm, covariates=list(pdiag2="IDC", sex=0)))
    expect_error(qmatrix.msm(cavfaccov.msm, covariates=list(pdiag2="Nonexistent")), "Level .+ unknown")
    expect_warning(qmatrix.msm(cavfaccov.msm, covariates=list(nonex="Nonexistent")), "Covariate .+ unknown")
    expect_warning(qmatrix.msm(cavfaccov.msm, covariates=list(pdiag2="IDC", foo="bar")), "Covariate .+ unknown")
    expect_error(qmatrix.msm(cavfaccov2.msm, covariates=list("IDC")), "Expected covariate list")
    expect_error(qmatrix.msm(cavfaccov2.msm, covariates=list("IDC", 0)), "Level .+ unknown")
})

test_that("pearson test",{
    ## panel data
    pears <- pearson.msm(psor.msm, pval=FALSE)
    expect_equal(4027, pears$test$stat, tol=1)
    ## exact death times
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, death=c(4), constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
    pears <- pearson.msm(psor.msm)
    expect_equal(1799, pears$test$stat, tol=1)
    pearson.msm(psor.msm, boot=TRUE, B=3)
    ## multiple death states
    psor2.q <- rbind(c(0,0.1,0,0,0),c(0,0,0.1,0,0),c(0,0,0,0.1,0.1),c(0,0,0,0,0),c(0,0,0,0,0))
    psor2 <- psor; psor2$state[psor2$state==4][1:10] <- 5
    psor.msm <- msm(state ~ months, subject=ptnum, data=psor2, qmatrix = psor2.q, covariates = ~ollwsdrt+hieffusn, death=c(4,5),  constraint = list(hieffusn=c(1,1,1,1),ollwsdrt=c(1,1,2,1)))
    pears <- pearson.msm(psor.msm)
    expect_equal(2047, pears$test$stat, tol=1)
})
