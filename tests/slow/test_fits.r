### Slow model fits: only run locally

context("msm simple model fits")

cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE)

test_that("cav model fit",{
    expect_equal(3968.7978930519, cav.msm$minus2loglik)
})

test_that("cav model fit with covariates",{
    cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE, covariates = ~ sex)
    expect_equal(3954.77700438268, cavsex.msm$minus2loglik)
    ## previously tested output functions here, but unnecessary if lik is right
})

test_that("cav model fit with qconstraints",{
    cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE,qconstraint = c(1,1,2,2,2,3,3))
    expect_equal(4116.22686367935, cav.msm$minus2loglik)
})

test_that("cav model fit with covariate constraints",{
    cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE,qconstraint = c(1,1,2,2,2,3,3))
    cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, covariates = ~ sex, covinits = list(sex=rep(0.01, 7)), constraint=list(sex=c(1,2,3,1,2,3,2)))
    expect_equal(3959.35551766943, cav.msm$minus2loglik)
})

test_that("some covariate effects constrained to be minus others",{
    psor.constr.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,-1,1),ollwsdrt=c(-1,1,2)))
    expect_equal(psor.constr.msm$Qmatrices$hieffusn[1,2], -psor.constr.msm$Qmatrices$hieffusn[2,3])
    expect_equal(psor.constr.msm$Qmatrices$ollwsdrt[1,2], -psor.constr.msm$Qmatrices$ollwsdrt[2,3])
    psor.constr.msm <- msm(state ~ months, subject=ptnum, data=psor,  qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(-10, 10, -10),ollwsdrt=c(-2,2,1)))
    expect_equal(1129.538, psor.constr.msm$minus2loglik, tol=1e-05)
    expect_equal(psor.constr.msm$Qmatrices$hieffusn[1,2], -psor.constr.msm$Qmatrices$hieffusn[2,3])
    expect_equal(psor.constr.msm$Qmatrices$ollwsdrt[1,2], -psor.constr.msm$Qmatrices$ollwsdrt[2,3])
})

test_that("qconstraints and constraints together: bug in <=0.9.2",{
    psorqc.msm <- msm(state ~ months, subject=ptnum, data=psor,
                      qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                      qconstraint = c(1,2,2),
                       constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
    stopifnot(isTRUE(all.equal(1120.04586677803, psorqc.msm$minus2loglik, tol=1e-06)))

    psorfix.msm <- msm(state ~ months, subject=ptnum, data=psor,
                       qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                       fixedpars=6,
                       constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
    stopifnot(isTRUE(all.equal(1121.71012633891, psorfix.msm$minus2loglik, tol=1e-06)))
})

test_that("MLE for exact transition times is same as crudeinits",{
    fiveq.i <- fiveq; fiveq.i[fiveq.i!=0] <- 1
    (msmtest.fix <- msm(state ~ time, qmatrix = fiveq.i, gen.inits=TRUE, subject = ptnum, data = bos, exacttimes=TRUE, fixedpars=TRUE))
    (msmtest.nofix <- msm(state ~ time, qmatrix = fiveq.i, gen.inits=TRUE, subject = ptnum, data = bos, exacttimes=TRUE))
    expect_equal(msmtest.nofix$minus2loglik, msmtest.fix$minus2loglik, tol=1e-06)
})

test_that("simple model, covariates on specific intensities",{
    expect_equal(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, fixedpars=c(4,5,7))$minus2loglik,
                 msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn))$minus2loglik)
    expect_error(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "4-1" = ~hieffusn)), "not permitted by the qmatrix")
})

test_that("simple model, covariates on specific intensities and additional fixedpars",{
    expect_equal(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, fixedpars=c(4,5,7,9))$minus2loglik,
                 msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), fixedpars=6)$minus2loglik)
    expect_error(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), fixedpars=7), "Elements of fixedpars should")
})

test_that("simple model, covariates on specific intensities, constraints",{
    ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), constraint=list(hieffusn=c(1,1)))
    expect_equal(ps$Qmatrices$hieffusn[2,3], ps$Qmatrices$hieffusn[3,4])
    expect_error(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), constraint=list(hieffusn=c(1,1,1))),"constraint of length")
    expect_error(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), constraint="foo"),"constraint should be a list")
    expect_error(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), constraint=list(hieffusn=c("foo",1)),"constraint should be a list of numeric vectors"))
})

test_that("simple model, covariates on specific intensities, covinits",{
    expect_equal(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), covinits=list(hieffusn=c(0.1,0.1)), fixedpars=TRUE)$minus2loglik,
                 msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, covinits=list(hieffusn=c(0,0.1,0.1)), fixedpars=TRUE)$minus2loglik)
    expect_warning(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), covinits=list(hieffusn=c(0, 0.1,0.1)), fixedpars=TRUE),"initial values of length 3, should be 2") # covinits wrong length
    expect_warning(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), covinits=list(hieffusn=c("foo")), fixedpars=TRUE),"covinits should be a list of numeric")
})

test_that("simple model, covariates on specific intensities, factors",{
    set.seed(12082012)
    psor$faccov <- factor(sample(1:3, size=nrow(psor), replace=TRUE))
    ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn+faccov))
    expect_equal(ps$Qmatrices$faccov2[1,2],0)
    expect_equal(ps$Qmatrices$faccov2[1,2],ps$Qmatrices$faccov3[1,2])
    expect_error(ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn+faccov, "2-3" = ~hieffusn+faccov), constraint=list(faccov=c(1,1))),"Covariate .+ not in model")
    ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn+faccov, "2-3" = ~hieffusn+faccov), constraint=list(faccov2=c(1,1)))
    expect_equal(ps$Qmatrices$faccov2[2,3],ps$Qmatrices$faccov2[3,4])
})

test_that("simple model, covariates on specific intensities, interactions",{
    ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt*hieffusn, "2-3" = ~hieffusn*ollwsdrt), method="BFGS")
    expect_equal(-0.780889267961125, ps$Qmatrices$"ollwsdrt:hieffusn"[2,3], tol=1e-05)
})

test_that("model fits with versus without derivatives",{
    sd <- system.time(cav.d.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE, use.deriv=TRUE))
    sn <- system.time(cav.nd.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE, use.deriv=FALSE))
    expect_true(sd["elapsed"] < sn["elapsed"])
    expect_equal(cav.nd.msm$minus2loglik, cav.d.msm$minus2loglik)
    expect_equal(cav.nd.msm$estimates, cav.d.msm$estimates, tol=1e-04)
    expect_equal(cav.nd.msm$ci, cav.d.msm$ci, tol=1e-04)
    sd <- system.time(cav.d.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, use.deriv=TRUE, covariates = ~ sex, fixedpars=c(5,12)))
    sn <- system.time(cav.nd.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, use.deriv=FALSE, covariates = ~ sex, fixedpars=c(5,12)))
    expect_true(sd["elapsed"] < sn["elapsed"])
    expect_equal(cav.nd.msm$minus2loglik, cav.d.msm$minus2loglik)
    expect_equal(cav.nd.msm$estimates, cav.d.msm$estimates, tol=1e-06)
    expect_equal(cav.nd.msm$ci, cav.d.msm$ci, tol=1e-06)
})

test_that("nlm optimisation method",{
    sn <- system.time(cav.nd.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE,use.deriv=FALSE, opt.method="nlm")) # couple of warnings: ignore.   also slower with derivs
    so <- system.time(cav.od.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE))
    expect_equal(cav.nd.msm$minus2loglik, cav.od.msm$minus2loglik)
    expect_equal(cav.nd.msm$estimates, cav.od.msm$estimates, tol=1e-03)
    expect_equal(cav.nd.msm$ci, cav.od.msm$ci, tol=1e-03)
})

test_that("fisher scoring method",{
    cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q)
    cav.fish.msm <- msm( state ~ years, subject=PTNUM, data = cav, opt.method="fisher", qmatrix = twoway4.q)
    expect_equal(cav.fish.msm$minus2loglik, cav.msm$minus2loglik)
    expect_equal(cav.fish.msm$paramdata$estimates, cav.msm$paramdata$estimates, tol=1e-04)
})

if (require("minqa")){
test_that("bobyqa method",{
    system.time(psor.opt.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,
                                covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2))))
    system.time(psor.bob.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, opt.method="bobyqa",
                                covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2))))
    expect_equal(psor.opt.msm$minus2loglik, psor.bob.msm$minus2loglik)
    expect_error(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, opt.method="foo"),"Unknown optimisation method")
})
}

test_that("piecewise constant intensities",{
    cav5.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                    qmatrix = twoway4.q, death = TRUE, pci = c(5), covinits = list("timeperiod[5,Inf)"=rep(0.01,7)))
    expect_equal(3919.99601396194, cav5.msm$minus2loglik, tol=1e-06)
    cav10.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     qmatrix = twoway4.q, death = TRUE, pci = c(5,10),
                     covinits = list("timeperiod[5,10)"=rep(0.01,7), "timeperiod[10,Inf)"=rep(0.01,7)))
    expect_equal(3882.08834017773, cav10.msm$minus2loglik, tol=1e-06)
})

test_that("plot.survfit.msm",{
    if (interactive()){ 
        plot.survfit.msm(cav.msm)
        plot.survfit.msm(cav.msm, from=2)
        plot.survfit.msm(cav.msm, from=3)
        plot.survfit.msm(cav.msm, interp="midpoint")
        sd <- plot.survfit.msm(cav.msm, survdata=TRUE)
        plot(survfit(Surv(survtime, died) ~ 1, data=sd))
    }
})
