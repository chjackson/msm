test_that("subject weights",{
  cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                  qmatrix = twoway4.q, deathexact = TRUE, fixedpars=TRUE)
  cav$swt <- 1
  cavwt1.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                  qmatrix = twoway4.q, deathexact = TRUE, fixedpars=TRUE,
                  subject.weights = swt)
  
  expect_equal(cav.msm$minus2loglik, cavwt1.msm$minus2loglik)
  cav$swt[cav$PTNUM <= 100020] <- 1.2
  
  cavwt2.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     qmatrix = twoway4.q, deathexact = TRUE, fixedpars=TRUE,
                     subject.weights = swt)
  
  expect_true(cav.msm$minus2loglik != cavwt2.msm$minus2loglik)

  lik <- logLik(cav.msm, by.subject=TRUE)
  w <- cav$swt[!duplicated(cav$PTNUM)]
  expect_equivalent(cavwt2.msm$minus2loglik, - 2*sum(w*lik))
  
  cav$swt[1] <- "bad"
  expect_error(msm( state ~ years, subject=PTNUM, data = cav,
                    qmatrix = twoway4.q, deathexact = TRUE, fixedpars=TRUE,
                    subject.weights = swt), "must be numeric")

  cav$swt <- 1
  cav$swt[1:2] <- c(1.2, 1.3)
  expect_warning(msm( state ~ years, subject=PTNUM, data = cav,
                    qmatrix = twoway4.q, deathexact = TRUE, fixedpars=TRUE,
                    subject.weights = swt), "non-unique")

  cav$swt <- 1
  cav$swt[2] <- NA
  expect_no_warning(msm( state ~ years, subject=PTNUM, data = cav,
                        qmatrix = twoway4.q, deathexact = TRUE, fixedpars=TRUE,
                        subject.weights = swt))

})

test_that("subject weights, model fitting",{
  skip_on_cran()
  cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     qmatrix = twoway4.q, deathexact = TRUE, fixedpars=FALSE)

  cav$swt <- 1
  cavwt1.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     # work around boot not finding this in full test()
                     qmatrix = rbind(c(-0.5,0.25,0,0.25),c(0.166, -0.498, 0.166, 0.166), 
                                     c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0)),  
                     
                     deathexact = TRUE, fixedpars=FALSE,
                     subject.weights = swt)
  expect_equal(cav.msm$minus2loglik, cavwt1.msm$minus2loglik)
  pmatrix.msm(cavwt1.msm, ci="normal")
  set.seed(1)
  pmatrix.msm(cavwt1.msm, ci="boot", B=3)
  pearson.msm(cavwt1.msm, boot = TRUE, B=2)
  
  cav$swt[cav$PTNUM <= 100020] <- 1.2
  cavwt2.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     qmatrix = twoway4.q, deathexact = TRUE, fixedpars=FALSE,
                     subject.weights = swt)

  expect_true(cav.msm$minus2loglik != cavwt2.msm$minus2loglik)
})
