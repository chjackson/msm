
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
