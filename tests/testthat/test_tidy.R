test_that("tidy.msm",{
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  
                  covariates = ~ollwsdrt+hieffusn, 
                  constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
  x <- tidy(psor.msm)
  expect_equal(x$estimate[x$parclass=="hr" & x$state==2 & x$term=="hieffusn"],
               hazard.msm(psor.msm)$hieffusn[2,"HR"])
  
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q)
  x <- tidy(psor.msm)
  expect_equal(x$conf.low[x$parclass=="intens" & x$state==2],
               qmatrix.msm(psor.msm)[["L"]][2,3])
})

test_that("tidy.msm with covariates called baseline or Baseline",{
  psor$baseline <- psor$logbaseline <- psor$ollwsdrt
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  
                  covariates = ~baseline)
  x <- tidy(psor.msm)
  expect_equal(x$estimate[x$parclass=="hr" & x$state==2 & x$term=="baseline"],
               hazard.msm(psor.msm)$baseline[2,"HR"])

  psor$Baseline <- psor$ollwsdrt
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  
                  covariates = ~Baseline)
  psor.msm
  x <- tidy(psor.msm)
  expect_equal(x$estimate[x$parclass=="hr" & x$state==2 & x$term=="Baseline"],
               hazard.msm(psor.msm)$Baseline[2,"HR"])
})

test_that("tidy.msm with misclassification models",{
  misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                     qmatrix = oneway4.q, ematrix=ematrix, deathexact = 4, fixedpars=TRUE,
                     misccovariates = ~dage + sex, 
                     misccovinits = list(dage=c(0.01,0.02,0.03,0.04),
                                         sex=c(-0.013,-0.014,-0.015,-0.016)))
  x <- tidy(misccov.msm)
  expect_equal(x$estimate[x$parclass=="misc" & x$state==2 & x$tostate==1],
               ematrix.msm(misccov.msm)[2,1][["estimate"]])
  
  misccov.msm <- suppressWarnings(msm(state ~ years, subject = PTNUM, data = cav[1:500,],
                     qmatrix = oneway4.q, ematrix=ematrix, 
                     deathexact = 4, fixedpars=FALSE, control=list(maxit=20),
                     misccovariates = ~dage + sex))
  x <- tidy(misccov.msm)
  expect_equal(x$conf.high[x$parclass=="misc" & x$state==2 & x$tostate==1],
               ematrix.msm(misccov.msm)[["U"]][2,1])
  
  x <- ematrix.msm(misccov.msm)
  tx <- tidy(x)
  expect_equal(tx$conf.high[tx$state==3 & tx$tostate==2], x[["U"]][3,2])
  
  x <- ematrix.msm(misccov.msm, ci="none")
  tx <- tidy(x)
  expect_equal(tx$estimate[tx$state==3 & tx$tostate==2], x[3,2])
})

psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  
                covariates = ~ollwsdrt+hieffusn, 
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))

test_that("Tidying extractor function output",{
  x <- qmatrix.msm(psor.msm)
  tx <- tidy(x)
  expect_equal(tx$conf.low[tx$state==3], x[["L"]][3,4])

  x <- qmatrix.msm(psor.msm, ci="none")
  tx <- tidy(x)
  expect_equal(tx$estimate[tx$state==3], x[3,4])
  
  x <- pmatrix.msm(psor.msm)
  tx <- tidy(x)
  expect_equal(tx$estimate[tx$state==3 & tx$tostate==4], x[3,4])
  
  x <- pmatrix.msm(psor.msm, ci="normal", B=10)
  tx <- tidy(x)
  expect_equal(tx$conf.low[tx$state==3 & tx$tostate==4], x[["L"]][3,4])
  
  x <- pnext.msm(psor.msm)
  tx <- tidy(x)
  expect_equal(tx$estimate[tx$state==3 & tx$tostate==4], x[["estimates"]][3,4])
  
  x <- pnext.msm(psor.msm, ci="none")
  tx <- tidy(x)
  expect_equal(tx$estimate[tx$state==3 & tx$tostate==4], x[["estimates"]][3,4])
  
  x <- ppass.msm(psor.msm, tot=3)
  tx <- tidy(x)
  expect_equal(tx$estimate[tx$state==3 & tx$tostate==4], x[3,4])
  
  x <- ppass.msm(psor.msm, tot=3, ci="normal", B=10)
  tx <- tidy(x)
  expect_equal(tx$conf.low[tx$state==3 & tx$tostate==4], x[["L"]][3,4])
})



test_that("tidy for HMMs",{
  three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))
  hmodel3 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))
  fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev[1:1000,], 
                   qmatrix=three.q, 
                   deathexact=3, hmodel=hmodel3,
                   hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                   fixedpars=FALSE, center=FALSE)
  tx <- tidy(fev3.hid)
  expect_equal(tx$estimate[tx$parclass=="hmm" & tx$state==2 & tx$term=="mean"],
               fev3.hid$hmodel$pars[[3]])
  expect_equal(tx$conf.high[tx$parclass=="hmm" & tx$state==2 & tx$term=="sd"],
               fev3.hid$hmodel$ci[4,][2])
  
  # initprobs 
  fev3.hid <- suppressWarnings(msm(fev ~ days, subject=ptnum, data=fev[1:2000,], 
                   qmatrix=three.q, deathexact=3, hmodel=hmodel3, 
                   est.initprobs = TRUE, control=list(maxit=10)))
  tx <- tidy(fev3.hid)
  expect_equal(tx$estimate[tx$parclass=="initp" & tx$state==2],
               fev3.hid$hmodel$initprobs["State 2","Estimate"])
  expect_equal(tx$conf.low[tx$parclass=="initp" & tx$state==2],
               fev3.hid$hmodel$initprobs["State 2","LCL"])
  
  # initcovariates
  fev3.hid <- suppressWarnings(msm(fev ~ days, subject=ptnum, data=fev[1:3000,], 
                   qmatrix=three.q, exacttimes=TRUE, 
                   deathexact=3, hmodel=hmodel3, initcovariates = ~acute,
                   est.initprobs = TRUE, control=list(maxit=100)))
  tx <- tidy(fev3.hid)
  expect_equal(tx$estimate[tx$parclass=="initpcov" & tx$state==2],
               fev3.hid$hmodel$icoveffect[[1]])
})

test_that("tidy prevalence.msm",{
  x <- prevalence.msm(psor.msm)
  tx <- tidy(x)
  expect_equal(tx$expected[tx$output=="count" & tx$time==tx$time[8] & tx$state==2],
               x$Expected[8,"State 2"])
  
  x <- prevalence.msm(psor.msm, ci="normal", B=3)
  tx <- tidy(x)
  expect_equal(tx$conf.high[tx$output=="count" & tx$time==tx$time[8] & tx$state==2],
               unname(x$Expected$ci[8,2,"97.5%"]))
})

test_that("tidy totlos.msm",{
  x <- totlos.msm(psor.msm)
  tx <- tidy(x)
  expect_equal(tx$estimate[tx$state==3], x[["State 3"]])

  x <- totlos.msm(psor.msm, ci="normal", B=3)
  tx <- tidy(x)
  expect_equal(tx$conf.high[tx$state==3], x["97.5%",][["State 3"]])
  
  x <- efpt.msm(psor.msm, tostate = 3)
  tx <- tidy(x)
  expect_equal(x[3], tx$estimate[tx$state==3])
    
  x <- efpt.msm(psor.msm, tostate = 3, ci="normal", B=3)
  tx <- tidy(x)
  expect_equal(unname(x[3,2]), tx$conf.high[tx$state==2])
  
  x <- envisits.msm(psor.msm)
  tx <- tidy(x)
  expect_equal(tx$estimate[tx$state==3], x[["State 3"]])
})


