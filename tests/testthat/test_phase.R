test_that("one phased state", {
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, fixedpars=TRUE,
                  phase.states=c(1))
  expect_equal(psor.msm$minus2loglik, 1308.35590680014, tol=1e-06)
}
)

test_that("working example",{
  skip_on_cran()
  ## SIMULATE
  mst1 <- 5 # Short-stay mean 
  mst2 <- 30 # Long-stay mean 
  p2 <- 0.9 # Long-stay probability 
  q23 <- 1/5 # Transition rate between phases
  
  l1 <- (p2/mst1)
  mu1 <- (1-p2)/mst1
  mu2 <- 1/(mst2-mst1)
  Q <- rbind(c(0,l1,mu1*0.4,mu1*0.6),
             c(0,0,mu2*0.4,mu2*0.6),
             c(0,0,0,q23),
             c(0,0,0,0))
  # Given the hidden state, the observed state is deterministic
  E <- rbind(c(1,0,0,0),
             c(1,0,0,0),
             c(0,1,0,0),
             c(0,0,1,0))
  nsubj <- 1000; nobspt <- 10
  set.seed(1)
  sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 100, length=nobspt))
  sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=Q, ematrix=E)
  statetable.msm(obs, subject, sim2.df)
  
  ## FIT
  Q3 <- rbind(c(0,0.5,0.5),c(0,0,0.5),c(0,0,0))
  s.msm <- msm(obs ~ time, subject=subject, data=sim2.df, phase.states=1, qmatrix=Q3, # default inits
               phase.inits=list(list(trans=0.05, 
                                     exit=matrix(c(0.1,0.1,0.1,0.1),nrow=2))),
               control=list(trace=1,REPORT=1,fnscale=50000,maxit=10000),method="BFGS")
  
  # Parameter estimates should agree with the values used for simulation
  s.msm 
  c(l1, mu1*0.4, mu1*0.6, mu2*0.4, mu2*0.6, q23)
  means <- phasemeans.msm(s.msm)
  
  expect_equal(means[[1,"Short stay mean"]], 5, tol=1)
  expect_equal(means[[1,"Long stay mean"]], 30, tol=5)
  expect_equal(means[[1,"Long stay probability"]], 0.9, tol=0.1)
})


test_that("two phased states", {
  ## note phase transition within state 2 appears hard to estimate.
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, fixedpars=TRUE,
                  phase.states=c(1,2))
  expect_equal(psor.msm$minus2loglik, 1303.768173132893, tol=1e-06)
})


test_that("supplying initial values",{
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1,2),
                  phase.inits = list(
                    list(trans=c(0.42032), exit=c(0.17404,0.04809)),
                    list(trans=c(0.42032), exit=c(0.17404,0.04809))),                 
                  fixedpars=TRUE)
  expect_equal(psor.msm$minus2loglik, 1280.788860486683, tol=1e-06)
})

test_that("errors in initial values",{
  expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1,2), phase.inits = "foo"), "phase.inits should be a list")
  expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1,2), phase.inits = list(1,2,3)), "phase.inits of length 3, but there are 2 phased states")
  expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1,2), phase.inits = list(1,2)), "phase.inits.+ list of length 1, should be 2")
  expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1), phase.inits = list(list(trans=c(0.42032, 1), exit=c(0.17404,0.04809)))), "phase.inits.+trans of length 2, should be 1")
  expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1), phase.inits = list(list(trans=c(0.42032), exit=c(0.17404,0.04809,1,1)))), "phase.inits.+exit has 4 columns, should be 2")
  expect_error(psor2.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, phase.states=c(1), phase.inits = list(list(trans=c(0.42032), exit=matrix(c(0.17404,0.04809,1,1),ncol=2,byrow=TRUE)))), "phase.inits.+exit has 2 rows, but there are 1 exit states from this state")
})


test_that("with censoring", {
  psor2 <- psor
  psor2$state[c(16,18)] <- 99
  psorc.msm <- msm(state ~ months, subject=ptnum, data=psor2, qmatrix = psor.q, fixedpars=TRUE,
                   censor = 99, censor.states=c(3,4))
  expect_equal(psorc.msm$minus2loglik, 1290.189703452163, tol=1e-06)
})


test_that("HMM on top", {
  ### Note 1.6 clarifies obstrue not supported for hmmCat: so test result in 1.5 wrong
  miscnew.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                     qmatrix = oneway4.q, death = 5, #obstrue=firstobs,
                     fixedpars=TRUE,
                     # could try to fit, but weakly identifiable.  fit seems improved though (-2LL 3864 vs 3951)
                     # control=list(trace=1,REPORT=1,fnscale=4000, maxit=10000),
                     phase.states = 1, 
                     hmodel=list(
                       hmmCat(prob=c(0.9, 0.1, 0, 0)),
                       hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                       hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent(4))
  )
  expect_equal(miscnew.msm$minus2loglik, 4357.70908245511, tol=1e-06)  
})
