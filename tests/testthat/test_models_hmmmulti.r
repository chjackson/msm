context("HMMs with multivariate responses")

## Simulate data from a Markov model 
nsubj <- 30; nobspt <- 5
sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt),
                     time = seq(0, 20, length=nobspt))
set.seed(1)
two.q <- rbind(c(-0.1, 0.1), c(0, 0))
dat <- simmulti.msm(sim.df[,1:2], qmatrix=two.q, drop.absorb=FALSE)

## Bin(40, 0.1) for state 1, Bin(40, 0.5) for state 2
dat$obs1 <- dat$obs2 <- NA
set.seed(1)
dat$obs1[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.1)
dat$obs2[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.1)
dat$obs1[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.5)
dat$obs2[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.5)
dat$obs <- cbind(obs1 = dat$obs1, obs2 = dat$obs2)

dat$dobs1 <- dat$dobs2 <- NA
set.seed(1)
## Bin(40, 0.1) and Bin(40, 0.2) for state 1, 
dat$dobs1[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.1)
dat$dobs2[dat$state==1] <- rbinom(sum(dat$state==1), 40, 0.2)
## Bin(40, 0.5) and Bin(40, 0.6) for state 2
dat$dobs1[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.6)
dat$dobs2[dat$state==2] <- rbinom(sum(dat$state==2), 40, 0.5)
dat$dobs <- cbind(dobs1 = dat$dobs1, dobs2 = dat$dobs2)

options(msm.test.analytic.derivatives=TRUE)
err <- 1e-04

test_that("HMMs with multiple responses from the same distribution",{
    hmm <- msm(obs ~ time, subject=subject, data=dat, qmatrix=two.q,
               hmodel = list(hmmBinom(size=40, prob=0.2),
               hmmBinom(size=40, prob=0.2)), fixedpars=TRUE)
    expect_equal(hmm$minus2loglik, 4387.58552977954, tol=1e-06)
    expect_lt(deriv_error(hmm), err)
})

test_that("HMMs with multiple responses from different distributions",{
    hmm <- msm(dobs ~ time, subject=subject, data=dat, qmatrix=two.q,   
               hmodel = list(hmmMV(hmmBinom(size=40, prob=0.3),
               hmmBinom(size=40, prob=0.3)),                 
               hmmMV(hmmBinom(size=40, prob=0.3),
                     hmmBinom(size=40, prob=0.3))),
               fixedpars=TRUE)
    expect_equal(hmm$minus2loglik, 3767.11569380418, tol=1e-06)
    expect_lt(deriv_error(hmm), err)
})

test_that("HMMs with multiple responses from different distributions: non-default initprobs, different probs",{
    hmm <- msm(dobs ~ time, subject=subject, data=dat, qmatrix=two.q,
               initprobs=c(0.6, 0.4),
               hmodel = list(hmmMV(hmmBinom(size=40, prob=0.3),
                                   hmmBinom(size=40, prob=0.3)),                 
                             hmmMV(hmmBinom(size=40, prob=0.4),
                                   hmmBinom(size=40, prob=0.4))),
               fixedpars=TRUE)
    expect_lt(deriv_error(hmm), err)
})

dat$dobsmiss <- dat$dobs
dat$dobsmiss[1:10,2] <- NA

test_that("HMMs with multiple responses from different distributions: missing data",{
    hmm <- msm(dobsmiss ~ time, subject=subject, data=dat, qmatrix=two.q,   
               hmodel = list(hmmMV(hmmBinom(size=40, prob=0.3),
               hmmBinom(size=40, prob=0.3)),                 
               hmmMV(hmmBinom(size=40, prob=0.3),
                     hmmBinom(size=40, prob=0.3))),
               fixedpars=TRUE)
    expect_lt(deriv_error(hmm), err)
})

## not tested: different number of models by design for different obs 

dat$obstrue <- NA
obstimes <- seq(2, 147, by=5)  # times when true state is known
dat$obstrue[obstimes] <- dat$state[obstimes]

test_that("HMMs with multiple responses: true state known sometimes",{
    hmm <- msm(dobs ~ time, subject=subject, data=dat, qmatrix=two.q,
               obstrue=obstrue,              
               hmodel = list(hmmMV(hmmBinom(size=40, prob=0.3),
               hmmBinom(size=40, prob=0.3)),
               hmmMV(hmmBinom(size=40, prob=0.3),
                     hmmBinom(size=40, prob=0.3))),
               fixedpars=TRUE)
    expect_equal(hmm$minus2loglik, 3804.03972787726, tol=1e-06)
    expect_lt(deriv_error(hmm), err)
})

## analytic derivatives

options(msm.test.analytic.derivatives=NULL)


## censoring? 
