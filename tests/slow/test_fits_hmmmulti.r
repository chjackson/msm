context("HMMs with multivariate responses: model fitting")

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

test_that("viterbi with multiple responses from different distributions",{
    hmm <- msm(dobs ~ time, subject=subject, data=dat, qmatrix=two.q,   
               hmodel = list(hmmMV(hmmBinom(size=40, prob=0.3),
               hmmBinom(size=40, prob=0.3)),                 
               hmmMV(hmmBinom(size=40, prob=0.3),
                     hmmBinom(size=40, prob=0.3))),
               control=list(maxit=10000))
    expect_equal(hmm$minus2loglik, 1467.06582421925, tol=1e-06)
    vit <- viterbi.msm(hmm)
    expect_equal(as.numeric(table(dat$state, vit$fitted)), c(73, 0, 0, 77))
})
