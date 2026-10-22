test_that("simfitted.msm",{
  skip_on_cran()
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt)
  expect_equal(simfitted.msm(psor.msm)$time[1], psor$months[1])
 
  psorc <- psor
  psorc$state[2] <- 99
  psorc.msm <- msm(state ~ months, subject=ptnum, data=psorc, qmatrix = psor.q,
                   covariates = ~ollwsdrt, censor=99)
  set.seed(1)
  simdat <- simfitted.msm(psorc.msm)$state[2]
  expect_equal(simdat, psorc$state[2])
})

simul.q3 <- rbind(c(0, 0.1, 0.02), c(0.1, 0, 0.05), c(0, 0, 0))
simul_data <- function(N) {
  d <- data.frame(subject=rep(1:N, each=7), time=rep(c(0,1,2,4,7,10,15), N))
  d$z <- rbinom(nrow(d), 1, 0.4)
  d
}

## Compare the distribution of the simulated state at each time between the
## C engine and the original R implementation
expect_same_dist <- function(a, b, by=c("time","state")) {
  ka <- do.call(paste, a[by]); kb <- do.call(paste, b[by])
  lev <- union(ka, kb)
  tab <- rbind(table(factor(ka, lev)), table(factor(kb, lev)))
  tab <- tab[, colSums(tab) >= 10, drop=FALSE]
  expect_gt(suppressWarnings(chisq.test(tab)$p.value), 0.001)
}

test_that("simmulti.msm C engine matches original R implementation", {
  skip_on_cran()
  set.seed(1)
  d <- simul_data(2000)
  cv <- list(z=c(1, 0.5, -1, 0.7))
  expect_same_dist(simmulti.msm(d, simul.q3, covariates=cv),
                   simmulti_legacy(d, simul.q3, covariates=cv))
  expect_same_dist(simmulti.msm(d, simul.q3, covariates=cv, death=3),
                   simmulti_legacy(d, simul.q3, covariates=cv, death=3))
  expect_same_dist(simmulti.msm(d[,1:2], simul.q3, drop.absorb=FALSE),
                   simmulti_legacy(d[,1:2], simul.q3, drop.absorb=FALSE))
  em <- rbind(c(0,.1,0), c(.1,0,0), c(0,0,0))
  d$x <- rnorm(nrow(d))
  expect_same_dist(simmulti.msm(d, simul.q3, ematrix=em, misccovariates=list(x=c(0.5,-0.5))),
                   simmulti_legacy(d, simul.q3, ematrix=em, misccovariates=list(x=c(0.5,-0.5))),
                   by=c("time","state","obs"))
})

test_that("simmulti.msm output structure", {
  set.seed(1)
  d <- simul_data(20)
  d$other <- letters[1:7]
  sim <- simmulti.msm(d, simul.q3, covariates=list(z=c(1, 0.5, -1, 0.7)))
  expect_equal(names(sim), c("subject","time","state","z","other","keep"))
  expect_equal(sim$time, d$time[sim$keep])
  expect_equal(sim$other, d$other[sim$keep])
  expect_true(all(sim$state %in% 1:3))
  ## only the first observation in the absorbing state is kept
  dead <- tapply(sim$state == 3, sim$subject, sum)
  expect_true(all(dead <= 1))
  ## no transitions possible: all states equal start
  sim0 <- simmulti.msm(d, matrix(0, 3, 3), start=2, drop.absorb=FALSE)
  expect_equal(nrow(sim0), nrow(d))
  expect_true(all(sim0$state == 2))
  ## exact death times
  simd <- simmulti.msm(d, simul.q3, death=3, drop.absorb=FALSE)
  expect_equal(nrow(simd), nrow(d))
  ## starting in the absorbing state
  sima <- simmulti.msm(d, simul.q3, start=3)
  expect_equal(nrow(sima), 2 * 20)
})

test_that("simmulti.msm data checks", {
  d <- simul_data(5)
  expect_warning(simmulti.msm(d[nrow(d):1,], simul.q3), "not ordered")
  expect_warning(simmulti.msm(d[c(1,1:nrow(d)),], simul.q3), "duplicated")
  d2 <- d; d2$time[2:3] <- d2$time[3:2]
  expect_error(simmulti.msm(d2, simul.q3), "not ordered by time")
  expect_error(msm.check.times(c(0,1,0,1), c(1,2,1,2)), "not adjacent")
})

test_that("sim.msm C engine matches original R implementation", {
  skip_on_cran()
  z <- cbind(z=c(0,1,1,0,1,0,0,1,1,0))
  summ <- function(f) t(replicate(3000, { s <- f()
    c(n=length(s$states), last=s$states[length(s$states)], tend=s$times[length(s$times)]) }))
  set.seed(1)
  a <- summ(function() sim.msm(simul.q3, 12, covs=z, beta=c(1,0.5,-1,0.7), obstimes=0:9))
  b <- summ(function() sim_msm_legacy(simul.q3, 12, covs=z, beta=c(1,0.5,-1,0.7), obstimes=0:9))
  expect_same_dist(as.data.frame(a), as.data.frame(b), by=c("n","last"))
  expect_gt(suppressWarnings(ks.test(a[,"tend"], b[,"tend"])$p.value), 0.001)
})

test_that("sim.msm output structure", {
  set.seed(1)
  s <- sim.msm(simul.q3, 30)
  expect_equal(s$times[1], 0)
  expect_false(is.unsorted(s$times))
  expect_true(s$states[length(s$states)] == 3 || s$times[length(s$times)] == 30)
  expect_equal(s$qmatrix, simul.q3)
  ## censoring row repeats the last state
  s <- sim.msm(rbind(c(0, 1e-10), c(0, 0)), 5, mintime=2)
  expect_equal(s$states, c(1, 1))
  expect_equal(s$times, c(2, 5))
  ## starting in an absorbing state
  expect_equal(sim.msm(simul.q3, 30, start=3)[c("states","times")], list(states=3, times=0))
  expect_equal(sim.msm(simul.q3, 5, mintime=5)[c("states","times")], list(states=1, times=5))
  ## covariates only change the intensities from their obstimes
  s <- sim.msm(rbind(c(0, 1), c(0, 0)), 10, covs=cbind(x=c(-100, 100)), beta=1, obstimes=c(0, 4))
  expect_gte(s$times[2], 4)
  expect_lt(s$times[2], 4.1)
})
