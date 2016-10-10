context("analytic transition probability matrices")

fixq <- function(Q){ diag(Q) <- 0;  diag(Q) <- - rowSums(Q); Q } # avoid namespace faff

nsubj <- 50; nobspt <- 6
sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length=nobspt), x = rnorm(nsubj*nobspt), y = rnorm(nsubj*nobspt)* 5 + 2 )
set.seed(22061976)

test_that("2 state analytic P matrices",{
    (two.q <- fixq(rbind(c(0, exp(-2)), c(0, 0))))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=two.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-3)), c(0, 0)), analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-3)), c(0, 0)), analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    ## 1,2
    (two.q <- fixq(rbind(c(0, exp(-2)), c(exp(-2), 0))))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=two.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-3)), c(exp(-1), 0)), analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-3)), c(exp(-1), 0)), analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
})

test_that("3 state analytic P matrices",{
    ## 1,2
    (three.q <- fixq(rbind(c(0, exp(-3), exp(-3)), c(0, 0, 0), c(0, 0, 0))))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-1)), c(0, 0, 0), c(0, 0, 0)), analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-1)), c(0, 0, 0), c(0, 0, 0)), analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    ## 1,4
    (three.q <- fixq(rbind(c(0, exp(-3), 0), c(0, 0, exp(-3)), c(0, 0, 0))))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0), c(0, 0, exp(-1)), c(0, 0, 0)), analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0), c(0, 0, exp(-1)), c(0, 0, 0)), analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

                                        # 4,5 (== 1,4)
    nsubj <- 500; nobspt <- 6
    sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length=nobspt),
                         x = rnorm(nsubj*nobspt), y = rnorm(nsubj*nobspt)* 5 + 2 )
    (three.q <- fixq(rbind(c(0, 0, 0), c(0, 0, exp(-2)), c(exp(-3), 0, 0))))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, start=rep(2,500))
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, 0, 0), c(0, 0, exp(-3)), c(exp(-2), 0, 0)), analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, 0, 0), c(0, 0, exp(-3)), c(exp(-2), 0, 0)), analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    ## 1,6
    (three.q <- fixq(rbind(c(0, exp(-3), 0), c(0, 0, 0), c(0, exp(-3), 0))))
    nsubj <- 50; nobspt <- 6
    sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length=nobspt), x = rnorm(nsubj*nobspt), y = rnorm(nsubj*nobspt)* 5 + 2 )
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, start=c(rep(3,25),rep(1,25)))
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0), c(0, 0, 0), c(0, exp(-1), 0)), analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0), c(0, 0, 0), c(0, exp(-1), 0)), analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

                                        #1,2,4 (= 3,4,5)
    (three.q <- fixq(rbind(c(0, 0, 0), c(exp(-3), 0, exp(-3)), c(exp(-3), 0, 0))))
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, start=rep(2, 50))
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, 0, 0), c(exp(-2), 0, exp(-2)), c(2*exp(-2), 0, 0)), analyticp=TRUE))

    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, 0, 0), c(exp(-2), 0, exp(-2)), c(2*exp(-2), 0, 0)), analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

                                        # 1,2,4
    (three.q <- fixq(rbind(c(0, exp(-3), exp(-6)), c(0, 0, exp(-3)), c(0, 0, 0))))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, exp(-1)), c(0, 0, 0)), analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, exp(-1)), c(0, 0, 0)), analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

                                        # 1,2,4 (=1,2,6)
    (three.q <- fixq(rbind(c(0, exp(-3), exp(-6)), c(0, 0, 0), c(0, exp(-3), 0))))
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, 0), c(0, exp(-2), 0)), analyticp=TRUE))


    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, 0), c(0, exp(-2), 0)), analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

                                        #1,3,5  (= 2,4,5)
    (three.q <- fixq(rbind(c(0, 0, exp(-3)), c(0, 0, exp(-3)), c(exp(-3), 0, 0))))
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, start=c(rep(1, 25), rep(2,25)))
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, 0, exp(-1)), c(0, 0, exp(-1)), c(exp(-2), 0, 0)), fixedpars=TRUE, analyticp=TRUE))

    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, 0, exp(-1)), c(0, 0, exp(-1)), c(exp(-2), 0, 0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

                                        #1,2,4,6
    (three.q <- fixq(rbind(c(0, exp(-3), exp(-3)), c(0, 0, exp(-3)), c(0, exp(-3), 0))))
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, start=c(rep(1, 25), rep(2,25)))
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-1)), c(0, 0, exp(-1)), c(0, exp(-2), 0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-1)), c(0, 0, exp(-1)), c(0, exp(-2), 0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-1)), c(0, 0, exp(-1)), c(0, exp(-1), 0)), fixedpars=TRUE, analyticp=TRUE))

    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), exp(-1)), c(0, 0, exp(-1)), c(0, exp(-1), 0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
})

## FOUR STATES
                                        #1,5,9
test_that("4 state analytic P matrices",{
    (four.q <- fixq(rbind(c(0, exp(-3), 0, 0), c(0, 0, exp(-3), 0), c(0, 0, 0, exp(-3)), c(0,0,0,0))))
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=four.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0, 0), c(0, 0, exp(-1), 0), c(0, 0, 0, exp(-1)), c(0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0, 0), c(0, 0, exp(-1), 0), c(0, 0, 0, exp(-1)), c(0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0, 0), c(0, 0, exp(-1), 0), c(0, 0, 0, exp(-2)), c(0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0, 0), c(0, 0, exp(-1), 0), c(0, 0, 0, exp(-2)), c(0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0, 0), c(0, 0, exp(-2), 0), c(0, 0, 0, exp(-1)), c(0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-1), 0, 0), c(0, 0, exp(-2), 0), c(0, 0, 0, exp(-1)), c(0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, 0), c(0, 0, exp(-1), 0), c(0, 0, 0, exp(-1)), c(0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, 0), c(0, 0, exp(-1), 0), c(0, 0, 0, exp(-1)), c(0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, 0), c(0, 0, exp(-1), 0), c(0, 0, 0, exp(-3)), c(0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, 0), c(0, 0, exp(-1), 0), c(0, 0, 0, exp(-3)), c(0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

                                        #13569
    (four.q <- fixq(rbind(c(0, exp(-3), 0, exp(-3)), c(0, 0, exp(-3), exp(-3)), c(0, 0, 0, exp(-3)), c(0,0,0,0))))
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=four.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, exp(-2)), c(0, 0, exp(-2), exp(-2)), c(0, 0, 0, 2*exp(-2)), c(0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, exp(-2)), c(0, 0, exp(-2), exp(-2)), c(0, 0, 0, 2*exp(-2)), c(0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik, tol=1e-06)
    ## no convergence with analytic, in 1.2 or 1.3, but matches with fixed
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, exp(-2)), c(0, 0, exp(-2), exp(-2)), c(0, 0, 0, exp(-2)), c(0,0,0,0)), fixedpars=FALSE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, exp(-2)), c(0, 0, exp(-2), exp(-2)), c(0, 0, 0, exp(-2)), c(0,0,0,0)), fixedpars=FALSE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, exp(-2)), c(0, 0, exp(-3), exp(-3)), c(0, 0, 0, 2*exp(-2)), c(0,0,0,0)), fixedpars=FALSE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, exp(-2)), c(0, 0, exp(-3), exp(-3)), c(0, 0, 0, 2*exp(-2)), c(0,0,0,0)), fixedpars=FALSE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-3), 0, exp(-3)), c(0, 0, exp(-2), exp(-2)), c(0, 0, 0, 2*exp(-2)), c(0,0,0,0)), fixedpars=FALSE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-3), 0, exp(-3)), c(0, 0, exp(-2), exp(-2)), c(0, 0, 0, 2*exp(-2)), c(0,0,0,0)), fixedpars=FALSE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, exp(-2)), c(0, 0, exp(-3), exp(-3)), c(0, 0, 0, exp(-2)), c(0,0,0,0)), fixedpars=FALSE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0, exp(-2), 0, exp(-2)), c(0, 0, exp(-3), exp(-3)), c(0, 0, 0, exp(-2)), c(0,0,0,0)), fixedpars=FALSE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
})

## FIVE STATES
                                        #1_6_11_16
test_that("5 state analytic P matrices",{
    (five.q <- fixq(rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0))))
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=five.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),0,0), c(0,0,0,exp(-4),0), c(0,0,0,0,exp(-5)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),0,0), c(0,0,0,exp(-4),0), c(0,0,0,0,exp(-5)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-3),0), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-3),0), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),0,0), c(0,0,0,exp(-4),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),0,0), c(0,0,0,exp(-4),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-4),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-4),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-4),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-4),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-3)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-3)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-3),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-3),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),0,0), c(0,0,0,exp(-2),0), c(0,0,0,0,exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

                                        #1_4_6_8_11_12_16
    (five.q <- fixq(rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,exp(-2)), c(0,0,0,0,0))))
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=five.q)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-3),0,exp(-3)), c(0,0,0,exp(-4),exp(-4)), c(0,0,0,0,exp(-5)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-3),0,exp(-3)), c(0,0,0,exp(-4),exp(-4)), c(0,0,0,0,exp(-5)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-3),exp(-3)), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-3),exp(-3)), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-3),0,exp(-3)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-3),0,exp(-3)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-3),0,exp(-3)), c(0,0,0,exp(-4),exp(-4)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-3),0,exp(-3)), c(0,0,0,exp(-4),exp(-4)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,exp(-3)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,exp(-3)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,exp(-4)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,exp(-3)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-4),exp(-4)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,exp(-3)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-4),exp(-4)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,exp(-3)), c(0,0,exp(-4),0,exp(-4)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,exp(-3)), c(0,0,exp(-4),0,exp(-4)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,exp(-3)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,exp(-3)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-3),exp(-3)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-3),exp(-3)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-3),0,exp(-3)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-3),0,exp(-3)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,exp(-3)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,exp(-3)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,exp(-2)), c(0,0,exp(-2),0,exp(-2)), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,2*exp(-2)), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

                                        #1_6_7_11_12
    (five.q <- fixq(rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),exp(-2),0), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,0), c(0,0,0,0,0))))
    set.seed(22061976)
    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=five.q)

    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),exp(-3),0), c(0,0,0,exp(-4),exp(-4)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-3),exp(-3),0), c(0,0,0,exp(-4),exp(-4)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,2*exp(-2),0,0,0), c(0,0,exp(-2),exp(-2),0), c(0,0,0,exp(-3),exp(-3)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,2*exp(-2),0,0,0), c(0,0,exp(-2),exp(-2),0), c(0,0,0,exp(-3),exp(-3)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,2*exp(-2),0,0,0), c(0,0,exp(-3),exp(-3),0), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,2*exp(-2),0,0,0), c(0,0,exp(-3),exp(-3),0), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-2),exp(-2),0), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-3),0,0,0), c(0,0,exp(-2),exp(-2),0), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)

    (sim.mod1 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),exp(-2),0), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=TRUE))
    (sim.mod2 <- msm(state ~ time, subject=subject, data=sim2.df, qmatrix = rbind(c(0,exp(-2),0,0,0), c(0,0,exp(-2),exp(-2),0), c(0,0,0,exp(-2),exp(-2)), c(0,0,0,0,0), c(0,0,0,0,0)), fixedpars=TRUE, analyticp=FALSE))
    expect_equal(sim.mod1$minus2loglik, sim.mod2$minus2loglik)
})
