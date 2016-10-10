context("msm data summaries")

test_that("statetable.msm", {
    stab <- statetable.msm(state, PTNUM, data=cav)
    expect_that(stab, equals(structure(c(1367L, 46L, 4L, 204L, 134L, 13L, 44L, 54L, 107L, 148L, 48L, 55L), .Dim = 3:4, .Dimnames = structure(list(from = c("1","2", "3"), to = c("1", "2", "3", "4")), .Names = c("from", "to")), class = "table")))
    expect_equal(as.numeric(stab), c(1367, 46, 4, 204, 134, 13, 44, 54, 107, 148, 48, 55))
    expect_error(statetable.msm(state,PTNUM), "not found")
    stabc <- statetable.msm(state, PTNUM, cav.cens)
    expect_equal(as.numeric(stabc), c(1367, 46, 4, 204, 134, 13, 44, 54, 107, 127, 40, 34, 21, 8, 21))
})

test_that("crudeinits.msm", {
    cinits <- crudeinits.msm(state ~ years, PTNUM, data=cav, qmatrix=twoway4.q)
    expect_equal(as.numeric(cinits), c(-0.117314905786477, 0.116817878212849, 0, 0, 0.067989320398981, -0.375848825554382, 0.049084006577444, 0, 0, 0.137134030945518, -0.256747111328168, 0, 0.049325585387496, 0.121896916396016, 0.207663104750724, 0))
    expect_error(crudeinits.msm(state ~ years, PTNUM, qmatrix=twoway4.q), "not found")
})

test_that("crudeinits.msm ignores inconsistent transitions unless exact times", {
    cav.wrong <- cav; cav.wrong$state[4] <- 1
    expect_error(crudeinits.msm(state ~ years, PTNUM, oneway4.q, cav), NA)
    expect_warning(msm(state ~ years, subject=PTNUM, data = cav.wrong, qmatrix = twoway4.q, exacttimes=TRUE, fixedpars=TRUE), "inconsistent with intensity")
    obstype <- rep(2, nrow(cav)); obstype[10] <- 1
    expect_warning(msm(state ~ years, subject=PTNUM, data = cav.wrong, qmatrix = twoway4.q, obstype=obstype, fixedpars=TRUE), "inconsistent with intensity")
})

test_that("crudeinits.msm handles NAs", {
    psor2 <- psor;  psor2$ptnum[13:14] <- psor2$months[7:8] <- psor2$state[7:8] <- NA
    expect_equal(crudeinits.msm(state ~ months, ptnum, data=psor2, qmatrix=psor.q),
                 crudeinits.msm(state ~ months, ptnum, data=psor2[-c(7,8,13,14),], qmatrix=psor.q))
})

test_that("crudeinits.msm handles censoring",{
    expect_error(crudeinits.msm(state ~ years, PTNUM, twoway4.q, cav.cens), "elements not in 1, 2")
    cru <- crudeinits.msm(state~ years, PTNUM, twoway4.q, cav.cens, censor=99)
    expect_equal(as.numeric(cru), c(-0.111798558088660, 0.122878533946307, 0, 0, 0.0689030388220138, -0.373978146793108, 0.0618112185064827, 0, 0, 0.144248713763056, -0.223471328446514, 0, 0.0428955192666458, 0.106850899083745, 0.161660109940032, 0))
    cru <- crudeinits.msm(state~ years, PTNUM, twoway4.q, cav.cens2, censor=c(99,999), censor.states=list(c(1,2,3),c(2,3)))
    expect_equal(as.numeric(cru), c(-0.107299472349819, 0.134927714425074, 0, 0, 0.0697104852209013, -0.369584609077378, 0.0789635719132074, 0, 0, 0.158393403890305, -0.170075385659216, 0, 0.0375889871289174, 0.0762634907619986, 0.0911118137460085, 0), tol=1e-06)
    cru <- crudeinits.msm(state~ years, PTNUM, twoway4.q, cav.cens3, censor=c(99,999), censor.states=list(c(2,3),c(1,2,3)))
    expect_equal(as.numeric(cru), c(-0.112107245394208, 0.124370575094641, 0, 0, 0.0686998668421821, -0.370347934726264, 0.0659650781282531, 0, 0, 0.135425737325276, -0.238489128617530, 0, 0.0434073785520255, 0.110551622306348, 0.172524050489277, 0), tol=1e-06)
})
