test_that("msm2Surv example works",{
  msmdat <- data.frame(
    subj = c(1, 1, 1, 1, 1, 2, 2, 2),
    days = c(0, 27, 75, 97, 1106, 0, 90, 1037),
    status = c(1, 2, 3, 4, 4, 1, 2, 2),
    age = c(66, 66, 66, 66, 69, 49, 49, 51),
    treat = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
  # transitions only allowed to next state up or state 4
  Q <- rbind(c(1, 1, 0, 1), 
             c(0, 1, 1, 1),
             c(0, 0, 1, 1),
             c(0, 0, 0, 0))
  dat <- msm2Surv(data=msmdat, subject="subj", time="days", state="status", 
                  Q=Q)
  expect_true(inherits(dat, "msdata"))
  expect_equivalent(na.omit(as.vector(t(attr(dat, "trans")))),
                    1:5)
})
