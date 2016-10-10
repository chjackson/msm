twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
twoway4.q2 <- rbind(c(-0.51, 0.24, 0, 0.25), c(0.162, -0.498, 0.168, 0.166), c(0, 0.26, -0.5, 0.25), c(0, 0, 0, 0))
twoway3.q <- rbind(c(-0.5, 0.25, 0), c(0.166, -0.498, 0.166), c(0, 0.25, -0.5))
oneway4.q <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), c(0, 0, 0, 0.126), c(0, 0, 0, 0))
rownames(twoway4.q) <- colnames(twoway4.q) <- rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well","Mild","Severe","Death")
twoway4.i <- twoway4.q; twoway4.i[twoway4.i!=0] <- 1
oneway4.i <- oneway4.q; oneway4.i[oneway4.i!=0] <- 1
psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
fiveq <- rbind(c(0,0.01,0,0,0.002), c(0,0,0.07,0,0.01), c(0,0,0,0.07,0.02), c(0,0,0,0,0.03), c(0,0,0,0,0))
ematrix <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.1, 0),c(0, 0.1, 0, 0),c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- rownames(ematrix) <- colnames(ematrix) <- c("Well","Mild","Severe","Death")

cav.cens <- cav
cav.cens$state[cav$state==4][1:50] <- 99
cav.cens2 <- cav
cav.cens2$state[cav$state==4][1:50] <- 99
cav.cens2$state[cav$state==4][51:100] <- 999
cav.cens3 <- cav
ns <- c(cav$state[2:nrow(cav)], 0)
cav.cens3$state[cav$state==4][1:50] <- 99
cav.cens3$state[ns==4][1:50] <- 999

deriv_error <- function(object){
    if (!isTRUE(getOption("msm.test.analytic.derivatives")))
        stop("msm.test.analytic.derivatives option not set")
    object$paramdata$deriv.test$error["nd"]    
}
