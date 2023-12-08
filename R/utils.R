### msm PACKAGE
### USEFUL FUNCTIONS NOT SPECIFIC TO MULTI-STATE MODELS

## Tests for a valid continuous-time Markov model transition intensity matrix

is.qmatrix <- function(Q) {
    Q <- unclass(Q)
    Q2 <- Q; diag(Q2) <- 0
    isTRUE(all.equal(-diag(Q), rowSums(Q2))) && isTRUE(all(diag(Q)<=0)) && isTRUE(all(Q2>=0))
}


## Transform vector of parameters constrained on [a, b] to real line.
## Vectorised.  a=-Inf or b=Inf represent unbounded below or above.

glogit <- function(x, a, b) {
    if (is.null(a)) a <- -Inf
    if (is.null(b)) b <- Inf
    ret <- numeric(length(x))
    attributes(ret) <- attributes(x)
    nn <- is.infinite(a) & is.infinite(b)
    nb <- is.infinite(a) & is.finite(b)
    an <- is.finite(a) & is.infinite(b)
    ab <- is.finite(a) & is.finite(b)
    ret[nn] <- x[nn]
    ret[nb] <- log(b[nb] - x[nb])
    ret[an] <- log(x[an] - a[an])
    ret[ab] <- log((x[ab] - a[ab]) / (b[ab] - x[ab]))
    ret
}

dglogit <- function(x, a, b) {
    if (is.null(a)) a <- -Inf
    if (is.null(b)) b <- Inf
    ret <- numeric(length(x))
    attributes(ret) <- attributes(x)
    nn <- is.infinite(a) & is.infinite(b)
    nb <- is.infinite(a) & is.finite(b)
    an <- is.finite(a) & is.infinite(b)
    ab <- is.finite(a) & is.finite(b)
    ret[nn] <- 1
    ret[nb] <- -1 / (b[nb] - x[nb])
    ret[an] <- 1 / (x[an] - a[an])
    ret[ab] <- 1/(x[ab] - a[ab]) + 1/(b[ab] - x[ab])
    ret
}

# d/dx log( (x-a)/(b-x) )    = (b-x)/(x-a) * (1/(b-x) + (x-a)/(b-x)^2)
# = 1/(x-a) + 1/(b-x)
# = d/dx   log(x-a) - log(b-x)

## Inverse transform vector of parameters constrained on [a, b]: back
## from real line to constrained scale.  Vectorised.  a=-Inf or b=Inf
## represent unbounded below or above.

gexpit <- function(x, a, b) {
    if (is.null(a)) a <- -Inf
    if (is.null(b)) b <- Inf
    ret <- numeric(length(x))
    attributes(ret) <- attributes(x)
    nn <- is.infinite(a) & is.infinite(b)
    nb <- is.infinite(a) & is.finite(b)
    an <- is.finite(a) & is.infinite(b)
    ab <- is.finite(a) & is.finite(b)
    ret[nn] <- x[nn]
    ret[nb] <- b[nb] - exp(x[nb])
    ret[an] <- exp(x[an]) + a[an]
    ret[ab] <- (b[ab]*exp(x[ab]) + a[ab]) / (1 + exp(x[ab]))
    ret
}

## Derivative of gexpit w.r.t. x

dgexpit <- function(x, a, b) {
    if (is.null(a)) a <- -Inf
    if (is.null(b)) b <- Inf
    ret <- numeric(length(x))
    attributes(ret) <- attributes(x)
    nn <- is.infinite(a) & is.infinite(b)
    nb <- is.infinite(a) & is.finite(b)
    an <- is.finite(a) & is.infinite(b)
    ab <- is.finite(a) & is.finite(b)
    ret[nn] <- 1
    ret[nb] <- - exp(x[nb])
    ret[an] <- exp(x[an])
    ret[ab] <- (b[ab] - a[ab])*exp(x[ab]) / (1 + exp(x[ab]))^2
    ret
}
