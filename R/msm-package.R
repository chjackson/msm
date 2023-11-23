#' Multi-State Markov and Hidden Markov Models in Continuous Time
#'
#' msm: Functions for fitting continuous-time Markov and hidden Markov
#' multi-state models to longitudinal data.  Designed for processes
#' observed at arbitrary times in continuous time (intermittently
#' observed or panel data) but some other observation schemes are
#' supported. Both Markov transition rates and the hidden Markov
#' output process can be modelled in terms of covariates, which may be
#' constant or piecewise-constant in time.
#' 
#' @name msm-package
#' @aliases msm-package
#' @docType package
#' @useDynLib msm, .registration=TRUE
#'
#' @importFrom graphics plot persp contour image filled.contour legend lines par text
#' @importFrom grDevices rainbow 
#' @importFrom stats coef as.formula deriv dexp dnorm integrate logLik model.extract model.frame model.matrix na.fail na.omit na.pass numericDeriv optimHess pchisq pexp plogis pnorm qlogis qnorm quantile rbeta rbinom reformulate rexp rgamma rlnorm rnbinom rnorm rpois rt runif rweibull sd setNames terms uniroot
#' @importFrom utils head tail
#' @importFrom mvtnorm rmvnorm
#' @importFrom survival Surv survfit
#' @importFrom expm expm
#' 
"_PACKAGE"
