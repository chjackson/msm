# Could also get CDF simply by making tostate absorbing and calculating pmatrix.
# TODO time-dependent covariates?.  Unclear if the expectation has a solution for piecewise-constant rate.

#' Expected first passage time
#' 
#' Expected time until first reaching a particular state or set of states in a
#' Markov model.
#' 
#' The expected first passage times from each of a set of states
#' \eqn{\mathbf{i}}{i} to to the remaining set of states
#' \eqn{\overline{\mathbf{i}}}{ibar} in the state space, for a model with
#' transition intensity matrix \eqn{Q}, are
#' 
#' \deqn{-Q_{\mathbf{i},\mathbf{i}}^{-1} \mathbf{1}}{-Q_{i,i}^{-1} 1}
#' 
#' where \eqn{\mathbf{1}}{1} is a vector of ones, and
#' \eqn{Q_{\mathbf{i},\mathbf{i}}}{Q_{i,i}} is the square subset of \eqn{Q}
#' pertaining to states \eqn{\mathbf{i}}{i}.
#' 
#' It is equal to the sum of mean sojourn times for all states between the
#' "from" and "to" states in a unidirectional model.  If there is non-zero
#' chance of reaching an absorbing state before reaching \code{tostate}, then
#' it is infinite.  It is trivially zero if the "from" state equals
#' \code{tostate}.
#' 
#' This function currently only handles time-homogeneous Markov models.  For
#' time-inhomogeneous models it will assume that \eqn{Q} equals the average
#' intensity matrix over all times and observed covariates.  Simulation might
#' be used to handle time dependence.
#' 
#' Note this is the \emph{expectation} of first passage time, and the
#' confidence intervals are CIs for this mean, not predictive intervals for the
#' first passage time.  The full distribution of the first passage time to a
#' set of states can be obtained by setting the rows of the intensity matrix
#' \eqn{Q} corresponding to that set of states to zero to make a model where
#' those states are absorbing.  The corresponding transition probability matrix
#' \eqn{Exp(Qt)} then gives the probabilities of having hit or passed that
#' state by a time \eqn{t} (see the example below). This is implemented in
#' \code{\link{ppass.msm}}.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param qmatrix Instead of \code{x}, you can simply supply a transition
#' intensity matrix in \code{qmatrix}.
#' @param tostate State, or set of states supplied as a vector, for which to
#' estimate the first passage time into.  Can be integer, or character matched
#' to the row names of the Q matrix.
#' @param start Starting state (integer).  By default (\code{start="all"}),
#' this will return a vector of expected passage times from each state in turn.
#' 
#' Alternatively, this can be used to obtain the expected first passage time
#' from a \emph{set} of states, rather than single states.  To achieve this,
#' \code{state} is set to a vector of weights, with length equal to the number
#' of states in the model.  These weights should be proportional to the
#' probability of starting in each of the states in the desired set, so that
#' weights of zero are supplied for other states.  The function will calculate
#' the weighted average of the expected passage times from each of the
#' corresponding states.
#' @param covariates Covariate values defining the intensity matrix for the
#' fitted model \code{x}, as supplied to \code{\link{qmatrix.msm}}.
#' @param ci If \code{"normal"}, then calculate a confidence interval by
#' simulating \code{B} random vectors from the asymptotic multivariate normal
#' distribution implied by the maximum likelihood estimates (and covariance
#' matrix) of the log transition intensities and covariate effects.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' 
#' If \code{"none"} (the default) then no confidence interval is calculated.
#' @param cl Width of the symmetric confidence interval, relative to 1.
#' @param B Number of bootstrap replicates.
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @param ... Arguments to pass to \code{\link{MatrixExp}}.
#' @return A vector of expected first passage times, or "hitting times", from
#' each state to the desired state.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{sojourn.msm}}, \code{\link{totlos.msm}},
#' \code{\link{boot.msm}}.
#' @references Norris, J. R. (1997) Markov Chains. Cambridge University Press.
#' @keywords models
#' @examples
#' 
#' twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
#'              c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
#' efpt.msm(qmatrix=twoway4.q, tostate=3)
#' # given in state 1, expected time to reaching state 3 is infinite
#' # since may die (state 4) before entering state 3
#' 
#' # If we remove the death state from the model, EFPTs become finite
#' Q <- twoway4.q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
#' efpt.msm(qmatrix=Q, tostate=3)
#' 
#' # Suppose we cannot die or regress while in state 2, can only go to state 3
#' Q <- twoway4.q; Q[2,4] <- Q[2,1] <- 0; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
#' efpt.msm(qmatrix=Q, tostate=3)
#' # The expected time from 2 to 3 now equals the mean sojourn time in 2.
#' -1/Q[2,2]
#' 
#' # Calculate cumulative distribution of the first passage time
#' # into state 3 for the following three-state model
#' Q <- twoway4.q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
#' # Firstly form a model where the desired hitting state is absorbing
#' Q[3,] <- 0
#' MatrixExp(Q, t=10)[,3]
#' ppass.msm(qmatrix=Q, tot=10)
#' # Given in state 1 at time 0, P(hit 3 by time 10) = 0.479
#' MatrixExp(Q, t=50)[,3]  # P(hit 3 by time 50) = 0.98
#' ppass.msm(qmatrix=Q, tot=50)
#' 
#' 
#' @export efpt.msm
efpt.msm <- function(x=NULL, qmatrix=NULL, tostate, start="all", covariates="mean",
                     ci=c("none","normal","bootstrap"), cl = 0.95, B = 1000, cores=NULL, ...)
{
  ci <- match.arg(ci)
  if (!is.null(x)) {
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    qmatrix <- qmatrix.msm(x, covariates=covariates, ci="none")
  }
  else if (!is.null(qmatrix)) {
    if (!is.matrix(qmatrix) || (nrow(qmatrix) != ncol(qmatrix)))
      stop("expected qmatrix to be a square matrix")
    if (ci != "none") {warning("No fitted model supplied: not calculating confidence intervals."); ci <- "none"}
  }
  if (is.character(tostate)) {
    if (!tostate %in% rownames(qmatrix)) stop(sprintf("state \"%s\" unknown", tostate))
    tostate <- match(tostate, rownames(qmatrix))
  }
  est <- rep(NA, nrow(qmatrix))
  ## EFPT is zero if we're already in tostate
  est[tostate] <- 0
  abstate <- absorbing.msm(qmatrix=qmatrix)
  ## EFPT is infinite for other absorbing states
  est[setdiff(abstate,tostate)] <- Inf
  fromstate <- setdiff(1:nrow(qmatrix), union(abstate,tostate))
  ## EFPT is infinite if any chance of absorbing elsewhere before
  ## hitting tostate.  To calculate this, form Q matrix with tostate
  ## made absorbing, and look at P matrix in unit time.
  Qred <- qmatrix; Qred[tostate,] <- 0
  Pmat <- MatrixExp(Qred, ...)
  Pmat[Pmat < 1e-16] <- 0
  p.abs <- rowSums(Pmat[fromstate,setdiff(abstate,tostate),drop=FALSE])
  est[fromstate][p.abs>0] <- Inf

  ## Any states left from which EFPT to tostate is nonzero and finite.
  ## Use standard linear equation solution
  ## see, e.g. equation (3) of Harrison and Knottenbelt (2001)
  if (any(is.na(est))){
    fromstate <- which(is.na(est))
    Q <- as.matrix(qmatrix[fromstate, fromstate])
    est[fromstate] <- solve(-Q, rep(1,nrow(Q)))
  }
  if (!is.character(start)) {
    if (!is.numeric(start) || (length(start)!=nrow(qmatrix)))
      stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmatrix))
    start <- start / sum(start)
    est <- est %*% start
  }
  else if (any(start!="all"))
    stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmatrix))
  e.ci <- switch(ci,
                 bootstrap = efpt.ci.msm(x=x, qmatrix=qmatrix, tostate=tostate, start=start, covariates=covariates, cl=cl, B=B, cores=cores),
                 normal = efpt.normci.msm(x=x, qmatrix=qmatrix, tostate=tostate, start=start, covariates=covariates, cl=cl, B=B),
                 none = NULL)
  res <- if (ci=="none") est else rbind(est, e.ci)
  class(res) <- "msm.estbystate"
  res
}
