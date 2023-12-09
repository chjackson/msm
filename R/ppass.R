#' Passage probabilities
#' 
#' Probabilities of having visited each state by a particular time in a
#' continuous time Markov model.
#' 
#' The passage probabilities to state \eqn{s} are computed by setting the
#' \eqn{s}th row of the transition intensity matrix \eqn{Q} to zero, giving an
#' intensity matrix \eqn{Q^*}{Q*} for a simplified model structure where state
#' \eqn{s} is absorbing.  The probabilities of passage are then equivalent to
#' row \eqn{s} of the transition probability matrix \eqn{Exp(tQ^*)}{Exp(tQ*)} 
#' (\code{\link{pmatrix.msm}}) under this
#' simplified model for \eqn{t=}\code{tot}.   
#' 
#' For time-inhomogenous models, 
#' this process is generalised by calculating an intensity matrix for each
#' time period, zeroing the appropriate row of each, and calculating and multiplying
#' transition probability matrices as in \code{\link{pmatrix.piecewise.msm}}.
#' 
#' Note this is different from the probability of occupying each state at
#' exactly time \eqn{t}, given by \code{\link{pmatrix.msm}}.  The passage
#' probability allows for the possibility of having visited the state before
#' \eqn{t}, but then occupying a different state at \eqn{t}.
#' 
#' The mean of the passage distribution is the expected first passage time,
#' \code{\link{efpt.msm}}.
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#'
#' @param qmatrix Instead of \code{x}, you can simply supply a transition
#' intensity matrix in \code{qmatrix}.
#'
#' @param tot Finite time to forecast the passage probabilites for.
#'
#' @param start Starting state (integer).  By default (\code{start="all"}),
#' this will return a matrix one row for each starting state.
#' 
#' Alternatively, this can be used to obtain passage probabilities from a
#' \emph{set} of states, rather than single states.  To achieve this,
#' \code{state} is set to a vector of weights, with length equal to the number
#' of states in the model.  These weights should be proportional to the
#' probability of starting in each of the states in the desired set, so that
#' weights of zero are supplied for other states.  The function will calculate
#' the weighted average of the passage probabilities from each of the
#' corresponding states.
#' 
#' @param covariates Covariate values defining the intensity matrix for the
#' fitted model \code{x}, as supplied to \code{\link{qmatrix.msm}}.
#'
#' @param piecewise.times For models with time-dependent covariates,
#'   this defines the cut points in time at which the transition
#'   intensity matrix changes.  This is not required for models fitted
#'   with the \code{pci} option to \code{\link{msm}}, which are
#'   handled automatically.
#'
#' @param piecewise.covariates For models with time-dependent
#'   covariates, this is the list of covariates for each time period
#'   defined by \code{piecewise.times}, in the format documented for
#'   the \code{covariates} argument to
#'   \code{\link{pmatrix.piecewise.msm}}.  This is not required for
#'   models fitted with the \code{pci} option to \code{\link{msm}},
#'   which are handled automatically.
#'
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
#'
#' @param cl Width of the symmetric confidence interval, relative to 1.
#'
#' @param B Number of bootstrap replicates.
#'
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#'
#' @param ... Arguments to pass to \code{\link{MatrixExp}}.
#'
#' @return A matrix whose \eqn{r, s} entry is the probability of having visited
#' state \eqn{s} at least once before time \eqn{t}, given the state at time
#' \eqn{0} is \eqn{r}.  The diagonal entries should all be 1.
#'
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#'
#' @seealso \code{\link{efpt.msm}}, \code{\link{totlos.msm}},
#' \code{\link{boot.msm}}.
#'
#' @references Norris, J. R. (1997) Markov Chains. Cambridge University Press.
#'
#' @keywords models
#'
#' @examples
#' 
#' Q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
#'            c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
#' 
#' ## ppass[1,2](t) converges to 0.5 with t, since given in state 1, the
#' ## probability of going to the absorbing state 4 before visiting state
#' ## 2 is 0.5, and the chance of still being in state 1 at t decreases.
#' 
#' ppass.msm(qmatrix=Q, tot=2)
#' ppass.msm(qmatrix=Q, tot=20)
#' ppass.msm(qmatrix=Q, tot=100)
#' 
#' Q <- Q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
#' 
#' ## Probability of about 1/2 of visiting state 3 by time 10.5, the
#' ## median first passage time
#' 
#' ppass.msm(qmatrix=Q, tot=10.5)
#' 
#' ## Mean first passage time from state 2 to state 3 is 10.02: similar
#' ## to the median
#' 
#' efpt.msm(qmatrix=Q, tostate=3)
#' 
#' @export ppass.msm
ppass.msm <- function(x=NULL, qmatrix=NULL, tot, start="all", covariates="mean",
                      piecewise.times=NULL, piecewise.covariates=NULL,
                      ci=c("none","normal","bootstrap"), cl = 0.95, B = 1000, cores=NULL, ...)
{
  ci <- match.arg(ci)
  if (!is.null(x) && !inherits(x, "msm"))
    stop("expected x to be a msm model")
  if (is.null(x$pci) && is.null(piecewise.times)){
    if (!is.null(x)) {
      qmatrix <- qmatrix.msm(x, covariates=covariates, ci="none")
    }
    else if (!is.null(qmatrix)) {
      if (!is.matrix(qmatrix) || (nrow(qmatrix) != ncol(qmatrix)))
        stop("expected qmatrix to be a square matrix")
      if (ci != "none") {warning("No fitted model supplied: not calculating confidence intervals."); ci <- "none"}
    }
    res <- array(dim=dim(qmatrix))
    if (!is.null(dimnames(qmatrix))) {
      dimnames(res) <- dimnames(qmatrix)
      names(dimnames(res)) <- c("from","to")
    }
    states <- 1:nrow(qmatrix)
    for (i in states) {
      Qred <- qmatrix; Qred[states[i],] <- 0
      res[,i] <- MatrixExp(Qred*tot, ...)[,i]
    }
    if (!is.character(start)) {
      if (!is.numeric(start) || (length(start)!=nrow(qmatrix)))
        stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmatrix))
      start <- start / sum(start)
      res <- matrix(start %*% res, nrow=1, dimnames=list(from="start",to=colnames(res)))
    }
    else if (any(start!="all"))
      stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmatrix))
  } else {
    res <- ppass.td.msm(x=x, tot=tot, start=start, covariates=covariates,
                         piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates, ...)
  }
  
  p.ci <- switch(ci,
                 bootstrap = ppass.ci.msm(x=x, qmatrix=qmatrix, tot=tot, start=start,
                                          covariates=covariates,
                                          piecewise.times=piecewise.times,
                                          piecewise.covariates=piecewise.covariates,
                                          cl=cl, B=B, cores=cores),
                 normal = ppass.normci.msm(x=x, qmatrix=qmatrix, tot=tot, start=start,
                                           covariates=covariates, 
                                           piecewise.times=piecewise.times,
                                           piecewise.covariates=piecewise.covariates,
                                           cl=cl, B=B),
                 none = NULL)
  if (ci != "none") {
    res <- list(estimates=res, L=p.ci$L, U=p.ci$U)
  }
  class(res) <- "msm.est"
  res
}


## Thanks to Jon Fintzi <https://github.com/fintzij>

ppass.td.msm <- function(x=NULL, tot, start="all", covariates="mean",
                          piecewise.times=NULL, piecewise.covariates=NULL, ...) {
  
  if (!is.null(x$pci)){
    ## this goes to the times argument in pmatrix.piecewise.msm
    piecewise.times <- x$pci 
    
    ## for defining piecewise constant intensity matrices
    piecewise.covariates <- msm.fill.pci.covs(x, covariates)
  } else {
    ## user-supplied times and covariates for non-pci time-dependent models
    validate_piecewise(piecewise.times, piecewise.covariates)
  }
  
  ## get temporary qmatrix
  qmat_temp <- qmatrix.msm(x, covariates = "mean", ci = "none")
  
  ## array for storing the first passage probabilities
  res <- array(dim=dim(qmat_temp))
  if (!is.null(dimnames(qmat_temp))) {
    dimnames(res) <- dimnames(qmat_temp)
    names(dimnames(res)) <- c("from","to")
  }
  
  ## generate list of intensity matrices
  Q_list = lapply(piecewise.covariates,
                  function(y) qmatrix.msm(x, covariates = y, ci = "none"))
  states <- seq_len(nrow(qmat_temp))
  for (i in states) {
    
    ## zero out the appropriate row in the intensity matrices
    Q_ppass <- 
      lapply(Q_list, 
             function(M) {Q = M; Q[i,] = 0; return(Q)})
    
    ## supply Q_ppass to pmatrix.piecewise.msm
    res[,i] <-
      pmatrix.piecewise.msm(qlist = Q_ppass,
                            times = piecewise.times, 
                            covariates = piecewise.covariates,
                            t1 = 0, t2 = tot, 
                            ci = "none")[,i]
  }
  
  if (!is.character(start)) {
    if (!is.numeric(start) || (length(start)!=nrow(qmat_temp)))
      stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmat_temp))
    start <- start / sum(start)
    res <- matrix(start %*% res, nrow=1, dimnames=list(from="start",to=colnames(res)))
  }
  else if (any(start!="all"))
    stop("Expected \"start\" to be \"all\" or a numeric vector of length ", nrow(qmat_temp))
  
  res
}
