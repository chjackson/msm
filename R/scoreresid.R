#' Score residuals
#' 
#' Score residuals for detecting outlying subjects.
#' 
#' The score residual for a single subject is
#' 
#' \deqn{U(\theta)^T I(\theta)^{-1} U(\theta)}{U(theta)^T I(theta)^{-1}
#' U(theta)}
#' 
#' where \eqn{U(\theta)}{U(theta)} is the vector of first derivatives of the
#' log-likelihood for that subject at maximum likelihood estimates
#' \eqn{\theta}{theta}, and \eqn{I(\theta)}{theta} is the observed Fisher
#' information matrix, that is, the matrix of second derivatives of minus the
#' log-likelihood for that subject at theta.
#' 
#' Subjects with a higher influence on the maximum likelihood estimates will
#' have higher score residuals.
#' 
#' These are only available for models with analytic derivatives (which
#' includes all non-hidden and most hidden Markov models).
#' 
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param plot If \code{TRUE}, display a simple plot of the residuals in
#' subject order, labelled by subject identifiers
#' @return Vector of the residuals, named by subject identifiers.
#' @author Andrew Titman \email{a.titman@@lancaster.ac.uk} (theory), Chris
#' Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk} (code)
#' @keywords models
#' @export scoreresid.msm
scoreresid.msm <- function(x, plot=FALSE){
  if (!inherits(x, "msm")) stop("expected x to be a msm model")
  if (!deriv_supported(x$data, x$hmodel, x$cmodel))
    stop("Score residuals not available, since analytic derivatives not implemented for this model")
  derivs <- Ccall.msm(x$paramdata$opt$par, do.what="deriv.subj", expand.data(x), x$qmodel, x$qcmodel, x$cmodel, x$hmodel, x$paramdata)
  cov <- x$paramdata$covmat[x$paramdata$optpars,x$paramdata$optpars]
  sres <- colSums(t(derivs) * cov %*% t(derivs))
  names(sres) <- unique(x$data$mf$"(subject)")
  if (plot) {
    plot(sres, type="n")
    text(seq_along(sres), sres, names(sres))
  }
  sres
}
