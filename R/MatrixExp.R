#' Matrix exponential
#' 
#' Calculates the exponential of a square matrix.
#' 
#' See the \code{\link[expm]{expm}} documentation for details of the algorithms
#' it uses.
#' 
#' Generally the exponential \eqn{E} of a square matrix \eqn{M} can often be
#' calculated as
#' 
#' \deqn{E = U \exp(D) U^{-1}}{E = U exp(D) U^{-1}}
#' 
#' where \eqn{D} is a diagonal matrix with the eigenvalues of \eqn{M} on the
#' diagonal, \eqn{\exp(D)}{exp(D)} is a diagonal matrix with the exponentiated
#' eigenvalues of \eqn{M} on the diagonal, and \eqn{U} is a matrix whose
#' columns are the eigenvectors of \eqn{M}.
#' 
#' This method of calculation is used if \code{"pade"} or \code{"series"} is
#' supplied but \eqn{M} has distinct eigenvalues.  I If \eqn{M} has repeated
#' eigenvalues, then its eigenvector matrix may be non-invertible. In this
#' case, the matrix exponential is calculated using the Pade approximation
#' defined by Moler and van Loan (2003), or the less robust power series
#' approximation,
#' 
#' \deqn{\exp(M) = I + M + M^2/2 + M^3 / 3! + M^4 / 4! + ...}{exp(M) = I + M +
#' M^2/2 + M^3 / 3! + M^4 / 4! + ...}
#' 
#' For a continuous-time homogeneous Markov process with transition intensity
#' matrix \eqn{Q}, the probability of occupying state \eqn{s} at time \eqn{u +
#' t} conditional on occupying state \eqn{r} at time \eqn{u} is given by the
#' \eqn{(r,s)} entry of the matrix \eqn{\exp(tQ)}{exp(tQ)}.
#' 
#' If \code{mat} is a valid transition intensity matrix for a continuous-time
#' Markov model (i.e. diagonal entries non-positive, off-diagonal entries
#' non-negative, rows sum to zero), then for certain simpler model structures,
#' there are analytic formulae for the individual entries of the exponential of
#' \code{mat}.  These structures are listed in the PDF manual and the formulae
#' are coded in the \pkg{msm} source file \code{src/analyticp.c}.  These
#' formulae are only used if \code{method="analytic"}.  This is more efficient,
#' but it is not the default in \code{MatrixExp} because the code is not robust
#' to extreme values.  However it is the default when calculating likelihoods
#' for models fitted by \code{\link{msm}}.
#' 
#' The implementation of the Pade approximation used by \code{method="pade"}
#' was taken from JAGS by Martyn Plummer
#' (\url{https://mcmc-jags.sourceforge.io}).
#' 
#' @param mat A square matrix
#' 
#' @param t An optional scaling factor for \code{mat}.  If a vector is supplied,
#' then an array of matrices is returned with different scaling factors. 
#' 
#' @param method Under the default of \code{NULL}, this simply wraps the
#' \code{\link[expm]{expm}} function from the \pkg{expm} package.  This is
#' recommended.  Options to \code{\link{expm}} can be supplied to
#' \code{\link{MatrixExp}}, including \code{method}.
#' 
#' Otherwise, for backwards compatibility, the following options, which use
#' code in the \pkg{msm} package, are available: \code{"pade"} for a Pade
#' approximation method, \code{"series"} for the power series approximation, or
#' \code{"analytic"} for the analytic formulae for simpler Markov model
#' intensity matrices (see below).  These options are only used if \code{mat}
#' has repeated eigenvalues, thus the usual eigen-decomposition method cannot
#' be used.
#' @param ... Arguments to pass to \code{\link{expm}}.
#' @return The exponentiated matrix \eqn{\exp(mat)}{exp(mat)}. Or, if \code{t}
#' is a vector of length 2 or more, an array of exponentiated matrices.
#' @references Cox, D. R. and Miller, H. D. \emph{The theory of stochastic
#' processes}, Chapman and Hall, London (1965)
#' 
#' Moler, C and van Loan, C (2003).  Nineteen dubious ways to compute the
#' exponential of a matrix, twenty-five years later.  \emph{SIAM Review}
#' \bold{45}, 3--49.
#' @keywords math
#' @export MatrixExp
MatrixExp <- function(mat, t = 1, method=NULL,...){
    if (!is.matrix(mat) || (nrow(mat)!= ncol(mat)))
        stop("\"mat\" must be a square matrix")
    qmodel <- if (is.qmatrix(mat) && !is.null(method) && method=="analytic") msm.form.qmodel(mat) else list(iso=0, perm=0, qperm=0)
    if (!is.null(method) && method=="analytic") {
        if (!is.qmatrix(mat))
            warning("Analytic method not available since matrix is not a Markov model intensity matrix. Using \"pade\".")
        else if (qmodel$iso==0) warning("Analytic method not available for this Markov model structure. Using \"pade\".")

    }
    if (length(t) > 1) res <- array(dim=c(dim(mat), length(t)))
    for (i in seq_along(t)) {
        if (is.null(method) || !(method %in% c("pade","series","analytic"))) {
            if (is.null(method)) method <- eval(formals(expm::expm)$method)
            resi <- expm::expm(t[i]*mat, method=method, ...)
        } else {
            ccall <- .C("MatrixExpR", as.double(mat), as.integer(nrow(mat)), res=double(length(mat)), as.double(t[i]),
                        as.integer(match(method, c("pade","series"))), # must match macro constants in pijt.c
                        as.integer(qmodel$iso), as.integer(qmodel$perm), as.integer(qmodel$qperm),
                        as.integer(0), NAOK=TRUE)
            resi <- matrix(ccall$res, nrow=nrow(mat))
        }
        if (length(t)==1) res <- resi
        else res[,,i] <- resi
    }
    res
}
