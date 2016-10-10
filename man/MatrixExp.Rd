\name{MatrixExp}
\title{Matrix exponential}
\alias{MatrixExp}
\concept{Matrix exponential}
\description{
  Calculates the exponential of a square matrix.
}
\usage{
MatrixExp(mat, t = 1, method=NULL, ...)
}
\arguments{
  \item{mat}{A square matrix}
  \item{t}{An optional scaling factor for \code{mat}.}
  \item{method}{ Under the default of \code{NULL}, this simply wraps
    the \code{\link[expm]{expm}} function from the \pkg{expm} package.
    This is recommended.  Options to \code{\link{expm}} can be supplied to
    \code{\link{MatrixExp}}, including \code{method}.

    Otherwise, for backwards compatibility, the following options, which
    use code in the \pkg{msm} package, are available:  
    \code{"pade"} for a Pade approximation method, \code{"series"} for the
    power series approximation, or \code{"analytic"} for the analytic
    formulae for simpler Markov model intensity matrices (see below).
    These options are only used if \code{mat} has repeated
    eigenvalues, thus the usual eigen-decomposition method cannot be used.
  }
  \item{...}{Arguments to pass to \code{\link{expm}}.}
}
\value{
  The exponentiated matrix \eqn{\exp(mat)}{exp(mat)}. Or, if \code{t} is
  a vector of length 2 or more, an array of exponentiated matrices.
}
\details{
  See the \code{\link[expm]{expm}} documentation for details of the
    algorithms it uses. 

  Generally the exponential \eqn{E} of a square matrix \eqn{M} can
  often be calculated as

  \deqn{E = U \exp(D) U^{-1}}{E = U exp(D) U^{-1}}

  where \eqn{D} is a diagonal matrix with the eigenvalues
  of \eqn{M} on the diagonal, \eqn{\exp(D)}{exp(D)} is a diagonal matrix with the
  exponentiated eigenvalues
  of \eqn{M} on the diagonal, and \eqn{U} is a matrix whose columns are the
  eigenvectors of \eqn{M}.

  This method of calculation is used if \code{"pade"} or \code{"series"}
  is supplied but \eqn{M} has distinct
  eigenvalues.  I If \eqn{M} has repeated eigenvalues, then its
  eigenvector matrix may be non-invertible. In this case, the matrix
  exponential is calculated using the Pade approximation defined by
  Moler and van Loan (2003), or the less robust power series
  approximation,

  \deqn{\exp(M) = I + M + M^2/2 + M^3 / 3! + M^4 / 4! + ...}{exp(M) = I + M + M^2/2 + M^3 / 3! + M^4 / 4! + ...}

  For a continuous-time homogeneous Markov process with transition
intensity matrix
  \eqn{Q},  the probability of occupying state \eqn{s} at time \eqn{u + t}
  conditional on occupying state \eqn{r} at time \eqn{u} is given by the
  \eqn{(r,s)} entry of the matrix \eqn{\exp(tQ)}{exp(tQ)}.

  If \code{mat} is a valid transition intensity matrix for a
  continuous-time Markov model (i.e. diagonal entries non-positive,
  off-diagonal entries non-negative, rows sum to zero), then for certain
  simpler model structures, there are
  analytic formulae for the individual entries of the exponential
  of \code{mat}.  These structures are
  listed in the PDF manual and the formulae are coded in the \pkg{msm}
  source file \code{src/analyticp.c}.   These formulae are only used if
  \code{method="analytic"}.  This is more efficient, but it is not the
  default in \code{MatrixExp} because the code is not robust to extreme
  values.  However it is the default when calculating likelihoods for models fitted
  by \code{\link{msm}}.

  The implementation of the Pade approximation used by
  \code{method="pade"} was taken from
  JAGS by Martyn Plummer (\url{http://mcmc-jags.sourceforge.net}).
}
\references{
  Cox, D. R. and Miller, H. D. \emph{The theory of stochastic
    processes}, Chapman and Hall, London (1965)

  Moler, C and van Loan, C (2003).  Nineteen dubious ways to compute
  the exponential of a matrix, twenty-five years later.
  \emph{SIAM Review} \bold{45}, 3--49.
}
\keyword{math}
