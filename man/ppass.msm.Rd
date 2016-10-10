\name{ppass.msm}
\alias{ppass.msm}
\title{Passage probabilities}
\description{
  Probabilities of having visited each state by a particular time
  in a continuous time Markov model.
}
\usage{
ppass.msm(x=NULL, qmatrix=NULL, tot, start="all", covariates="mean",
           piecewise.times=NULL, piecewise.covariates=NULL,
           ci=c("none","normal","bootstrap"), cl=0.95, B=1000,
           cores=NULL, \dots)
}
\arguments{
  \item{x}{A fitted multi-state model, as returned by
    \code{\link{msm}}. }

  \item{qmatrix}{Instead of \code{x}, you can simply supply a transition
    intensity matrix in \code{qmatrix}.}

  \item{tot}{Finite time to forecast the passage probabilites for.}

  \item{start}{Starting state (integer).  By default (\code{start="all"}),
    this will return a matrix one row for each starting state.

    Alternatively, this can be used to obtain passage probabilities
    from a \emph{set} of states, rather than single states.  To
    achieve this, \code{state} is set to a vector of weights, with
    length equal to the number of states in the model.  These weights
    should be proportional to the probability of starting in each of the
    states in the desired set, so that weights of zero are supplied for
    other states.  The function will calculate the weighted average of
    the passage probabilities from each of the corresponding states.
  }

  \item{covariates}{Covariate values defining the intensity matrix for
    the fitted model \code{x}, as supplied to \code{\link{qmatrix.msm}}.}

  \item{piecewise.times}{Currently ignored: not implemented for
    time-inhomogeneous models.}

  \item{piecewise.covariates}{Currently ignored: not implemented for
    time-inhomogeneous models.}

  \item{ci}{If \code{"normal"}, then calculate a confidence interval
    by simulating \code{B} random vectors
    from the asymptotic multivariate normal distribution implied by the
    maximum likelihood estimates (and covariance matrix) of the log
    transition intensities and covariate effects.

    If \code{"bootstrap"} then calculate a confidence interval by
    non-parametric bootstrap refitting.  This is 1-2 orders of magnitude
    slower than the \code{"normal"} method, but is expected to be more
    accurate. See \code{\link{boot.msm}} for more details of
    bootstrapping in \pkg{msm}.

    If \code{"none"} (the default) then no confidence interval is
    calculated.}

  \item{cl}{Width of the symmetric confidence interval, relative to 1.}

  \item{B}{Number of bootstrap replicates.}

  \item{cores}{Number of cores to use for bootstrapping using parallel
    processing. See \code{\link{boot.msm}} for more details.}

  \item{...}{Arguments to pass to \code{\link{MatrixExp}}.}
}
\value{
  A matrix whose \eqn{r, s} entry is the probability of having visited
  state \eqn{s} at least once before time \eqn{t}, given the state at
  time \eqn{0} is \eqn{r}.    The diagonal entries should all be 1.
}
\details{
  The passage probabilities to state \eqn{s} are computed by
  setting the \eqn{s}th row of the transition intensity matrix \eqn{Q}
  to zero, giving an intensity matrix \eqn{Q*} for a simplified model structure
  where state \eqn{s} is absorbing.  The probabilities of passage are
  then equivalent to row \eqn{s} of the transition probability matrix
  \eqn{Exp(tQ*)} under this simplified model for \eqn{t=}\code{tot}.

  Note this is different from the probability of occupying each state at
  exactly time \eqn{t}, given by \code{\link{pmatrix.msm}}.  The passage
  probability allows for the possibility of having visited the state
  before \eqn{t}, but then occupying a different state at \eqn{t}.

  The mean of the passage distribution is the expected first passage
  time, \code{\link{efpt.msm}}.

  This function currently only handles time-homogeneous Markov models.
  For time-inhomogeneous models the covariates are held constant at the
  value supplied, by default the column means of the design matrix over
  all observations.
}
\references{
Norris, J. R. (1997) Markov Chains. Cambridge University Press.
}
\seealso{
  \code{\link{efpt.msm}}, \code{\link{totlos.msm}}, \code{\link{boot.msm}}.
}
\examples{
Q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
           c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))

## ppass[1,2](t) converges to 0.5 with t, since given in state 1, the
## probability of going to the absorbing state 4 before visiting state
## 2 is 0.5, and the chance of still being in state 1 at t decreases.

ppass.msm(qmatrix=Q, tot=2)
ppass.msm(qmatrix=Q, tot=20)
ppass.msm(qmatrix=Q, tot=100)

Q <- Q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)

## Probability of about 1/2 of visiting state 3 by time 10.5, the
## median first passage time

ppass.msm(qmatrix=Q, tot=10.5)

## Mean first passage time from state 2 to state 3 is 10.02: similar
## to the median

efpt.msm(qmatrix=Q, tostate=3)
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
