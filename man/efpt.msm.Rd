\name{efpt.msm}
\alias{efpt.msm}
\title{Expected first passage time}
\description{
  Expected time until first reaching a particular state or set of states
  in a Markov model.
}
\usage{
efpt.msm(x=NULL, qmatrix=NULL, tostate, start="all", covariates="mean",
           ci=c("none","normal","bootstrap"), cl=0.95, B=1000,
           cores=NULL, \dots)
}
\arguments{
  \item{x}{A fitted multi-state model, as returned by
    \code{\link{msm}}. }

  \item{qmatrix}{Instead of \code{x}, you can simply supply a transition
  intensity matrix in \code{qmatrix}.}

  \item{tostate}{State, or set of states supplied as a vector, for which
    to estimate the first passage time into.  Can be integer, or character
    matched to the row names of the Q matrix.}

  \item{start}{Starting state (integer).  By default (\code{start="all"}),
    this will return a vector of expected passage times from each
    state in turn.

    Alternatively, this can be used to obtain the expected first passage
    time from a \emph{set} of states, rather than single states.  To
    achieve this, \code{state} is set to a vector of weights, with
    length equal to the number of states in the model.  These weights
    should be proportional to the probability of starting in each of the
    states in the desired set, so that weights of zero are supplied for
    other states.  The function will calculate the weighted average of
    the expected passage times from each of the corresponding states.
  }

  \item{covariates}{Covariate values defining the intensity matrix for
  the fitted model \code{x}, as supplied to \code{\link{qmatrix.msm}}.}

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
  A vector of expected first passage times, or "hitting times", from each state to the
  desired state.
}
\details{
  The expected first passage times from each of a set of states
  \eqn{\mathbf{i}}{i}
  to to the remaining set of states \eqn{\overline{\mathbf{i}}}{ibar} in the
  state space, for a model with transition intensity matrix \eqn{Q}, are

  \deqn{-Q_{\mathbf{i},\mathbf{i}}^{-1} \mathbf{1}}{-Q_{i,i}^{-1} 1}

  where \eqn{\mathbf{1}}{1} is a vector of ones, and \eqn{Q_{\mathbf{i},\mathbf{i}}}{Q_{i,i}}
  is the square subset of \eqn{Q} pertaining to states \eqn{\mathbf{i}}{i}.

  It is equal to the sum of mean sojourn times for all states
  between the "from" and "to" states in a unidirectional model.
  If there is non-zero chance of reaching an absorbing state
  before reaching \code{tostate}, then it is infinite.  It is
  trivially zero if the "from" state equals \code{tostate}.

  This function currently only handles time-homogeneous Markov models.
  For time-inhomogeneous models it will assume that \eqn{Q} equals the
  average intensity matrix over all times and observed covariates.
  Simulation might be used to handle time dependence.

  Note this is the \emph{expectation} of first passage time, and the
  confidence intervals are CIs for this mean, not predictive intervals
  for the first passage time.  The full distribution of the first
  passage time to a set of states can be obtained by setting the rows of
  the intensity matrix \eqn{Q} corresponding to that set of states to
  zero to make a model where those states are absorbing.  The
  corresponding transition probability matrix \eqn{Exp(Qt)} then gives
  the probabilities of having hit or passed that state by a time \eqn{t}
  (see the example below). This is implemented in
  \code{\link{ppass.msm}}.
  
}
\references{
Norris, J. R. (1997) Markov Chains. Cambridge University Press.
}
\seealso{
  \code{\link{sojourn.msm}}, \code{\link{totlos.msm}}, \code{\link{boot.msm}}.
}
\examples{
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
             c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
efpt.msm(qmatrix=twoway4.q, tostate=3)
# given in state 1, expected time to reaching state 3 is infinite
# since may die (state 4) before entering state 3

# If we remove the death state from the model, EFPTs become finite
Q <- twoway4.q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
efpt.msm(qmatrix=Q, tostate=3)

# Suppose we cannot die or regress while in state 2, can only go to state 3
Q <- twoway4.q; Q[2,4] <- Q[2,1] <- 0; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
efpt.msm(qmatrix=Q, tostate=3)
# The expected time from 2 to 3 now equals the mean sojourn time in 2.
-1/Q[2,2]

# Calculate cumulative distribution of the first passage time
# into state 3 for the following three-state model
Q <- twoway4.q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
# Firstly form a model where the desired hitting state is absorbing
Q[3,] <- 0
MatrixExp(Q, t=10)[,3]
ppass.msm(qmatrix=Q, tot=10)
# Given in state 1 at time 0, P(hit 3 by time 10) = 0.479
MatrixExp(Q, t=50)[,3]  # P(hit 3 by time 50) = 0.98
ppass.msm(qmatrix=Q, tot=50)

}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
