\name{viterbi.msm}
\title{Calculate the probabilities of underlying states and the most
  likely path through them}
\alias{viterbi.msm}

\description{
  For a fitted hidden Markov model, or a model with censored state
  observations, the Viterbi algorithm recursively constructs the path
  with the highest probability through the underlying states.  The
  probability of each hidden state is also computed for hidden Markov
  models.
}

\usage{
viterbi.msm(x)
}
\arguments{
  \item{x}{A fitted hidden Markov multi-state model, or a model with
    censored state observations, as produced by \code{\link{msm}}}
}
\value{
    A data frame with columns:

    \code{subject} = subject identification numbers

    \code{time} = times of observations

    \code{observed} = corresponding observed states

    \code{fitted} = corresponding fitted states found by Viterbi
    recursion. If the model is not a hidden Markov model and there are
    no censored state observations, this is just the observed states. 

    For hidden Markov models, an additional matrix \code{pstate} is also
    returned inside the data frame, giving the probability of each
    hidden state at each point, conditionally on all the data.  This is
    computed by the forward/backward algorithm.
  }
\seealso{
  \code{\link{msm}}
}
\references{
  Durbin, R., Eddy, S., Krogh, A. and Mitchison, G. \emph{Biological
    sequence analysis}, Cambridge University Press, 1998.
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
