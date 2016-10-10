\name{boot.msm}
\alias{boot.msm}
\title{Bootstrap resampling for multi-state models}
\description{
  Draw a number of bootstrap resamples, refit a \code{\link{msm}} model to the
  resamples, and calculate statistics on the refitted models.
}
\usage{
boot.msm(x, stat=pmatrix.msm, B=1000, file=NULL, cores=NULL)
}
\arguments{
  \item{x}{A fitted msm model, as output by \code{\link{msm}}.}
  \item{stat}{A function to call on each refitted msm model. By default
    this is \code{\link{pmatrix.msm}}, returning the transition
    probability matrix in one time unit. If
    \code{NULL} then no function is computed.}

  \item{B}{Number of bootstrap resamples.}

  \item{file}{Name of a file in which to save partial results after each
    replicate. This is saved using \code{\link{save}} and can be
    restored using \code{\link{load}}, producing an object called
    \code{boot.list} containing the partial results.  Not 
    supported when using parallel processing.}

  \item{cores}{Number of processor cores to use for parallel processing.
  Requires the \pkg{doParallel} package to be installed.  If not
  specified, parallel processing is not used. If \code{cores} is set to
  the string \code{"default"}, the default methods of
  \code{\link[parallel]{makeCluster}} (on Windows) or
  \code{\link[doParallel]{registerDoParallel}} (on Unix-like) are
  used.
}

}
\details{
  The bootstrap datasets are computed by resampling independent
  transitions between pairs of states (for non-hidden models without
  censoring), or independent individual series (for hidden models or
  models with censoring).  Therefore this approach doesn't work if,
  for example, the data for a HMM consist of a series of observations
  from just one individual, and is inaccurate for small numbers of
  independent transitions or individuals.

  Confidence intervals or standard errors for the corresponding
  statistic can be calculated by summarising the returned list of
  \code{B} replicated outputs.  This is currently implemented for most
  the output functions \code{\link{qmatrix.msm}}, \code{\link{ematrix.msm}},
  \code{\link{qratio.msm}},
  \code{\link{pmatrix.msm}},   \code{\link{pmatrix.piecewise.msm}},
  \code{\link{totlos.msm}} and \code{\link{prevalence.msm}}.
  For other outputs, users will have to write their own code to
  summarise the output of \code{\link{boot.msm}}.

  Most of \pkg{msm}'s output functions present confidence intervals
  based on asymptotic standard errors calculated from the
  Hessian. These are expected to be underestimates of the true standard
  errors (Cramer-Rao lower bound).  Some of these functions use a
  further approximation, the delta method (see
  \code{\link{deltamethod}}) to obtain standard errors of transformed
  parameters. Bootstrapping should give a more
  accurate estimate of the uncertainty.

  An alternative method which is less accurate though faster than
  bootstrapping, but more accurate than the delta method, is to draw a
  sample from the asymptotic multivariate normal distribution implied by
  the maximum likelihood estimates (and covariance matrix), and
  summarise the transformed estimates.  See \code{\link{pmatrix.msm}}.

  All objects used in the original call to \code{\link{msm}} which
  produced \code{x}, such as the \code{qmatrix}, should be in the
  working environment, or else \code{boot.msm} will produce an
  \dQuote{object not found} error.  This enables \code{boot.msm} to
  refit the original model to the replicate datasets.  However there is
  currently a limitation.  In the original call to \code{msm}, the
  \code{"formula"} argument should be specified directly, as, for
  example,

  \code{msm(state ~ time, data = ...)}

  and not, for example,

  \code{form = data$state ~ data$time}

  \code{msm(formula=form, data = ...)}

  otherwise \code{boot.msm} will be unable to draw the replicate
  datasets.

  \code{boot.msm} will also fail with an incomprehensible error if the
  original call to msm used a used-defined object whose name is the same
  as a built-in R object, or an object in any other loaded package.  For
  example, if you have called a Q matrix \code{q}, when \code{q()} is
  the built-in function for quitting R.

  If \code{stat} is \code{NULL}, then \code{B} different \code{msm} model
  objects will be stored in memory. This is unadvisable, as \code{msm}
  objects tend to be large, since they contain the original data used for
  the \code{msm} fit, so this will be wasteful of memory.

  To specify more than one statistic, write a function consisting of a
  list of different function calls, for example,

  \code{stat = function(x) list (pmatrix.msm(x, t=1), pmatrix.msm(x,
    t=2))}
}
\value{
  A list with \code{B} components, containing the result of calling
  function \code{stat} on each of the refitted models.  If \code{stat}
  is \code{NULL}, then each component just contains the refitted
  model.  If one of the \code{B} model fits was unsuccessful and
  resulted in an error, then the corresponding list component will
  contain the error message.
}
\references{
  Efron, B. and Tibshirani, R.J. (1993) \emph{An Introduction to the
    Bootstrap}, Chapman and Hall.
}
\seealso{
  \code{\link{qmatrix.msm}}, \code{\link{qratio.msm}},
  \code{\link{sojourn.msm}}, \code{\link{ematrix.msm}}, \code{\link{pmatrix.msm}}, \code{\link{pmatrix.piecewise.msm}}, \code{\link{totlos.msm}}, \code{\link{prevalence.msm}}.
}
\examples{
\dontrun{
  ## Psoriatic arthritis example
  data(psor)
  psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix =
    psor.q, covariates = ~ollwsdrt+hieffusn,
    constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)),
    control = list(REPORT=1,trace=2), method="BFGS")
  ## Bootstrap the baseline transition intensity matrix.  This will take a long time.
  q.list <- boot.msm(psor.msm, function(x)x$Qmatrices$baseline)
  ## Manipulate the resulting list of matrices to calculate bootstrap standard errors.
  apply(array(unlist(q.list), dim=c(4,4,5)), c(1,2), sd)
  ## Similarly calculate a bootstrap 95\% confidence interval
  apply(array(unlist(q.list), dim=c(4,4,5)), c(1,2),
        function(x)quantile(x, c(0.025, 0.975)))
  ## Bootstrap standard errors are larger than the asymptotic standard
  ## errors calculated from the Hessian
  psor.msm$QmatricesSE$baseline
}
}
\author{C.H.Jackson <chris.jackson@mrc-bsu.cam.ac.uk>}
\keyword{models}
