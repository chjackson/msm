\name{paramdata.object}
\alias{paramdata.object}
\title{Developer documentation: internal msm parameters object}
\description{
  An object giving information about the parameters of the multi-state
  model.  Used internally during maximum likelihood estimation and
  arranging results.  Returned as the \code{paramdata} component of a
  fitted \code{\link{msm}} model object.
}
\value{
  \item{inits}{Vector of initial values for distinct parameters which
    are being estimated.  These have been transformed to the real line
    (e.g. by log),
    and exclude parameters being fixed at their initial values,
    parameters defined to be always fixed (e.g. binomial denominators) and
    parameters constrained to equal previous ones.}
\item{plabs}{Names of parameters in \code{allinits}.}
\item{allinits}{Vector of parameter values before estimation, including
  those which are fixed or constrained to equal other parameters, and
  transformed to the real line.}
\item{hmmpars}{Indices of \code{allinits} which represent baseline parameters of
  hidden Markov outcome models (thus excluding covariate effects in HMMs
  and initial state occupancy probabilities). }
\item{fixed}{\code{TRUE} if all parameters are fixed, \code{FALSE} otherwise.}
\item{fixedpars}{Indices of parameters in \code{allinits} which are
  fixed, either by definition or as requested by the user in the
  \code{fixedpars} argument to \code{\link{msm}}. Excludes parameters
  fixed by constraining to equal other parameters.}
\item{notfixed}{Indices of parameters which are not fixed by the
  definition of \code{fixedpars}.}
\item{optpars}{Indices of parameters in \code{allinits} being estimated,
  thus those included in \code{inits}.}
\item{auxpars}{Indices of "auxiliary" parameters which are always fixed,
  for example, binomial denominators (\code{\link{hmmBinom}}) and the \code{which} parameter in
  \code{\link{hmmIdent}}.}
\item{constr}{Vector of integers, of length \code{npars}, indicating
  which sets of parameters are constrained to be equal to each other.
  If two of these integers are equal the corresponding parameters are
  equal.  A negative element indicates that parameter is defined to be
  minus some other parameter (this is used for covariate effects on
  transition intensities).}
\item{npars}{Total number of parameters, equal to  \code{length(allinits)}.}
\item{nfix}{Number of fixed parameters, equal to \code{length(fixedpars)}. }
\item{nopt}{Number of parameters being estimated, equal to
  \code{length(inits)} and \code{length(optpars)}.}
\item{ndup}{Number of parameters defined as duplicates of previous
  parameters by equality constraints (currently unused).}
\item{ranges}{Matrix of defined ranges for each parameter on the natural scale
  (e.g. 0 to infinity for rate parameters). }
\item{opt}{Object returned by the optimisation routine (such as \code{\link{optim}}).}
\item{foundse}{\code{TRUE} if standard errors are available after
  optimisation. If \code{FALSE} the optimisation probably hasn't converged. }
\item{lik}{Minus twice the log likelihood at the parameter estimates.}
\item{deriv}{Derivatives of the minus twice log likelihood at the
  parameter estimates, if available. }
\item{information}{Corresponding expected information matrix at the
  parameter estimates, if available.}
\item{params}{Vector of parameter values after maximum likelihood
  estimation, corresponding to \code{allinits}, still on the real-line
  transformed scale.}
\item{covmat}{Covariance matrix corresponding to \code{params}.}
\item{ci}{Matrix of confidence intervals corresponding to \code{params},
with nominal coverage (default 0.95) defined by the \code{cl} argument
of \code{\link{msm}}. }
\item{estimates.t}{Vector of parameter estimates, as \code{params} but
  with parameters on their natural scales.}
}
\seealso{
\code{\link{msm.object}}
}
