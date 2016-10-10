\name{draic.msm}
\alias{draic.msm}
\alias{drlcv.msm}
\title{Criteria for comparing two multi-state models with nested state spaces}
\description{
  A modification of Akaike's information criterion, and a leave-one-out likelihood
  cross-validation criterion, for comparing the predictive ability of two
  Markov multi-state models with nested state spaces.  This is evaluated based
  on the restricted or aggregated data which the models have in common.

  Note that standard AIC can be computed for one or more fitted
  \code{msm} models \code{x,y,...} using \code{\link{AIC}(x,y,...)}, and
  this can be used to compare models fitted to the same
  data. \code{draic.msm} and \code{drlcv.msm} are designed for models
  fitted to data with differently-aggregated state spaces.
}
\usage{
draic.msm(msm.full, msm.coarse, likelihood.only=FALSE,
         information=c("expected","observed"), tl=0.95)
drlcv.msm(msm.full, msm.coarse, tl=0.95, cores=NULL,
          verbose=TRUE,outfile=NULL)
}
\arguments{
  \item{msm.full}{Model on the bigger state space.}

  \item{msm.coarse}{Model on the smaller state space.

    The two models must both be non-hidden Markov models without
    censored states.
    
    The two models must be fitted to the same datasets, except that the
    state space of the coarse model must be an aggregated version of the
    state space of the full model.  That is, every state in the full
    dataset must correspond to a unique state in the coarse dataset.
    For example, for the full state variable \code{c(1,1,2,2,3,4)}, the
    corresponding coarse states could be \code{c(1,1,2,2,2,3)}, but not
    \code{c(1,2,3,4,4,4)}.

    The structure of allowed transitions in the coarse model must also
    be a collapsed version of the big model structure, but no check is
    currently made for this in the code.
    
    To use these functions, all objects
    which were used in the calls to fit \code{msm.full} and
    \code{msm.coarse} must be in the working environment, for example,
    datasets and definitions of transition matrices.
  }

  \item{likelihood.only}{Don't calculate Hessians and trace term (DRAIC).}

  \item{information}{Use observed or expected information in the DRAIC
  trace term.  Expected is the default, and much faster, though is only
  available for models fitted to pure panel data (all \code{obstype=1}
  in the call to \code{\link{msm}}, thus not exact transition times or
  exact death times) }

  \item{tl}{Width of symmetric tracking interval, by default 0.95 for a
    95\% interval.}

  \item{cores}{Number of processor cores to use in \code{drlcv} for
    cross-validation by parallel processing.  Requires the
    \pkg{doParallel} package to be installed.  If not specified,
    parallel processing is not used. If \code{cores} is set to the
    string \code{"default"}, the default methods of
    \code{\link[parallel]{makeCluster}} (on Windows) or
    \code{\link[doParallel]{registerDoParallel}} (on Unix-like) are
    used.}

  \item{verbose}{Print intermediate results of each iteration of
    cross-validation to the console while running. May not work with
    parallel processing. }

  \item{outfile}{Output file to print intermediate results of
    cross-validation.  Useful to track execution speed when using
    parallel processing, where output to the console may not work.}
}
\value{
  A list containing \eqn{D_{RAIC}}{D_RAIC} (\code{draic.msm}) or
    \eqn{D_{RLCV}}{D_RLCV} (\code{drlcv.msm}), its component terms, and tracking intervals.
}
\details{
  
  The difference in restricted AIC (Liquet and Commenges, 2011), as
  computed by this function, is defined as

  \deqn{D_{RAIC} = l(\gamma_n |\mathbf{x}'' ) - l(\theta_n |\mathbf{x}'' ) + trace ( J(\theta_n |\mathbf{x}'')J(\theta_n |\mathbf{x})^{-1} - J(\gamma_n |\mathbf{x}'' )J(\gamma_n |\mathbf{x}' )^{-1})}{D_RAIC = l(gamma_n |x'' ) - l(theta_n |x'' ) + trace ( J(theta_n |x'')J(theta_n |x)^{-1} - J(gamma_n |x'' )J(gamma_n |x' )^{-1})}

  where \eqn{\gamma}{gamma} and \eqn{\theta}{theta} are the maximum likelihood estimates of the smaller and bigger models, fitted to the smaller and bigger data, respectively.

  \eqn{l(\gamma_n |x'')}{l(gamma_n |x'')} represents the likelihood of the simpler model evaluated on the restricted data.

  \eqn{l(\theta_n |x'')}{l(theta_n |x'')} represents the likelihood of
  the complex model evaluated on the restricted data.  This is a hidden
  Markov model, with a misclassification matrix and initial state
  occupancy probabilities as described by Thom et al (2014).

  \eqn{J()} are the corresponding (expected or observed, as specified by the user) information matrices.

  \eqn{\mathbf{x}}{x} is the expanded data, to which the bigger model was originally fitted, and
  \eqn{\mathbf{x}'}{x'} is the data to which the smaller model was originally fitted.
  \eqn{\mathbf{x}''}{x''} is the restricted data which the two models
  have in common.  \eqn{\mathbf{x}'' = \mathbf{x}'}{x'' = x} in this implementation, so the models are nested.

  The difference in likelihood cross-validatory criteria  (Liquet and Commenges, 2011) is defined as

  \deqn{D_{RLCV} = 1/n \sum_{i=1}^n \log( h_{X''}(x_i'' | \gamma_{-i}) / g_{X''}(x_i''| \theta_{-i}))}{D_{RLCV} = 1/n \sum_{i=1}^n log( h_{X''}(x_i'' | gamma_{-i}) / g_{X''}(x_i''| theta_{-i}))}

  where \eqn{\gamma_{-i}} and \eqn{\theta_{-i}} are the maximum likelihood estimates from the smaller and bigger models fitted to datasets with subject \eqn{i} left out, \eqn{g()} and \eqn{h()} are the densities of the corresponding models, and \eqn{x_i''} is the restricted data from subject \eqn{i}.

  Tracking intervals are analogous to confidence intervals, but not strictly the same, since the quantity which D_RAIC aims to estimate, the difference in expected Kullback-Leibler discrepancy for predicting a replicate dataset, depends on the sample size.  See the references.

  Positive values for these criteria indicate the coarse model is
  preferred, while negative values indicate the full model is preferred.
}
\references{
  Thom, H. and Jackson, C. and Commenges, D. and Sharples, L. (2015)
  State selection in multistate models with application to quality of
  life in psoriatic arthritis.  Statistics In Medicine 34(16) 2381 -
  2480. 

  Liquet, B. and Commenges D. (2011) Choice of estimators based on
  different observations: Modified AIC and LCV criteria. Scandinavian
  Journal of Statistics; 38:268-287.
}
\seealso{
  \code{\link{logLik.msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk},
  H. H. Z. Thom \email{howard.thom@bristol.ac.uk}}
\keyword{models}
