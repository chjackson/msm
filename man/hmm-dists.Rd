\name{hmm-dists}
\alias{hmm-dists}
\alias{hmmCat}
\alias{hmmIdent}
\alias{hmmUnif}
\alias{hmmNorm}
\alias{hmmLNorm}
\alias{hmmExp}
\alias{hmmGamma}
\alias{hmmWeibull}
\alias{hmmPois}
\alias{hmmBinom}
\alias{hmmTNorm}
\alias{hmmMETNorm}
\alias{hmmMEUnif}
\alias{hmmNBinom}
\alias{hmmBeta}
\alias{hmmT}
\title{Hidden Markov model constructors}

\description{  
These functions are used to specify the distribution of the response
conditionally on the underlying state in a hidden Markov model.  A list
of these function calls, with one component for each state, should be 
used for the \code{hmodel} argument to \code{msm}. The initial values
for the parameters of the distribution should be given as arguments.
}
\usage{
hmmCat(prob, basecat)
hmmIdent(x)
hmmUnif(lower, upper)
hmmNorm(mean, sd)
hmmLNorm(meanlog, sdlog)
hmmExp(rate)
hmmGamma(shape, rate)
hmmWeibull(shape, scale)
hmmPois(rate)
hmmBinom(size, prob)
hmmTNorm(mean, sd, lower, upper)
hmmMETNorm(mean, sd, lower, upper, sderr, meanerr=0)
hmmMEUnif(lower, upper, sderr, meanerr=0)
hmmNBinom(disp, prob)
hmmBeta(shape1,shape2)
hmmT(mean,scale,df)
}
\arguments{
  \item{prob}{ (\code{hmmCat}) Vector of probabilities of observing
    category \code{1, 2, \dots, length(prob)} respectively.  Or
    the probability governing a binomial or negative binomial
    distribution.}

  \item{basecat}{ (\code{hmmCat}) Category which is considered to be the "baseline",
    so that during estimation, the probabilities are parameterised as
    probabilities relative to this baseline category. By default, the
    category with the greatest probability is
    used as the baseline. }
  
  \item{x}{ (\code{hmmIdent}) Code in the data which denotes the exactly-observed state.  } 
  
  \item{mean}{ (\code{hmmNorm,hmmLNorm,hmmTNorm}) Mean defining a Normal, or truncated Normal
    distribution. } 

  \item{sd}{ (\code{hmmNorm,hmmLNorm,hmmTNorm}) Standard deviation defining a
    Normal, or truncated Normal distribution. } 
  
  \item{meanlog}{ (\code{hmmNorm,hmmLNorm,hmmTNorm}) Mean on the log
    scale, for a log Normal distribution. } 

  \item{sdlog}{ (\code{hmmNorm,hmmLNorm,hmmTNorm}) Standard deviation on
    the log scale, for a log Normal distribution. } 

  \item{rate}{ (\code{hmmPois,hmmExp,hmmGamma}) Rate of a Poisson,
    Exponential or Gamma distribution (see \code{\link{dpois}}, \code{\link{dexp}}, \code{\link{dgamma}}). } 

  \item{shape}{ (\code{hmmPois,hmmExp,hmmGamma}) Shape parameter of a
    Gamma or Weibull distribution (see \code{\link{dgamma}},
    \code{\link{dweibull}}). } 

  \item{shape1,shape2}{First and second parameters of a beta
    distribution (see \code{\link{dbeta}}).}
  
  \item{scale}{ (\code{hmmGamma}) Scale parameter of a Gamma
    distribution (see \code{\link{dgamma}}), or unstandardised Student t
    distribution. }

  \item{df}{Degrees of freedom of the Student t distribution.}
  
  \item{size}{Order of a Binomial distribution (see \code{\link{dbinom}}).}
  
  \item{disp}{Dispersion parameter of a negative binomial distribution,
    also called \code{size} or \code{order}.  (see
    \code{\link{dnbinom}}).}
  
  \item{lower}{ (\code{hmmUnif,hmmTNorm,hmmMEUnif}) Lower limit for an Uniform or truncated Normal distribution. } 

  \item{upper}{ (\code{hmmUnif,hmmTNorm,hmmMEUnif}) Upper limit for an Uniform or truncated Normal
    distribution. } 

  \item{sderr}{ (\code{hmmMETNorm,hmmUnif}) Standard deviation of the Normal measurement error
    distribution. } 

  \item{meanerr}{ (\code{hmmMETNorm,hmmUnif}) Additional shift in the
    measurement error, fixed to 0 by default.  This may
    be modelled in terms of covariates. } 
}
\value{
   Each function returns an object of class \code{hmodel}, which is a
   list containing information about the model.  The only component
   which may be useful to end users is \code{r}, a function of one
   argument \code{n} which returns a random sample of size \code{n} from
   the given distribution. 
}
\details{
  \code{hmmCat} represents a categorical response distribution on the
  set \code{1, 2, \dots, length(prob)}.  The
  Markov model with misclassification is an example of this type of model. The
  categories in this case are (some subset of) the underlying states. 
  
  The \code{hmmIdent} distribution is used for underlying states which are
  observed exactly without error.

  \code{hmmUnif}, \code{hmmNorm}, \code{hmmLNorm}, \code{hmmExp},
  \code{hmmGamma}, \code{hmmWeibull}, \code{hmmPois}, \code{hmmBinom}, 
   \code{hmmTNorm}, \code{hmmNBinom} and \code{hmmBeta} represent Uniform, Normal,
   log-Normal, exponential, Gamma, Weibull, Poisson, Binomial, truncated Normal,
  negative binomial and beta distributions, respectively, with parameterisations
  the same as the default parameterisations in the corresponding base R
  distribution functions. 

  \code{hmmT} is the Student t distribution with general mean
  \eqn{\mu}{mu}, scale \eqn{\sigma}{sigma} and degrees of freedom
  \code{df}.  
  The variance is \eqn{\sigma^2 df/(df + 2)}{sigma^2 df/(df + 2)}.
  Note the t distribution in base R \code{\link{dt}} is a standardised one with
  mean 0 and scale 1.  These allow any positive (integer or non-integer)
  \code{df}.    By default, all three
  parameters, including \code{df}, are estimated when fitting a hidden
  Markov model, but in practice, \code{df} might need to be fixed for identifiability - this can be done
  using the \code{fixedpars} argument to \code{\link{msm}}. 
  
  The \code{hmmMETNorm} and \code{hmmMEUnif} distributions are
  truncated Normal and Uniform distributions, but with additional Normal measurement error on the
  response. These are generalisations of the distributions proposed by
Satten and Longini (1996) for modelling the progression of CD4 cell
counts in monitoring HIV disease.  See \code{\link{medists}} for
density, distribution, quantile and random generation functions for
these distributions.   See also \code{\link{tnorm}} for
density, distribution, quantile and random generation functions for
the truncated Normal distribution.

See the PDF manual  \file{msm-manual.pdf} in the \file{doc}
  subdirectory for algebraic definitions of all these distributions.
  New hidden Markov model response distributions can be added to
  \pkg{msm} by following the instructions in Section 2.17.1.
  
  Parameters which can be modelled in terms of covariates, on the scale
  of a link function, are as follows.

  \tabular{ll}{
    PARAMETER NAME \tab LINK FUNCTION \cr
    \code{mean} \tab identity \cr
    \code{meanlog} \tab identity \cr
    \code{rate} \tab log \cr
    \code{scale} \tab log \cr
    \code{meanerr} \tab identity \cr
    \code{prob} \tab (multinomial logistic regression)
  }

  Parameters \code{basecat, lower, upper, size, meanerr} are fixed at
  their initial values. All other parameters are estimated while fitting
  the hidden Markov model, unless the appropriate \code{fixedpars}
  argument is supplied to \code{msm}. 

  For categorical response distributions \code{(hmmCat)} the
  outcome probabilities initialized to zero are fixed at zero, and the
  probability corresponding to \code{basecat} is fixed to one minus the
  sum of the remaining probabilities.  These remaining probabilities are
  estimated, and can be modelled in terms of covariates via multinomial
  logistic regression (relative to \code{basecat}).
  
}
\references{
  Satten, G.A. and Longini, I.M.  Markov chains with measurement error:
  estimating the 'true' course of a marker of the progression of human
  immunodeficiency virus disease (with discussion) \emph{Applied
    Statistics} 45(3): 275-309 (1996).
    
  Jackson, C.H. and Sharples, L.D. Hidden Markov models for the
  onset and progresison of bronchiolitis obliterans syndrome in lung
  transplant recipients \emph{Statistics in Medicine}, 21(1): 113--128
  (2002).

  Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W. and
  Couto, E.  Multi-state Markov models for disease progression with
  classification error. \emph{The Statistician}, 52(2): 193--209 (2003).
}
\seealso{
  \code{\link{msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{distribution}
