#' Multi-State Markov and Hidden Markov Models in Continuous Time
#'
#' msm: Functions for fitting continuous-time Markov and hidden Markov
#' multi-state models to longitudinal data.  Designed for processes
#' observed at arbitrary times in continuous time (intermittently
#' observed or panel data) but some other observation schemes are
#' supported. Both Markov transition rates and the hidden Markov
#' output process can be modelled in terms of covariates, which may be
#' constant or piecewise-constant in time.
#' 
#' @name msm-package
#' @aliases msm-package
#' @docType package
#' @useDynLib msm, .registration=TRUE
#'
#' @importFrom graphics plot persp contour image filled.contour legend lines par text
#' @importFrom grDevices rainbow 
#' @importFrom stats coef as.formula deriv dexp dnorm integrate logLik model.extract model.frame model.matrix na.fail na.omit na.pass numericDeriv optimHess pchisq pexp plogis pnorm qlogis qnorm quantile rbeta rbinom reformulate rexp rgamma rlnorm rnbinom rnorm rpois rt runif rweibull sd setNames terms uniroot
#' @importFrom utils head tail
#' @importFrom mvtnorm rmvnorm
#' @importFrom survival Surv survfit
#' @importFrom expm expm
#' 
"_PACKAGE"


#' Coxian phase-type distribution with two phases
#' 
#' Density, distribution, quantile functions and other utilities for the Coxian
#' phase-type distribution with two phases.
#' 
#' This is the distribution of the time to reach state 3 in a continuous-time
#' Markov model with three states and transitions permitted from state 1 to
#' state 2 (with intensity \eqn{\lambda_1}{lambda1}) state 1 to state 3
#' (intensity \eqn{\mu_1}{mu1}) and state 2 to state 3 (intensity
#' \eqn{\mu_2}{mu2}).  States 1 and 2 are the two "phases" and state 3 is the
#' "exit" state.
#' 
#' The density is
#' 
#' \deqn{f(t | \lambda_1, \mu_1) = e^{-(\lambda_1+\mu_1)t}(\mu_1 +
#' (\lambda_1+\mu_1)\lambda_1 t)}{f(t | l1, mu1) = exp(-(l1+mu1)*t)*(mu1 +
#' (l1+mu1)*l1*t)}
#' 
#' if \eqn{\lambda_1 + \mu_1 = \mu_2}{l1 + mu1 = mu2}, and
#' 
#' \deqn{f(t | \lambda_1, \mu_1, \mu_2) =
#' \frac{(\lambda_1+\mu_1)e^{-(\lambda_1+\mu_1)t}(\mu_2-\mu_1) +
#' \mu_2\lambda_1e^{-\mu_2t}}{\lambda_1+\mu_1-\mu_2}}{f(t | l1, mu1, mu2) =
#' ((l1+mu1)*exp(-(l1+mu1)*t)*(mu2-mu1) + mu2*l1*exp(-mu2*t))/(l1+mu1-mu2)}
#' 
#' otherwise.  The distribution function is
#' 
#' \deqn{F(t | \lambda_1, \mu_1) = 1 - e^{-(\lambda_1+\mu_1) t} (1 + \lambda_1
#' t)}{F(t | l1, mu1) = 1 - exp(-(l1+mu1)*t)*(1 + l1*t)}
#' 
#' if \eqn{\lambda_1 + \mu_1 = \mu_2}{l1 + mu1 = mu2}, and
#' 
#' \deqn{F(t | \lambda_1, \mu_1, \mu_2) = 1 - \frac{e^{-(\lambda_1+\mu_1) }{F(t
#' | l1, mu1, mu2) = 1 - (exp(-(l1+mu1)*t)*(-mu1+mu2) +
#' l1*exp(-mu2*t))/(l1+mu1-mu2)}\deqn{ t} (\mu_2 - \mu_1) + \lambda_1 e^{-\mu_2
#' }{F(t | l1, mu1, mu2) = 1 - (exp(-(l1+mu1)*t)*(-mu1+mu2) +
#' l1*exp(-mu2*t))/(l1+mu1-mu2)}\deqn{ t}}{\lambda_1+\mu_1-\mu_2}}{F(t | l1,
#' mu1, mu2) = 1 - (exp(-(l1+mu1)*t)*(-mu1+mu2) + l1*exp(-mu2*t))/(l1+mu1-mu2)}
#' otherwise.  Quantiles are calculated by numerically inverting the
#' distribution function.
#' 
#' The mean is \eqn{(1 + \lambda_1/\mu_2) / (\lambda_1 + \mu_1)}{(1 + l1/mu2) /
#' (l1 + mu1)}.
#' 
#' The variance is \eqn{(2 + 2\lambda_1(\lambda_1+\mu_1+ \mu_2)/\mu_2^2 - (1 +
#' \lambda_1/\mu_2)^2)/(\lambda_1+\mu_1)^2}{(2 + 2*l1*(l1+mu1+ mu2)/mu2^2 - (1
#' + l1/mu2)^2)/(l1+mu1)^2}.
#' 
#' If \eqn{\mu_1=\mu_2}{mu1=mu2} it reduces to an exponential distribution with
#' rate \eqn{\mu_1}{mu1}, and the parameter \eqn{\lambda_1}{l1} is redundant.
#' Or also if \eqn{\lambda_1=0}{l1=0}.
#' 
#' The hazard at \eqn{x=0} is \eqn{\mu_1}, and smoothly increasing if
#' \eqn{\mu_1<\mu_2}{mu1<mu2}.  If \eqn{\lambda_1 + \mu_1 \geq \mu_2}{l1 + mu1
#' >= mu2} it increases to an asymptote of \eqn{\mu_2}{mu2}, and if
#' \eqn{\lambda_1 + \mu_1 \leq \mu_2}{l1 + mu1 <= mu2} it increases to an
#' asymptote of \eqn{\lambda_1 + \mu_1}{l1 + mu1}.  The hazard is decreasing if
#' \eqn{\mu_1>\mu_2}{mu1>mu2}, to an asymptote of \eqn{\mu_2}{mu2}.
#'
#' @name twophase
#' @aliases d2phase p2phase q2phase r2phase h2phase
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param l1 Intensity for transition between phase 1 and phase 2.
#' @param mu1 Intensity for transition from phase 1 to exit.
#' @param mu2 Intensity for transition from phase 2 to exit.
#' @param log logical; if TRUE, return log density or log hazard.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x],
#' otherwise, P[X > x].
#' @return \code{d2phase} gives the density, \code{p2phase} gives the
#' distribution function, \code{q2phase} gives the quantile function,
#' \code{r2phase} generates random deviates, and \code{h2phase} gives the
#' hazard.
#' @section Alternative parameterisation: An individual following this
#' distribution can be seen as coming from a mixture of two populations:
#' 
#' 1) "short stayers" whose mean sojourn time is \eqn{M_1 = }{M1 =
#' 1/(l1+mu1)}\eqn{ 1/(\lambda_1+\mu_1)}{M1 = 1/(l1+mu1)} and sojourn
#' distribution is exponential with rate \eqn{\lambda_1 + \mu_1}{l1+mu1}.
#' 
#' 2) "long stayers" whose mean sojourn time \eqn{M_2 = }{1/(l1+mu1) +
#' 1/mu2}\eqn{ 1/(\lambda_1+\mu_1) + 1/\mu_2}{1/(l1+mu1) + 1/mu2} and sojourn
#' distribution is the sum of two exponentials with rate \eqn{\lambda_1 +
#' }{l1+mu1}\eqn{ \mu_1}{l1+mu1} and \eqn{\mu_2}{mu2} respectively.  The
#' individual is a "long stayer" with probability \eqn{p=\lambda_1/(\lambda_1 +
#' \mu_1)}.
#' 
#' Thus a two-phase distribution can be more intuitively parameterised by the
#' short and long stay means \eqn{M_1 < M_2} and the long stay probability
#' \eqn{p}.  Given these parameters, the transition intensities are
#' \eqn{\lambda_1=p/M_1}{l1=p/M1}, \eqn{\mu_1=(1-p)/M_1}{mu1=(1-p)/M1}, and
#' \eqn{\mu_2=1/(M_2-M_1)}{mu2 = 1/(M2 - M1)}.  This can be useful for choosing
#' intuitively reasonable initial values for procedures to fit these models to
#' data.
#' 
#' The hazard is increasing at least if \eqn{M_2 < 2M_1}{M2 < 2M1}, and also
#' only if \eqn{(M_2 - 2M_1)/(M_2 - M_1) < p}{(M2 - 2M1)/(M2 - M1) < p}.
#' 
#' For increasing hazards with \eqn{\lambda_1 + \mu_1 \leq \mu_2}{l1 + mu1 <=
#' mu2}, the maximum hazard ratio between any time \eqn{t} and time 0 is
#' \eqn{1/(1-p)}.
#' 
#' For increasing hazards with \eqn{\lambda_1 + \mu_1 \geq \mu_2}{l1 + mu1 >=
#' mu2}, the maximum hazard ratio is \eqn{M_1/((1-p)(M_2 - }{M1/((1-p)(M2 -
#' M1))}\eqn{ M_1))}{M1/((1-p)(M2 - M1))}. This is the minimum hazard ratio for
#' decreasing hazards.
#' 
#' % Illustration of hazard ratio at short mean and long mean.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @references C. Dutang, V. Goulet and M. Pigeon (2008). actuar: An R Package
#' for Actuarial Science. Journal of Statistical Software, vol. 25, no. 7,
#' 1-37. URL http://www.jstatsoft.org/v25/i07
#' @keywords distribution
NULL





#' Aortic aneurysm progression data
#' 
#' This dataset contains longitudinal measurements of grades of aortic
#' aneurysms, measured by ultrasound examination of the diameter of the aorta.
#' 
#' 
#' @name aneur
#' @docType data
#' @format A data frame containing 4337 rows, with each row corresponding to an
#' ultrasound scan from one of 838 men over 65 years of age.
#' 
#' \tabular{rll}{ \code{ptnum} \tab (numeric) \tab Patient identification
#' number \cr \code{age} \tab (numeric) \tab Recipient age at examination
#' (years) \cr \code{diam} \tab (numeric) \tab Aortic diameter\cr \code{state}
#' \tab (numeric) \tab State of aneurysm. \cr }
#' 
#' The states represent successive degrees of aneurysm severity, as indicated
#' by the aortic diameter.
#' 
#' \tabular{rll}{ State 1 \tab Aneurysm-free \tab < 30 cm \cr State 2 \tab Mild
#' aneurysm \tab 30-44 cm \cr State 3 \tab Moderate aneurysm \tab 45-54 cm \cr
#' State 4 \tab Severe aneurysm \tab > 55 cm \cr }
#' 
#' 683 of these men were aneurysm-free at age 65 and were re-screened every two
#' years.  The remaining men were aneurysmal at entry and had successive
#' screens with frequency depending on the state of the aneurysm.  Severe
#' aneurysms are repaired by surgery.
#' @references Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W.
#' and Couto, E.  Multi-state Markov models for disease progression with
#' classification error. \emph{The Statistician}, 52(2): 193--209 (2003)
#' 
#' Couto, E. and Duffy, S. W. and Ashton, H. A. and Walker, N. M.  and Myles,
#' J. P. and Scott, R. A. P. and Thompson, S. G. (2002) \emph{Probabilities of
#' progression of aortic aneurysms: estimates and implications for screening
#' policy} Journal of Medical Screening 9(1):40--42
#' @source The Chichester, U.K. randomised controlled trial of screening for
#' abdominal aortic aneurysms by ultrasonography.
#' @keywords datasets
NULL





#' Bronchiolitis obliterans syndrome after lung transplants
#' 
#' A dataset containing histories of bronchiolitis obliterans syndrome (BOS)
#' from lung transplant recipients. BOS is a chronic decline in lung function,
#' often observed after lung transplantation.  The condition is classified into
#' four stages of severity: none, mild, moderate and severe.
#' 
#' The entry time of each patient into each stage of BOS was estimated by
#' clinicians, based on their history of lung function measurements and acute
#' rejection and infection episodes.  BOS is only assumed to occur beyond six
#' months after transplant.  In the first six months the function of each
#' patient's new lung stabilises.  Subsequently BOS is diagnosed by comparing
#' the lung function against the "baseline" value.
#' 
#' The objects \code{bos3} and \code{bos4} contain the same data, but with
#' mild/moderate/severe combined, and moderate/severe combined, to give 3 and
#' 4-state representations respectively.
#' 
#' @name bos
#' @aliases bos bos3 bos4
#' @docType data
#' @format A data frame containing 638 rows, grouped by patient, including
#' histories of 204 patients.  The first observation for each patient is
#' defined to be stage 1, no BOS, at six months after transplant.  Subsequent
#' observations denote the entry times into stages 2, 3, 4, representing mild,
#' moderate and severe BOS respectively, and stage 5, representing death.
#' \tabular{rll}{ \code{ptnum} \tab (numeric) \tab Patient identification
#' number \cr \code{time} \tab (numeric) \tab Months after transplant \cr
#' \code{state} \tab (numeric) \tab BOS state entered at this time \cr }
#' @references Heng. D. et al. (1998).  Bronchiolitis Obliterans Syndrome:
#' Incidence, Natural History, Prognosis, and Risk Factors.  Journal of Heart
#' and Lung Transplantation 17(12)1255--1263.
#' @source Papworth Hospital, U.K.
#' @keywords datasets
NULL





#' Heart transplant monitoring data
#' 
#' A series of approximately yearly angiographic examinations of heart
#' transplant recipients.  The state at each time is a grade of cardiac
#' allograft vasculopathy (CAV), a deterioration of the arterial walls.
#' 
#' 
#' @name cav
#' @docType data
#' @format A data frame containing 2846 rows.  There are 622 patients, the rows
#' are grouped by patient number and ordered by years after transplant, with
#' each row representing an examination and containing additional covariates.
#' 
#' \tabular{rll}{ \code{PTNUM} \tab (numeric) \tab Patient identification
#' number \cr \code{age} \tab (numeric) \tab Recipient age at examination
#' (years) \cr \code{years} \tab (numeric) \tab Examination time (years after
#' transplant)\cr \code{dage} \tab (numeric) \tab Age of heart donor (years)
#' \cr \code{sex} \tab (numeric) \tab sex (0=male, 1=female) \cr \code{pdiag}
#' \tab (factor) \tab Primary diagnosis (reason for transplant) \cr \tab \tab
#' IHD=ischaemic heart disease, IDC=idiopathic dilated cardiomyopathy. \cr
#' \code{cumrej} \tab (numeric) \tab Cumulative number of acute rejection
#' episodes \cr \code{state} \tab (numeric) \tab State at the examination. \cr
#' \tab \tab State 1 represents no CAV, state 2 is mild/moderate CAV \cr \tab
#' \tab and state 3 is severe CAV.  State 4 indicates death.  \cr
#' \code{firstobs} \tab (numeric) \tab 0 = record represents an angiogram or
#' date of death.\cr \tab \tab 1 = record represents transplant (patient's
#' first observation) \cr \code{statemax} \tab (numeric) \tab Maximum observed
#' state so far for this patient (added in version 1.5.1) }
#' @references Sharples, L.D. and Jackson, C.H. and Parameshwar, J. and
#' Wallwork, J. and Large, S.R. (2003). Diagnostic accuracy of coronary
#' angiopathy and risk factors for post-heart-transplant cardiac allograft
#' vasculopathy. Transplantation 76(4):679-82
#' @source Papworth Hospital, U.K.
#' @keywords datasets
NULL





#' Developer documentation: censoring model object
#' 
#' A list giving information about censored states, their labels in the data
#' and what true states they represent.
#' 
#' @name cmodel.object
#' @return \item{ncens}{The number of distinct values used for censored
#' observations in the \code{state} data supplied to \code{\link{msm}}.}
#' \item{censor}{A vector of length \code{ncens}, giving the labels used for
#' censored states in the data.} \item{states}{A vector obtained by
#' \code{unlist()}ing a list with \code{ncens} elements, each giving the set of
#' true states that an observation with this label could be.}
#' \item{index}{Index into \code{states} for the first state corresponding to
#' each \code{censor}, plus an extra \code{length(states)+1}.}
#' @seealso \code{\link{msm.object}}.
NULL





#' Developer documentation: model for covariates on misclassification
#' probabilities
#' 
#' A list representing the model for covariates on misclassification
#' probabilities.
#' 
#' @name ecmodel.object
#' @return \item{npars}{Number of covariate effect parameters.  This is defined
#' as the number of covariates on misclassification (with factors expanded as
#' contrasts) multiplied by the number of allowed misclassifications in the
#' model.  } \item{ndpars}{Number of distinct covariate effect parameters, as
#' \code{npars}, but after any equality constraints have been applied.}
#' \item{ncovs}{Number of covariates on misclassification, with factors
#' expanded as contrasts.} \item{constr}{List of equality constraints on these
#' covariate effects, as supplied in the \code{miscconstraint} argument to
#' \code{\link{msm}}.} \item{covlabels}{Names / labels of these covariates in
#' the model matrix (see \code{\link{model.matrix.msm}}).} \item{inits}{Initial
#' values for these covariate effects, as a vector formed from the
#' \code{misccovinits} list supplied to \code{\link{msm}}.}
#' \item{covmeans}{Means of these covariates in the data (excluding data not
#' required to fit the model, such as observations with missing data in other
#' elements or subjects' last observations).  This includes means of 0/1 factor
#' contrasts as well as continuous covariates (for historic reasons, which may
#' not be sensible).}
#' @seealso \code{\link{msm.object}}.
NULL





#' Developer documentation: misclassification model structure object
#' 
#' A list giving information about the misclassifications assumed in a
#' multi-state model fitted with the \code{ematrix} argument of
#' \code{\link{msm}}.  Returned in a fitted \code{\link{msm}} model object.
#' This information is converted internally to a \code{hmodel} object (see
#' \code{\link{hmodel.object}}) for use in likelihood computations.
#' 
#' @name emodel.object
#' @return \item{nstates}{Number of states (same as \code{qmodel$nstates}).}
#' \item{npars}{Number of allowed misclassifications, equal to
#' \code{sum(imatrix)}.} \item{imatrix}{Indicator matrix for allowed
#' misclassifications.  This has \eqn{(r,s)} entry 1 if misclassification of
#' true state \eqn{r} as observed state \eqn{s} is possible.  diagonal entries
#' are arbitrarily set to 0.} \item{ematrix}{Matrix of initial values for the
#' misclassification probabilities, supplied as the \code{ematrix} argument of
#' \code{\link{msm}}.} \item{inits}{Vector of these initial values, reading
#' across rows of \code{qmatrix} and excluding the diagonal and disallowed
#' transitions.} \item{constr}{Indicators for equality constraints on baseline
#' misclassification probabilities, taken from the \code{econstraint} argument
#' to \code{\link{msm}}, and mapped if necessary to the set (1,2,3,...)}
#' \item{ndpars}{Number of distinct misclassification probabilities, after
#' applying equality constraints.} \item{nipars}{Number of initial state
#' occupancy probabilities being estimated.  This is zero if
#' \code{est.initprobs=FALSE}, otherwise equal to the number of states.}
#' \item{initprobs}{Initial state occupancy probabilities, as supplied to
#' \code{\link{msm}} (initial values before estimation, if
#' \code{est.initprobs=TRUE}.)} \item{est.initprobs}{Are initial state
#' occupancy probabilities estimated (\code{TRUE} or \code{FALSE}), as supplied
#' in the \code{est.initprobs} argument of \code{\link{msm}}.}
#' @seealso \code{\link{msm.object}},\code{\link{qmodel.object}},
#' \code{\link{hmodel.object}}.
NULL





#' FEV1 measurements from lung transplant recipients
#' 
#' A series of measurements of the forced expiratory volume in one second
#' (FEV1) from lung transplant recipients, from six months onwards after their
#' transplant.
#' 
#' A baseline "normal" FEV1 for each individual is calculated using
#' measurements from the first six months after transplant. After six months,
#' as presented in this dataset, FEV1 is expressed as a percentage of the
#' baseline value.
#' 
#' FEV1 is monitored to diagnose bronchiolitis obliterans syndrome (BOS), a
#' long-term lung function decline, thought to be a form of chronic rejection.
#' Acute rejections and infections also affect the lung function in the short
#' term.
#' 
#' @name fev
#' @docType data
#' @format A data frame containing 5896 rows.  There are 204 patients, the rows
#' are grouped by patient number and ordered by days after transplant.  Each
#' row represents an examination and containing an additional covariate.
#' 
#' \tabular{rll}{ \code{ptnum} \tab (numeric) \tab Patient identification
#' number. \cr \code{days} \tab (numeric) \tab Examination time (days after
#' transplant). \cr \code{fev} \tab (numeric) \tab Percentage of baseline FEV1.
#' A code of 999 indicates the patient's date of death. \cr \code{acute} \tab
#' (numeric) \tab 0/1 indicator for whether the patient suffered an acute
#' infection or rejection \cr \tab \tab within 14 days of the visit.  \cr }
#' @references Jackson, C.H. and Sharples, L.D. Hidden Markov models for the
#' onset and progression of bronchiolitis obliterans syndrome in lung
#' transplant recipients \emph{Statistics in Medicine}, 21(1): 113--128 (2002).
#' @source Papworth Hospital, U.K.
#' @keywords datasets
NULL





#' Hidden Markov model constructors
#' 
#' These functions are used to specify the distribution of the response
#' conditionally on the underlying state in a hidden Markov model.  A list of
#' these function calls, with one component for each state, should be used for
#' the \code{hmodel} argument to \code{msm}. The initial values for the
#' parameters of the distribution should be given as arguments. Note the
#' initial values should be supplied as literal values - supplying them as
#' variables is currently not supported.
#' 
#' \code{hmmCat} represents a categorical response distribution on the set
#' \code{1, 2, \dots{}, length(prob)}.  The Markov model with misclassification
#' is an example of this type of model. The categories in this case are (some
#' subset of) the underlying states.
#' 
#' The \code{hmmIdent} distribution is used for underlying states which are
#' observed exactly without error.  For hidden Markov models with multiple
#' outcomes, (see \code{\link{hmmMV}}), the outcome in the data which takes the
#' special \code{hmmIdent} value must be the first of the multiple outcomes.
#' 
#' \code{hmmUnif}, \code{hmmNorm}, \code{hmmLNorm}, \code{hmmExp},
#' \code{hmmGamma}, \code{hmmWeibull}, \code{hmmPois}, \code{hmmBinom},
#' \code{hmmTNorm}, \code{hmmNBinom} and \code{hmmBeta} represent Uniform,
#' Normal, log-Normal, exponential, Gamma, Weibull, Poisson, Binomial,
#' truncated Normal, negative binomial and beta distributions, respectively,
#' with parameterisations the same as the default parameterisations in the
#' corresponding base R distribution functions.
#' 
#' \code{hmmT} is the Student t distribution with general mean \eqn{\mu}{mu},
#' scale \eqn{\sigma}{sigma} and degrees of freedom \code{df}.  The variance is
#' \eqn{\sigma^2 df/(df + 2)}{sigma^2 df/(df + 2)}.  Note the t distribution in
#' base R \code{\link{dt}} is a standardised one with mean 0 and scale 1.
#' These allow any positive (integer or non-integer) \code{df}.  By default,
#' all three parameters, including \code{df}, are estimated when fitting a
#' hidden Markov model, but in practice, \code{df} might need to be fixed for
#' identifiability - this can be done using the \code{fixedpars} argument to
#' \code{\link{msm}}.
#' 
#' The \code{hmmMETNorm} and \code{hmmMEUnif} distributions are truncated
#' Normal and Uniform distributions, but with additional Normal measurement
#' error on the response. These are generalisations of the distributions
#' proposed by Satten and Longini (1996) for modelling the progression of CD4
#' cell counts in monitoring HIV disease.  See \code{\link{medists}} for
#' density, distribution, quantile and random generation functions for these
#' distributions.  See also \code{\link{tnorm}} for density, distribution,
#' quantile and random generation functions for the truncated Normal
#' distribution.
#' 
#' See the PDF manual \file{msm-manual.pdf} in the \file{doc} subdirectory for
#' algebraic definitions of all these distributions.  New hidden Markov model
#' response distributions can be added to \pkg{msm} by following the
#' instructions in Section 2.17.1.
#' 
#' Parameters which can be modelled in terms of covariates, on the scale of a
#' link function, are as follows.
#' 
#' \tabular{ll}{ PARAMETER NAME \tab LINK FUNCTION \cr \code{mean} \tab
#' identity \cr \code{meanlog} \tab identity \cr \code{rate} \tab log \cr
#' \code{scale} \tab log \cr \code{meanerr} \tab identity \cr \code{meanp} \tab
#' logit \cr \code{prob} \tab logit or multinomial logit }
#' 
#' Parameters \code{basecat, lower, upper, size, meanerr} are fixed at their
#' initial values. All other parameters are estimated while fitting the hidden
#' Markov model, unless the appropriate \code{fixedpars} argument is supplied
#' to \code{msm}.
#' 
#' For categorical response distributions \code{(hmmCat)} the outcome
#' probabilities initialized to zero are fixed at zero, and the probability
#' corresponding to \code{basecat} is fixed to one minus the sum of the
#' remaining probabilities.  These remaining probabilities are estimated, and
#' can be modelled in terms of covariates via multinomial logistic regression
#' (relative to \code{basecat}).
#'
#' @name hmm-dists
#' @aliases hmm-dists hmmCat hmmIdent hmmUnif hmmNorm hmmLNorm hmmExp hmmGamma
#' hmmWeibull hmmPois hmmBinom hmmTNorm hmmMETNorm hmmMEUnif hmmNBinom
#' hmmBetaBinom hmmBeta hmmT
#' @param prob (\code{hmmCat}) Vector of probabilities of observing category
#' \code{1, 2, \dots{}, length(prob)} respectively.  Or the probability
#' governing a binomial or negative binomial distribution.
#' @param basecat (\code{hmmCat}) Category which is considered to be the
#' "baseline", so that during estimation, the probabilities are parameterised
#' as probabilities relative to this baseline category. By default, the
#' category with the greatest probability is used as the baseline.
#' @param x (\code{hmmIdent}) Code in the data which denotes the
#' exactly-observed state.
#' @param mean (\code{hmmNorm,hmmLNorm,hmmTNorm}) Mean defining a Normal, or
#' truncated Normal distribution.
#' @param sd (\code{hmmNorm,hmmLNorm,hmmTNorm}) Standard deviation defining a
#' Normal, or truncated Normal distribution.
#' @param meanlog (\code{hmmNorm,hmmLNorm,hmmTNorm}) Mean on the log scale, for
#' a log Normal distribution.
#' @param sdlog (\code{hmmNorm,hmmLNorm,hmmTNorm}) Standard deviation on the
#' log scale, for a log Normal distribution.
#' @param rate (\code{hmmPois,hmmExp,hmmGamma}) Rate of a Poisson, Exponential
#' or Gamma distribution (see \code{\link{dpois}}, \code{\link{dexp}},
#' \code{\link{dgamma}}).
#' @param shape (\code{hmmPois,hmmExp,hmmGamma}) Shape parameter of a Gamma or
#' Weibull distribution (see \code{\link{dgamma}}, \code{\link{dweibull}}).
#' @param shape1,shape2 First and second parameters of a beta distribution (see
#' \code{\link{dbeta}}).
#' @param scale (\code{hmmGamma}) Scale parameter of a Gamma distribution (see
#' \code{\link{dgamma}}), or unstandardised Student t distribution.
#' @param df Degrees of freedom of the Student t distribution.
#' @param size Order of a Binomial distribution (see \code{\link{dbinom}}).
#' @param disp Dispersion parameter of a negative binomial distribution, also
#' called \code{size} or \code{order}.  (see \code{\link{dnbinom}}).
#' @param meanp Mean outcome probability in a beta-binomial distribution
#' @param sdp Standard deviation describing the overdispersion of the outcome
#' probability in a beta-binomial distribution
#' @param lower (\code{hmmUnif,hmmTNorm,hmmMEUnif}) Lower limit for an Uniform
#' or truncated Normal distribution.
#' @param upper (\code{hmmUnif,hmmTNorm,hmmMEUnif}) Upper limit for an Uniform
#' or truncated Normal distribution.
#' @param sderr (\code{hmmMETNorm,hmmUnif}) Standard deviation of the Normal
#' measurement error distribution.
#' @param meanerr (\code{hmmMETNorm,hmmUnif}) Additional shift in the
#' measurement error, fixed to 0 by default.  This may be modelled in terms of
#' covariates.
#' @return Each function returns an object of class \code{hmodel}, which is a
#' list containing information about the model.  The only component which may
#' be useful to end users is \code{r}, a function of one argument \code{n}
#' which returns a random sample of size \code{n} from the given distribution.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}
#' @references Satten, G.A. and Longini, I.M.  Markov chains with measurement
#' error: estimating the 'true' course of a marker of the progression of human
#' immunodeficiency virus disease (with discussion) \emph{Applied Statistics}
#' 45(3): 275-309 (1996).
#' 
#' Jackson, C.H. and Sharples, L.D. Hidden Markov models for the onset and
#' progresison of bronchiolitis obliterans syndrome in lung transplant
#' recipients \emph{Statistics in Medicine}, 21(1): 113--128 (2002).
#' 
#' Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W. and Couto, E.
#' Multi-state Markov models for disease progression with classification error.
#' \emph{The Statistician}, 52(2): 193--209 (2003).
#' @keywords distribution
NULL





#' Developer documentation: hidden Markov model structure object
#' 
#' A list giving information about the models for the outcome data
#' conditionally on the states of a hidden Markov model.  Used in internal
#' computations, and returned in a fitted \code{\link{msm}} model object.
#' 
#' @name hmodel.object
#' @return \item{hidden}{\code{TRUE} for hidden Markov models, \code{FALSE}
#' otherwise.} \item{nstates}{Number of states, the same as
#' \code{qmodel$nstates}.} \item{fitted}{\code{TRUE} if the parameter values in
#' \code{pars} are the maximum likelihood estimates, \code{FALSE} if they are
#' the initial values.} \item{models}{The outcome distribution for each hidden
#' state.  A vector of length \code{nstates} whose \eqn{r}th entry is the index
#' of the state \eqn{r} outcome distributions in the vector of supported
#' distributions.  The vector of supported distributions is given in full by
#' \code{msm:::.msm.HMODELS}: the first few are 1 for categorical outcome, 2
#' for identity, 3 for uniform and 4 for normal. } \item{labels}{String
#' identifying each distribution in \code{models}.} \item{npars}{Vector of
#' length \code{nstates} giving the number of parameters in each outcome
#' distribution, excluding covariate effects.} \item{nipars}{Number of initial
#' state occupancy probabilities being estimated.  This is zero if
#' \code{est.initprobs=FALSE}, otherwise equal to the number of states.}
#' \item{totpars}{Total number of parameters, equal to \code{sum(npars)}. }
#' \item{pars}{A vector of length \code{totpars}, made from concatenating a
#' list of length \code{nstates} whose \eqn{r}th component is vector of the
#' parameters for the state \eqn{r} outcome distribution.  } \item{plabs}{List
#' with the names of the parameters in \code{pars}.} \item{parstate}{A vector
#' of length \code{totpars}, whose \eqn{i}th element is the state corresponding
#' to the \eqn{i}th parameter.} \item{firstpar}{A vector of length
#' \code{nstates} giving the index in \code{pars} of the first parameter for
#' each state.} \item{locpars}{Index in \code{pars} of parameters which can
#' have covariates on them. } \item{initprobs}{Initial state occupancy
#' probabilities, as supplied to \code{\link{msm}} (initial values before
#' estimation, if \code{est.initprobs=TRUE}.)} \item{est.initprobs}{Are initial
#' state occupancy probabilities estimated (\code{TRUE} or \code{FALSE}), as
#' supplied in the \code{est.initprobs} argument of \code{\link{msm}}.}
#' \item{ncovs}{Number of covariate effects per parameter in \code{pars}, with,
#' e.g. factor contrasts expanded.} \item{coveffect}{Vector of covariate
#' effects, of length \code{sum(ncovs)}.} \item{covlabels}{Labels of these
#' effects.} \item{coveffstate}{Vector indicating state corresponding to each
#' element of \code{coveffect}.} \item{ncoveffs}{Number of covariate effects on
#' HMM outcomes, equal to \code{sum(ncovs)}.} \item{nicovs}{Vector of length
#' \code{nstates-1} giving the number of covariate effects on each initial
#' state occupancy probability (log relative to the baseline probability).}
#' \item{icoveffect}{Vector of length \code{sum(nicovs)} giving covariate
#' effects on initial state occupancy probabilities.} \item{nicoveffs}{Number
#' of covariate effects on initial state occupancy probabilities, equal to
#' \code{sum(nicovs)}.} \item{constr}{Constraints on (baseline) hidden Markov
#' model outcome parameters, as supplied in the \code{hconstraint} argument of
#' \code{\link{msm}}, excluding covariate effects, converted to a vector and
#' mapped to the set 1,2,3,\ldots{} if necessary.} \item{covconstr}{Vector of
#' constraints on covariate effects in hidden Markov outcome models, as
#' supplied in the \code{hconstraint} argument of \code{\link{msm}}, excluding
#' baseline parameters, converted to a vector and mapped to the set
#' 1,2,3,\ldots{} if necessary.} \item{ranges}{Matrix of range restrictions for
#' HMM parameters, including those given to the \code{hranges} argument to
#' \code{\link{msm}}.} \item{foundse}{\code{TRUE} if standard errors are
#' available for the estimates.} \item{initpmat}{Matrix of initial state
#' occupancy probabilities with one row for each subject (estimated if
#' \code{est.initprobs=TRUE}).} \item{ci}{Confidence intervals for baseline HMM
#' outcome parameters.} \item{covci}{Confidence intervals for covariate effects
#' in HMM outcome models.}
#' @seealso \code{\link{msm.object}},\code{\link{qmodel.object}},
#' \code{\link{emodel.object}}.
NULL





#' Measurement error distributions
#' 
#' Truncated Normal and Uniform distributions, where the response is also
#' subject to a Normally distributed measurement error.
#' 
#' The normal distribution with measurement error has density
#' 
#' \deqn{ }{(Phi(upper, mu2, sigma3) - Phi(lower, mu2, sigma3)) / (Phi(upper,
#' mean, sd) - Phi(lower, mean, sd)) * phi(x, mean + meanerr, sigma2)}\deqn{
#' \frac{\Phi(u, \mu_2, \sigma_3) - \Phi(l, \mu_2, }{(Phi(upper, mu2, sigma3) -
#' Phi(lower, mu2, sigma3)) / (Phi(upper, mean, sd) - Phi(lower, mean, sd)) *
#' phi(x, mean + meanerr, sigma2)}\deqn{ \sigma_3)}{\Phi(u, \mu_0, \sigma_0) -
#' \Phi(l, \mu_0, \sigma_0)} }{(Phi(upper, mu2, sigma3) - Phi(lower, mu2,
#' sigma3)) / (Phi(upper, mean, sd) - Phi(lower, mean, sd)) * phi(x, mean +
#' meanerr, sigma2)}\deqn{ \phi(x, \mu_0 + \mu_\epsilon, \sigma_2)
#' }{(Phi(upper, mu2, sigma3) - Phi(lower, mu2, sigma3)) / (Phi(upper, mean,
#' sd) - Phi(lower, mean, sd)) * phi(x, mean + meanerr, sigma2)}\deqn{
#' }{(Phi(upper, mu2, sigma3) - Phi(lower, mu2, sigma3)) / (Phi(upper, mean,
#' sd) - Phi(lower, mean, sd)) * phi(x, mean + meanerr, sigma2)} where
#' \deqn{\sigma_2^2 = \sigma_0^2 + \sigma_\epsilon^2,}{sigma2*sigma2 = sd*sd +
#' sderr*sderr,} \deqn{\sigma_3 = \sigma_0 \sigma_\epsilon / \sigma_2,}{sigma3
#' = sd*sderr / sigma2,} \deqn{\mu_2 = (x - \mu_\epsilon) \sigma_0^2 + \mu_0
#' }{mu2 = (x - meanerr)*sd*sd + mean*sderr*sderr,}\deqn{
#' \sigma_\epsilon^2,}{mu2 = (x - meanerr)*sd*sd + mean*sderr*sderr,}
#' 
#' \eqn{\mu_0}{mean} is the mean of the original Normal distribution before
#' truncation, \cr \eqn{\sigma_0}{sd} is the corresponding standard deviation,
#' \cr \eqn{u} is the upper truncation point, \cr \eqn{l} is the lower
#' truncation point, \cr \eqn{\sigma_\epsilon}{sderr} is the standard deviation
#' of the additional measurement error, \cr \eqn{\mu_\epsilon}{meanerr} is the
#' mean of the measurement error (usually 0). \cr \eqn{\phi(x)}{phi(x)} is the
#' density of the corresponding normal distribution, and \cr
#' \eqn{\Phi(x)}{Phi(x)} is the distribution function of the corresponding
#' normal distribution.
#' 
#' The uniform distribution with measurement error has density
#' 
#' \deqn{(\Phi(x, \mu_\epsilon+l, \sigma_\epsilon) - \Phi(x, \mu_\epsilon+u,
#' \sigma_\epsilon)) }{(Phi(x, meanerr+l, sderr) - Phi(x, meanerr+u, sderr)) /
#' (upper - lower)}\deqn{ / (u - l)}{(Phi(x, meanerr+l, sderr) - Phi(x,
#' meanerr+u, sderr)) / (upper - lower)}
#' 
#' These are calculated from the original truncated Normal or Uniform density
#' functions \eqn{f(. | \mu, \sigma, l, u)}{f(. | mu, sd)} as
#' 
#' \deqn{\int f(y }{integral f(y | mu, sd, l, u) phi(x, y + meanerr, sderr)
#' dy}\deqn{ | \mu, \sigma, l, u) \phi(x, y }{integral f(y | mu, sd, l, u)
#' phi(x, y + meanerr, sderr) dy}\deqn{ + \mu_\epsilon, \sigma_\epsilon)
#' dy}{integral f(y | mu, sd, l, u) phi(x, y + meanerr, sderr) dy}
#' 
#' If \code{sderr} and \code{meanerr} are not specified they assume the default
#' values of 0, representing no measurement error variance, and no constant
#' shift in the measurement error, respectively.
#' 
#' Therefore, for example with no other arguments, \code{dmenorm(x)}, is simply
#' equivalent to \code{dtnorm(x)}, which in turn is equivalent to
#' \code{dnorm(x)}.
#' 
#' These distributions were used by Satten and Longini (1996) for CD4 cell
#' counts conditionally on hidden Markov states of HIV infection, and later by
#' Jackson and Sharples (2002) for FEV1 measurements conditionally on states of
#' chronic lung transplant rejection.
#' 
#' These distribution functions are just provided for convenience, and are not
#' optimised for numerical accuracy or speed.  To fit a hidden Markov model
#' with these response distributions, use a \code{\link{hmmMETNorm}} or
#' \code{\link{hmmMEUnif}} constructor. See the \code{\link{hmm-dists}} help
#' page for further details.
#'
#' @name medists
#' @aliases medists dmenorm pmenorm qmenorm rmenorm dmeunif pmeunif qmeunif
#' rmeunif
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @param lower lower truncation point.
#' @param upper upper truncation point.
#' @param sderr Standard deviation of measurement error distribution.
#' @param meanerr Optional shift for the measurement error distribution.
#' @param log,log.p logical; if TRUE, probabilities \eqn{p} are given as
#' \eqn{\log(p)}{log(p)}, or log density is returned.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <=
#' x]}, otherwise, \eqn{P[X > x]}.
#' @return \code{dmenorm}, \code{dmeunif} give the density, \code{pmenorm},
#' \code{pmeunif} give the distribution function, \code{qmenorm},
#' \code{qmeunif} give the quantile function, and \code{rmenorm},
#' \code{rmeunif} generate random deviates, for the Normal and Uniform versions
#' respectively.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{dnorm}}, \code{\link{dunif}}, \code{\link{dtnorm}}
#' @references Satten, G.A. and Longini, I.M.  Markov chains with measurement
#' error: estimating the 'true' course of a marker of the progression of human
#' immunodeficiency virus disease (with discussion) \emph{Applied Statistics}
#' 45(3): 275-309 (1996)
#' 
#' Jackson, C.H. and Sharples, L.D. Hidden Markov models for the onset and
#' progression of bronchiolitis obliterans syndrome in lung transplant
#' recipients \emph{Statistics in Medicine}, 21(1): 113--128 (2002).
#' @keywords distribution
#' @examples
#' 
#' ## what does the distribution look like?
#' x <- seq(50, 90, by=1)
#' plot(x, dnorm(x, 70, 10), type="l", ylim=c(0,0.06)) ## standard Normal
#' lines(x, dtnorm(x, 70, 10, 60, 80), type="l")       ## truncated Normal
#' ## truncated Normal with small measurement error
#' lines(x, dmenorm(x, 70, 10, 60, 80, sderr=3), type="l")
#' 
NULL





#' Fitted msm model objects
#' 
#' The \code{\link{msm}} function returns a list with the following components.
#' These are intended for developers and confident users.  To extract results
#' from fitted model objects, functions such as \code{\link{qmatrix.msm}} or
#' \code{\link{print.msm}} should be used instead.
#' 
#' @name msm.object
#' @return \item{call}{The original call to \code{\link{msm}}, as returned by
#' \code{\link{match.call}}.} \item{Qmatrices}{A list of matrices. The first
#' component, labelled \code{logbaseline}, is a matrix containing the estimated
#' transition intensities on the log scale with any covariates fixed at their
#' means in the data (or at zero, if \code{center=FALSE}). The component
#' labelled \code{baseline} is the equivalent on the untransformed scale. Each
#' remaining component is a matrix giving the linear effects of the labelled
#' covariate on the matrix of log intensities.  To extract an estimated
#' intensity matrix on the natural scale, at an arbitrary combination of
#' covariate values, use the function \code{\link{qmatrix.msm}}.  }
#' \item{QmatricesSE}{The standard error matrices corresponding to
#' \code{Qmatrices}.  } \item{QmatricesL,QmatricesU}{Corresponding lower and
#' upper symmetric confidence limits, of width 0.95 unless specified otherwise
#' by the \code{cl} argument.  } \item{Ematrices}{A list of matrices. The first
#' component, labelled \code{logitbaseline}, is the estimated misclassification
#' probability matrix (expressed as as log odds relative to the probability of
#' the true state) with any covariates fixed at their means in the data (or at
#' zero, if \code{center=FALSE}). The component labelled \code{baseline} is the
#' equivalent on the untransformed scale. Each remaining component is a matrix
#' giving the linear effects of the labelled covariate on the matrix of logit
#' misclassification probabilities.  To extract an estimated misclassification
#' probability matrix on the natural scale, at an arbitrary combination of
#' covariate values, use the function \code{\link{ematrix.msm}}.}
#' \item{EmatricesSE}{The standard error matrices corresponding to
#' \code{Ematrices}.} \item{EmatricesL,EmatricesU}{Corresponding lower and
#' upper symmetric confidence limits, of width 0.95 unless specified otherwise
#' by the \code{cl} argument.  }
#' 
#' \item{minus2loglik}{Minus twice the maximised log-likelihood.}
#' \item{deriv}{Derivatives of the minus twice log-likelihood at its maximum.}
#' 
#' \item{estimates}{Vector of untransformed maximum likelihood estimates
#' returned from \code{\link{optim}}.  Transition intensities are on the log
#' scale and misclassification probabilities are given as log odds relative to
#' the probability of the true state.} \item{estimates.t}{Vector of transformed
#' maximum likelihood estimates with intensities and probabilities on their
#' natural scales.}
#' 
#' \item{fixedpars}{Indices of \code{estimates} which were fixed during the
#' maximum likelihood estimation.} \item{center}{Indicator for whether the
#' estimation was performed with covariates centered on their means in the
#' data.} \item{covmat}{Covariance matrix corresponding to \code{estimates}.}
#' \item{ci}{Matrix of confidence intervals corresponding to
#' \code{estimates.t}}
#' 
#' \item{opt}{Return value from the optimisation routine (such as
#' \code{\link{optim}} or \code{\link{nlm}}), giving information about the
#' results of the optimisation.} \item{foundse}{Logical value indicating
#' whether the Hessian was positive-definite at the supposed maximum of the
#' likelihood.  If not, the covariance matrix of the parameters is unavailable.
#' In these cases the optimisation has probably not converged to a maximum.  }
#' \item{data}{A list giving the data used for the model fit, for use in
#' post-processing.  To extract it, use the methods
#' \code{\link{model.frame.msm}} or \code{\link{model.matrix.msm}}.
#' 
#' The format of this element changed in version 1.4 of \pkg{msm}, so that it
#' now contains a \code{\link{model.frame}} object \code{mf} with all the
#' variables used in the model.  The previous format (an ad-hoc list of vectors
#' and matrices) can be obtained with the function
#' \code{recreate.olddata(msmobject)}, where \code{msmobject} is the object
#' returned by \code{msm}.  } \item{qmodel}{A list of objects representing the
#' transition matrix structure and options for likelihood calculation.  See
#' \code{\link{qmodel.object}} for documentation of the components.}
#' \item{emodel}{A list of objects representing the misclassification model
#' structure, for models specified using the \code{ematrix} argument to
#' \code{\link{msm}}. See \code{\link{emodel.object}}.} \item{qcmodel}{A list
#' of objects representing the model for covariates on transition intensities.
#' See \code{\link{qcmodel.object}}.} \item{ecmodel}{A list of objects
#' representing the model for covariates on transition intensities.  See
#' \code{\link{ecmodel.object}}.} \item{hmodel}{A list of objects representing
#' the hidden Markov model structure. See \code{\link{hmodel.object}}.}
#' \item{cmodel}{A list giving information about censored states.  See
#' \code{\link{cmodel.object}}. } \item{pci}{Cut points for time-varying
#' intensities, as supplied to \code{\link{msm}}, but excluding any that are
#' outside the times observed in the data.} \item{paramdata}{A list giving
#' information about the parameters of the multi-state model.  See
#' \code{\link{paramdata.object}}.} \item{cl}{Confidence interval width, as
#' supplied to \code{\link{msm}}.} \item{covariates}{Formula for covariates on
#' intensities, as supplied to \code{\link{msm}}.}
#' \item{misccovariates}{Formula for covariates on misclassification
#' probabilities, as supplied to \code{\link{msm}}.} \item{hcovariates}{Formula
#' for covariates on hidden Markov model outcomes, as supplied to
#' \code{\link{msm}}.} \item{initcovariates}{Formula for covariates on initial
#' state occupancy probabilities in hidden Markov models, as supplied to
#' \code{\link{msm}}.} \item{sojourn}{ A list as returned by
#' \code{\link{sojourn.msm}}, with components:
#' 
#' \code{mean} = estimated mean sojourn times in the transient states, with
#' covariates fixed at their means (if center=TRUE) or at zero (if
#' center=FALSE).
#' 
#' \code{se} = corresponding standard errors. }
NULL





#' Developer documentation: internal msm parameters object
#' 
#' An object giving information about the parameters of the multi-state model.
#' Used internally during maximum likelihood estimation and arranging results.
#' Returned as the \code{paramdata} component of a fitted \code{\link{msm}}
#' model object.
#' 
#' @name paramdata.object
#' @return \item{inits}{Vector of initial values for distinct parameters which
#' are being estimated.  These have been transformed to the real line (e.g. by
#' log), and exclude parameters being fixed at their initial values, parameters
#' defined to be always fixed (e.g. binomial denominators) and parameters
#' constrained to equal previous ones.} \item{plabs}{Names of parameters in
#' \code{allinits}.} \item{allinits}{Vector of parameter values before
#' estimation, including those which are fixed or constrained to equal other
#' parameters, and transformed to the real line.} \item{hmmpars}{Indices of
#' \code{allinits} which represent baseline parameters of hidden Markov outcome
#' models (thus excluding covariate effects in HMMs and initial state occupancy
#' probabilities). } \item{fixed}{\code{TRUE} if all parameters are fixed,
#' \code{FALSE} otherwise.} \item{fixedpars}{Indices of parameters in
#' \code{allinits} which are fixed, either by definition or as requested by the
#' user in the \code{fixedpars} argument to \code{\link{msm}}. Excludes
#' parameters fixed by constraining to equal other parameters.}
#' \item{notfixed}{Indices of parameters which are not fixed by the definition
#' of \code{fixedpars}.} \item{optpars}{Indices of parameters in
#' \code{allinits} being estimated, thus those included in \code{inits}.}
#' \item{auxpars}{Indices of "auxiliary" parameters which are always fixed, for
#' example, binomial denominators (\code{\link{hmmBinom}}) and the \code{which}
#' parameter in \code{\link{hmmIdent}}.} \item{constr}{Vector of integers, of
#' length \code{npars}, indicating which sets of parameters are constrained to
#' be equal to each other.  If two of these integers are equal the
#' corresponding parameters are equal.  A negative element indicates that
#' parameter is defined to be minus some other parameter (this is used for
#' covariate effects on transition intensities).} \item{npars}{Total number of
#' parameters, equal to \code{length(allinits)}.} \item{nfix}{Number of fixed
#' parameters, equal to \code{length(fixedpars)}. } \item{nopt}{Number of
#' parameters being estimated, equal to \code{length(inits)} and
#' \code{length(optpars)}.} \item{ndup}{Number of parameters defined as
#' duplicates of previous parameters by equality constraints (currently
#' unused).} \item{ranges}{Matrix of defined ranges for each parameter on the
#' natural scale (e.g. 0 to infinity for rate parameters). } \item{opt}{Object
#' returned by the optimisation routine (such as \code{\link{optim}}).}
#' \item{foundse}{\code{TRUE} if standard errors are available after
#' optimisation. If \code{FALSE} the optimisation probably hasn't converged. }
#' \item{lik}{Minus twice the log likelihood at the parameter estimates.}
#' \item{deriv}{Derivatives of the minus twice log likelihood at the parameter
#' estimates, if available. } \item{information}{Corresponding expected
#' information matrix at the parameter estimates, if available.}
#' \item{params}{Vector of parameter values after maximum likelihood
#' estimation, corresponding to \code{allinits}, still on the real-line
#' transformed scale.} \item{covmat}{Covariance matrix corresponding to
#' \code{params}.} \item{ci}{Matrix of confidence intervals corresponding to
#' \code{params}, with nominal coverage (default 0.95) defined by the \code{cl}
#' argument of \code{\link{msm}}. } \item{estimates.t}{Vector of parameter
#' estimates, as \code{params} but with parameters on their natural scales.}
#' @seealso \code{\link{msm.object}}
NULL





#' Exponential distribution with piecewise-constant rate
#' 
#' Density, distribution function, quantile function and random generation for
#' a generalisation of the exponential distribution, in which the rate changes
#' at a series of times.
#' 
#' Consider the exponential distribution with rates \eqn{r_1, \ldots,
#' }{r1,\dots, rn}\eqn{ r_n}{r1,\dots, rn} changing at times \eqn{t_1, \ldots,
#' t_n}{t1, \dots, tn}, with \eqn{t_1 = 0}{t1 = 0}. Suppose \eqn{t_k}{tk} is
#' the maximum \eqn{t_i}{ti} such that \eqn{t_i < x}{ti < x}.  The density of
#' this distribution at \eqn{x > 0} is \eqn{f(x)} for \eqn{k = 1}, and
#' \deqn{\prod_{i=1}^k (1 - F(t_{i} - t_{i-1}, r_i)) f(x - t_{k},
#' r_{k})}{\prod{i=1 \dots k} (1 - F(ti - t{i-1}, r{i-1})) f(x - tk, rk)} for k
#' > 1.
#' 
#' where \eqn{F()} and \eqn{f()} are the distribution and density functions of
#' the standard exponential distribution.
#' 
#' If \code{rate} is of length 1, this is just the standard exponential
#' distribution.  Therefore, for example, \code{dpexp(x)}, with no other
#' arguments, is simply equivalent to \code{dexp(x)}.
#' 
#' Only \code{rpexp} is used in the \code{msm} package, to simulate from Markov
#' processes with piecewise-constant intensities depending on time-dependent
#' covariates.  These functions are merely provided for completion, and are not
#' optimized for numerical stability or speed.
#'
#' @name pexp
#' @aliases pexp dpexp ppexp qpexp rpexp
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param rate vector of rates.
#' @param t vector of the same length as \code{rate}, giving the times at which
#' the rate changes. The values of \code{t} should be in increasing order.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p), or
#' log density is returned.
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x],
#' otherwise, P[X > x].
#' @param start numeric scalar; delayed entry time. The random deviates will be
#' left truncated from this start time.
#' @return \code{dpexp} gives the density, \code{ppexp} gives the distribution
#' function, \code{qpexp} gives the quantile function, and \code{rpexp}
#' generates random deviates.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{dexp}}, \code{\link{sim.msm}}.
#' @keywords distribution
#' @examples
#' 
#' x <- seq(0.1, 50, by=0.1)
#' rate <- c(0.1, 0.2, 0.05, 0.3)
#' t <- c(0, 10, 20, 30)
#' ## standard exponential distribution
#' plot(x, dexp(x, 0.1), type="l")
#' ## distribution with piecewise constant rate
#' lines(x, dpexp(x, rate, t), type="l", lty=2)
#' ## standard exponential distribution
#' plot(x, pexp(x, 0.1), type="l")
#' ## distribution with piecewise constant rate
#' lines(x, ppexp(x, rate, t), type="l", lty=2)
#' 
NULL





#' Psoriatic arthritis data
#' 
#' A series of observations of grades of psoriatic arthritis, as indicated by
#' numbers of damaged joints.
#' 
#' 
#' @name psor
#' @docType data
#' @format A data frame containing 806 observations, representing visits to a
#' psoriatic arthritis (PsA) clinic from 305 patients.  The rows are grouped by
#' patient number and ordered by examination time. Each row represents an
#' examination and contains additional covariates.
#' 
#' \tabular{rll}{ \code{ptnum} \tab (numeric) \tab Patient identification
#' number \cr \code{months} \tab (numeric) \tab Examination time in months \cr
#' \code{state} \tab (numeric) \tab Clinical state of PsA.  Patients in states
#' 1, 2, 3 and 4 \cr \tab \tab have 0, 1 to 4, 5 to 9 and 10 or more damaged
#' joints, \cr \tab \tab respectively.  \cr \code{hieffusn} \tab (numeric) \tab
#' Presence of five or more effusions \cr \code{ollwsdrt} \tab (character) \tab
#' Erythrocyte sedimentation rate of less than 15 mm/h \cr }
#' @references Gladman, D. D. and Farewell, V.T. (1999) Progression in
#' psoriatic arthritis: role of time-varying clinical indicators.  J.
#' Rheumatol. 26(11):2409-13
#' @keywords datasets
#' @examples
#' 
#' ## Four-state progression-only model with high effusion and low
#' ## sedimentation rate as covariates on the progression rates.  High
#' ## effusion is assumed to have the same effect on the 1-2, 2-3, and 3-4
#' ## progression rates, while low sedimentation rate has the same effect
#' ## on the 1-2 and 2-3 intensities, but a different effect on the 3-4. 
#' 
#' data(psor)
#' psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
#' psor.msm <- msm(state ~ months, subject=ptnum, data=psor, 
#'                 qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
#'                 constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)),
#'                 fixedpars=FALSE, control = list(REPORT=1,trace=2), method="BFGS")
#' qmatrix.msm(psor.msm)
#' sojourn.msm(psor.msm)
#' hazard.msm(psor.msm)
#' 
NULL





#' Developer documentation: model for covariates on transition intensities
#' 
#' A list representing the model for covariates on transition intensities
#' 
#' @name qcmodel.object
#' @return \item{npars}{Number of covariate effect parameters.  This is defined
#' as the number of covariates on intensities (with factors expanded as
#' contrasts) multiplied by the number of allowed transitions in the model.
#' 
#' Note if \code{\link{msm}} was called with \code{covariates} set to a list of
#' different covariates for different intensities, then this will include
#' covariate effects that are implicitly defined as zero by this list.  The
#' information in \code{\link[=paramdata.object]{paramdata}} objects can be
#' used to identify wich ones are fixed at zero.
#' 
#' This also includes any \code{timeperiod} covariates in a time-inhomogeneous
#' model defined by the \code{pci} option to \code{\link{msm}}.  }
#' \item{ndpars}{Number of distinct covariate effect parameters, as
#' \code{npars}, but after any equality constraints have been applied.}
#' \item{ncovs}{Number of covariates on intensities, with factors expanded as
#' contrasts.} \item{constr}{List of equality constraints on these covariate
#' effects, as supplied in the \code{constraint} argument to
#' \code{\link{msm}}.} \item{covlabels}{Names / labels of these covariates in
#' the model matrix (see \code{\link{model.matrix.msm}}).} \item{inits}{Initial
#' values for these covariate effects, as a vector formed from the
#' \code{covinits} list supplied to \code{\link{msm}}.} \item{covmeans}{Means
#' of these covariates in the data (excluding data not required to fit the
#' model, such as observations with missing data in other elements or subjects'
#' last observations).  This includes means of 0/1 factor contrasts as well as
#' continuous covariates (for historic reasons, which may not be sensible).}
#' @seealso \code{\link{msm.object}}.
NULL





#' Developer documentation: transition model structure object
#' 
#' A list giving information about the structure of states and allowed
#' transitions in a multi-state model, and options for likelihood calculation.
#' Used in internal computations, and returned in a fitted \code{\link{msm}}
#' model object.
#' 
#' @name qmodel.object
#' @return \item{nstates}{Number of states} \item{iso}{Label for which basic
#' structure the model is isomorphic to in the list of structures for which
#' analytic formulae for the transition probabilities are implemented in the
#' source file \code{src/analyticp.c}.  This list is given by the internal
#' object \code{msm:::.msm.graphs} which is defined and documented in the
#' source file \code{R/constants.R}.
#' 
#' \code{iso} is 0 if the analytic P matrix is not implemented for this
#' structure, or if analytic P matrix calculations are disabled using
#' \code{use.analyticp=FALSE} in the call to \code{\link{msm}}.  }
#' \item{perm}{Permutation required to convert the base isomorphism into the
#' structure of this model. A vector of integers whose \eqn{r}th element is the
#' state number in the base structure representing state \eqn{r} in the current
#' structure.  } \item{qperm}{Inverse permutation: vector whose \eqn{r}th
#' element is the state number in the current structure representing the
#' \eqn{r}th state in the base structure.} \item{npars}{Number of allowed
#' instantaneous transitions, equal to \code{sum(imatrix)}.}
#' \item{imatrix}{Indicator matrix for allowed instantaneous transitions.  This
#' has \eqn{(r,s)} entry 1 if the transition from \eqn{r} to \eqn{s} is
#' permitted in continuous time, and 0 otherwise.  The diagonal entries are
#' arbitrarily set to 0.} \item{qmatrix}{Matrix of initial values for the
#' transition intensities, supplied as the \code{qmatrix} argument of
#' \code{\link{msm}}.} \item{inits}{Vector of these initial values, reading
#' across rows of \code{qmatrix} and excluding the diagonal and disallowed
#' transitions.} \item{constr}{Indicators for equality constraints on baseline
#' intensities, taken from the \code{qconstraint} argument to
#' \code{\link{msm}}, and mapped if necessary to the set (1,2,3,...).}
#' \item{ndpars}{Number of distinct allowed instantaneous transitions, after
#' applying equality constraints.} \item{expm}{Use \pkg{expm} package to
#' calculate matrix exponentials for likelihoods, as supplied to the
#' \code{use.expm} argument of \code{\link{msm}}.  \code{TRUE} or
#' \code{FALSE}.}
#' @seealso \code{\link{msm.object}},\code{\link{emodel.object}},
#' \code{\link{hmodel.object}}.
NULL





#' Summarise a fitted multi-state model
#' 
#' Summary method for fitted \code{\link{msm}} models. This is simply a wrapper
#' around \code{\link{prevalence.msm}} which produces a table of observed and
#' expected state prevalences for each time, and for models with covariates,
#' \code{\link{hazard.msm}} to print hazard ratios with 95\% confidence
#' intervals for covariate effects.
#' 
#' @name summary.msm
#' @aliases summary.msm print.summary.msm
#' @param object A fitted multi-state model object, as returned by
#' \code{\link{msm}}.
#' @param hazard.scale Vector with same elements as number of covariates on
#' transition rates. Corresponds to the increase in each covariate used to
#' calculate its hazard ratio. Defaults to all 1.
#' @param ... Further arguments passed to \code{\link{prevalence.msm}}.
#' @return A list of class \code{summary.msm}, with components:
#' 
#' \item{prevalences}{Output from \code{\link{prevalence.msm}}.}
#' 
#' \item{hazard}{Output from \code{\link{hazard.msm}}.}
#' 
#' \item{hazard.scale}{Value of the \code{hazard.scale} argument.}
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}},\code{\link{prevalence.msm}},
#' \code{\link{hazard.msm}}
#' @keywords models
NULL





#' Truncated Normal distribution
#' 
#' Density, distribution function, quantile function and random generation for
#' the truncated Normal distribution with mean equal to \code{mean} and
#' standard deviation equal to \code{sd} before truncation, and truncated on
#' the interval \code{[lower, upper]}.
#' 
#' The truncated normal distribution has density
#' 
#' \deqn{ f(x, \mu, \sigma) = \phi(x, \mu, \sigma) / (\Phi(u, \mu, \sigma) -
#' \Phi(l, \mu, \sigma)) }{ f(x, mu, sigma) = phi(x, mu, sigma) / (Phi(upper,
#' mu, sigma) - Phi(lower, mu, sigma)) }\deqn{ }{ f(x, mu, sigma) = phi(x, mu,
#' sigma) / (Phi(upper, mu, sigma) - Phi(lower, mu, sigma)) } for \eqn{l <= x
#' <= u}{lower <= x <= upper}, and 0 otherwise.
#' 
#' \eqn{\mu}{mean} is the mean of the original Normal distribution before
#' truncation, \cr \eqn{\sigma}{sd} is the corresponding standard deviation,
#' \cr \eqn{u} is the upper truncation point, \cr \eqn{l} is the lower
#' truncation point, \cr \eqn{\phi(x)}{phi(x)} is the density of the
#' corresponding normal distribution, and \cr \eqn{\Phi(x)}{Phi(x)} is the
#' distribution function of the corresponding normal distribution.
#' 
#' If \code{mean} or \code{sd} are not specified they assume the default values
#' of \code{0} and \code{1}, respectively.
#' 
#' If \code{lower} or \code{upper} are not specified they assume the default
#' values of \code{-Inf} and \code{Inf}, respectively, corresponding to no
#' lower or no upper truncation.
#' 
#' Therefore, for example, \code{dtnorm(x)}, with no other arguments, is simply
#' equivalent to \code{dnorm(x)}.
#' 
#' Only \code{rtnorm} is used in the \code{msm} package, to simulate from
#' hidden Markov models with truncated normal distributions. This uses the
#' rejection sampling algorithms described by Robert (1995).
#' 
#' These functions are merely provided for completion, and are not optimized
#' for numerical stability or speed.  To fit a hidden Markov model with a
#' truncated Normal response distribution, use a \code{\link{hmmTNorm}}
#' constructor. See the \code{\link{hmm-dists}} help page for further details.
#'
#' @name tnorm
#' @aliases tnorm dtnorm ptnorm qtnorm rtnorm
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @param lower lower truncation point.
#' @param upper upper truncation point.
#' @param log logical; if TRUE, return log density or log hazard.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x],
#' otherwise, P[X > x].
#' @return \code{dtnorm} gives the density, \code{ptnorm} gives the
#' distribution function, \code{qtnorm} gives the quantile function, and
#' \code{rtnorm} generates random deviates.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{dnorm}}
#' @references Robert, C. P. Simulation of truncated normal variables.
#' Statistics and Computing (1995) 5, 121--125
#' @keywords distribution
#' @examples
#' 
#' x <- seq(50, 90, by=1)
#' plot(x, dnorm(x, 70, 10), type="l", ylim=c(0,0.06)) ## standard Normal distribution
#' lines(x, dtnorm(x, 70, 10, 60, 80), type="l")       ## truncated Normal distribution
#' 
NULL



