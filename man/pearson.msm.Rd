\name{pearson.msm}
\alias{pearson.msm}
\title{Pearson-type goodness-of-fit test}
\description{Pearson-type goodness-of-fit test for multi-state models
  fitted to panel-observed data.
}
\usage{
pearson.msm(x, transitions=NULL, timegroups=3, intervalgroups=3,
            covgroups=3, groups=NULL, boot=FALSE, B=500,
            next.obstime=NULL, N=100, indep.cens=TRUE,
            maxtimes=NULL, pval=TRUE)
}
\arguments{
  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}.}
  \item{transitions}{This should be an integer vector indicating which
    interval transitions should be grouped together in the contingency
    table.  Its length should be the number of allowed interval
    transitions, excluding transitions from absorbing states
    to absorbing states.

    The allowed interval transitions are the set of pairs of states
    \eqn{(a,b)} for which it is possible to observe \eqn{a} at one time
    and \eqn{b} at any later time.  For example, in a "well-disease-death"
    model with allowed \emph{instantaneous} 1-2, 2-3 transitions, there are 5 allowed
    \emph{interval} transitions. In numerical order, these are 1-1, 1-2,
    1-3, 2-2 and 2-3, excluding absorbing-absorbing transitions.

    Then, to group transitions 1-1,1-2 together, and transitions
    2-2,2-3 together, specify

    \code{transitions = c(1,1,2,3,3)}.

    Only transitions from the same state may be grouped.  By default,
    each interval transition forms a separate group.
  }
  \item{timegroups}{Number of groups based on quantiles of the time since the start of
    the process.}
  \item{intervalgroups}{Number of groups based on quantiles of the time interval between
    observations, within time groups}
  \item{covgroups}{Number of groups based on quantiles of \eqn{\sum_r q_{irr}}{sum_r
      q_{irr}}, where \eqn{q_{irr}} are the diagonal entries of the
    transition intensity matrix for the \emph{i}th transition.  These are
    a function of the covariate effects and the covariate values at the
    \emph{i}th transition: \eqn{q_{irr}} is minus the sum of the
    off-diagonal entries \eqn{q_{rs}^{(0)} exp (\beta_{rs}^T z_i)} on
    the \emph{r}th row.

    Thus \code{covgroups} summarises the impact of covariates at each
    observation, by calculating the overall rate of progression through
    states at that observation.

    For time-inhomogeneous models specified using the \code{pci}
    argument to \code{\link{msm}}, if the only covariate is the time
    period, \code{covgroups} is set to 1, since \code{timegroups}
    ensures that transitions are grouped by time.
  }
  \item{groups}{A vector of arbitrary groups in which to categorise each
    transition. This can be an integer vector or a factor.
    This can be used to diagnose specific areas of poor fit.
    For example, the contingency table might be grouped by arbitrary
    combinations of covariates to detect types of individual for whom
    the model fits poorly.

    The length of \code{groups} should be \code{x$data$n}, the number of
    observations used in the model fit, which is the number of observations
    in the original dataset with any missing values excluded.  The value
    of \code{groups} at observation \eqn{i} is used to categorise the
    transition which \emph{ends} at observation i. Values of
    \code{groups} at the first observation for each subject are
    ignored. }

  \item{boot}{Estimate an "exact" p-value using a parametric bootstrap.

    All objects used in the original call to \code{\link{msm}} which
  produced \code{x}, such as the \code{qmatrix}, should be in the
  working environment, or else an
  \dQuote{object not found} error will be given.  This enables the
  original model to be refitted to the replicate datasets.

  Note that \code{groups} cannot be used with bootstrapping, as the
  simulated observations will not be in the same categories as the
  original observations.
  }

  \item{B}{Number of bootstrap replicates.}

  \item{next.obstime}{This is a vector of length \code{x$data$n} (the
    number of observations used in the model fit) giving the time to the
    next \emph{scheduled} observation following each time point.
    This is only used when times to death are known exactly.

    For individuals who died (entered an absorbing state) before the
    next scheduled observation, and the time of death is known exactly,
    \code{next.obstime} would be \emph{greater} than the observed death time.

    If the individual did not die, and a scheduled observation did
    follow that time point, \code{next.obstime} should just be the same
    as the time to that observation.

    \code{next.obstime} is used to determine a grouping of the time
    interval between observations, which should be based on scheduled
    observations. If exact times to death were used in the grouping, then shorter
    intervals would contain excess deaths, and the goodness-of-fit
    statistic would be biased.

    If \code{next.obstime} is unknown, it is multiply-imputed using a
    product-limit estimate based on the intervals to observations other
    than deaths. The resulting tables of transitions are averaged over
    these imputations.    This may be slow.
  }

  \item{N}{Number of imputations for the estimation of the distribution of
    the next scheduled observation time, when there are exact death
    times.}

  \item{indep.cens}{If \code{TRUE}, then times to censoring are included
    in the estimation of the distribution to the next scheduled
    observation time.  If \code{FALSE}, times to censoring are assumed to
    be systematically different from other observation times. }

  \item{maxtimes}{A vector of length \code{x$data$n}, or a common
    scalar, giving an upper bound for the next scheduled observation
    time.  Used in the multiple imputation when times to death are known exactly.
    If a value greater than \code{maxtimes} is simulated, then the next
    scheduled observation is taken as censored.   This should be
    supplied, if known.  If not supplied, this is taken to be the
    maximum interval occurring in the data, plus one time unit.
    For observations which are not exact death times, this should be
    the time since the previous observation.
  }
  \item{pval}{Calculate a p-value using the improved approximation of
    Titman (2009).  This is
    optional since it is not needed during bootstrapping, and it is
    computationally non-trivial.  Only available currently for non-hidden Markov models for
    panel data without exact death times.  Also not available for
    models with censoring, including time-homogeneous models fitted with
    the \code{pci} option to \code{\link{msm}}.
  }
}
\value{
  A list whose first two elements are contingency tables of observed
  transitions \eqn{O} and expected transitions \eqn{E}, respectively,
  for each combination of groups.  The third element is a table of the
  deviances \eqn{(O - E)^2 / E} multiplied by the sign of \eqn{O - E}.
  If the expected number of transitions is zero then the deviance is zero.
  Entries in the third matrix will be bigger in magnitude for groups for
  which the model fits poorly.  \cr



  \item{\code{"test"}}{the fourth element of the list, is a data frame with one row containing
  the Pearson-type goodness-of-fit test statistic \code{stat}.  The test statistic is the sum of the deviances.
  For panel-observed data without exact death times,
  misclassification or censored observations, \code{p} is the
  p-value for the test statistic calculated using the improved
  approximation of Titman (2009).

  For these models, for comparison with older versions of the package,
  \code{test} also presents \code{p.lower} and \code{p.upper}, which are
  theoretical lower and upper limits for the p-value of the test
  statistic, based on \ifelse{latex}{\eqn{\chi^2}}{chi-squared} distributions with \code{df.lower} and
  \code{df.upper} degrees of freedom, respectively.  \code{df.upper} is
  the number of independent cells in the contingency table, and
  \code{df.lower} is \code{df.upper} minus the number of estimated
  parameters in the model.}

  \item{\code{"intervalq"}}{(not printed by default) contains the
  definition of the grouping of the intervals between
  observations.  These groups are defined by quantiles within the groups
  corresponding to the time since the start of the process.}

  \item{\code{"sim"}}{If there are exact death times, this contains
  simulations of the contingency tables and test statistics for each
  imputation of the next scheduled sampling time.  These are averaged
  over to produce the presented tables and test statistic. This element
  is not printed by default.

  With exact death times, the null variance of the test statistic (formed
  by taking mean of simulated test statistics) is less than twice the
  mean (Titman, 2008), and the null distribution is not \ifelse{latex}{\eqn{\chi^2}}{chi-squared}.
  In this case, \code{p.upper} is an upper limit for the true asymptotic
  p-value, but \code{p.lower} is not a lower limit, and is not presented.}

  \item{\code{"boot"}}{If the bootstrap has been used, the element will contain
  the bootstrap replicates of the test statistics (not printed by
  default).}

  \item{\code{"lambda"}}{If the Titman (2009) p-value has been
  calculated, this contains the weights defining the null distribution
  of the test statistic as a weighted sum of \ifelse{latex}{\eqn{\chi^2_1}}{chi-squared(1)} random variables (not printed by default).}
}
\details{
  This method (Aguirre-Hernandez and Farewell, 2002) is intended for
  data which represent observations of the process at arbitrary times
  ("snapshots", or "panel-observed" data). For data which represent the
  exact transition times of the process, \code{\link{prevalence.msm}}
  can be used to assess fit, though without a formal test.

  When times of death are known exactly, states are
  misclassified, or an individual's final observation is a censored
  state, the modification by Titman and Sharples (2008) is
  used. The only form of censoring supported is a state at the end
  of an individual's series which represents an unknown transient state
  (i.e. the individual is only known to be alive at this time). Other
  types of censoring are omitted from the data before performing the
  test.

  See the references for further details of the methods.
  The method used for censored states is a modification of the method
  in the appendix to Titman and Sharples (2008), described at
  \url{http://www.mrc-bsu.cam.ac.uk/wp-content/uploads/robustcensoring.pdf}
  (Titman, 2007).

  Groupings of the time since initiation, the time interval and the
  impact of covariates are based on equally-spaced quantiles.
  The number of groups should be chosen that there are not many cells
  with small expected numbers of transitions, since the deviance
  statistic will be unstable for sparse contingency tables.
  Ideally, the expected numbers of transitions in each cell of the
  table should be no less than about 5.
  Conversely, the power of the test is reduced if there are too few
  groups. Therefore, some sensitivity analysis of the test results to
  the grouping is advisable.

  Saved model objects fitted with previous versions of R (versions less
  than 1.2) will need to be refitted under the current R for use with
  \code{pearson.msm}.

}
\references{
  Aguirre-Hernandez, R. and Farewell, V. (2002) A Pearson-type
  goodness-of-fit test for stationary and time-continuous Markov
  regression models. \emph{Statistics in Medicine} 21:1899-1911.

  Titman, A. and Sharples, L. (2008) A general goodness-of-fit test for
  Markov and hidden Markov models. \emph{Statistics in Medicine} 27(12):2177-2195

  Titman, A. (2009) Computation of the asymptotic null distribution of
  goodness-of-fit tests for multi-state models. \emph{Lifetime Data
  Analysis} 15(4):519-533.

  Titman, A. (2008) Model diagnostics in multi-state models of biological
  systems. PhD thesis, University of Cambridge.
}
\seealso{
  \code{\link{msm}}, \code{\link{prevalence.msm}}, \code{\link{scoreresid.msm}},
}
\examples{
psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
psor.msm <- msm(state ~ months, subject=ptnum, data=psor,
                qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
pearson.msm(psor.msm, timegroups=2, intervalgroups=2, covgroups=2)
# More 1-2, 1-3 and 1-4 observations than expected in shorter time
# intervals - the model fits poorly.
# A random effects model might accommodate such fast progressors.
}
\author{
  Andrew Titman \email{a.titman@lancaster.ac.uk},
  Chris Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}
}
\keyword{models}
