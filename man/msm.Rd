\name{msm}
\title{Multi-state Markov and hidden Markov models in continuous time}
\alias{msm}
\description{
  Fit a continuous-time Markov or hidden Markov multi-state model
  by maximum likelihood. Observations of the process
  can be made at arbitrary times, or the exact times of
  transition between states can be known.
  Covariates can be fitted to the Markov chain transition intensities or
  to the hidden Markov observation process.
}
\usage{
msm ( formula, subject=NULL, data = list(), qmatrix, gen.inits = FALSE,
      ematrix=NULL, hmodel=NULL, obstype=NULL, obstrue=NULL,
      covariates = NULL, covinits = NULL, constraint = NULL,
      misccovariates = NULL, misccovinits = NULL, miscconstraint = NULL,
      hcovariates = NULL, hcovinits = NULL, hconstraint = NULL, hranges=NULL,
      qconstraint=NULL, econstraint=NULL, initprobs = NULL,
      est.initprobs=FALSE, initcovariates = NULL, initcovinits = NULL,
      deathexact = NULL, death=NULL, exacttimes = FALSE, censor=NULL,
      censor.states=NULL, pci=NULL, phase.states=NULL, phase.inits=NULL,
      cl = 0.95, fixedpars = NULL, center=TRUE,
      opt.method="optim", hessian=NULL, use.deriv=TRUE,
      use.expm=TRUE, analyticp=TRUE, na.action=na.omit, ... )
}
\arguments{
  \item{formula}{ A formula giving the vectors containing the observed
    states and the corresponding observation times. For example,

    \code{state ~ time}

    Observed states should be numeric variables
    in the set \code{1, \dots, n}, where
    \code{n} is the number of states.  Factors are allowed only if their
    levels are called \code{"1", \dots, "n"}.

    The times can indicate different types of observation scheme, so be
    careful to choose the correct \code{obstype}.

    For hidden Markov models, \code{state} refers to the outcome
    variable, which need not be a discrete state.   It may also be a
    matrix, giving multiple observations at each time (see
    \code{\link{hmmMV}}). 
  }

  \item{subject}{Vector of subject identification numbers for the data
    specified by \code{formula}. If missing, then all observations
    are assumed to be on the same subject. These must be sorted so that
    all observations on the same subject are adjacent.}

  \item{data}{Optional data frame in which to interpret the variables
    supplied in \code{formula}, \code{subject}, \code{covariates},
    \code{misccovariates}, \code{hcovariates}, \code{obstype} and
    \code{obstrue}.}

  \item{qmatrix}{Matrix which indicates the allowed transitions in
    the continuous-time Markov chain, and optionally also the initial
    values of those transitions. If an instantaneous transition is not allowed  from state
    \eqn{r} to state \eqn{s}, then \code{qmatrix} should have
    \eqn{(r,s)} entry 0, otherwise it should be non-zero.

    If supplying initial values yourself, then the non-zero entries
    should be those values.  If using \code{gen.inits=TRUE} then the
    non-zero entries can be anything you like (conventionally 1).
    Any diagonal entry of \code{qmatrix} is ignored, as it is
    constrained to be equal to minus the sum of the rest of the row.

    For example,\cr

    \code{
    rbind(
    c( 0, 0.1, 0.01 ),
    c( 0.1, 0, 0.2 ),
    c( 0, 0, 0 )
    )
  }\cr

    represents a 'health - disease - death' model, with initial transition intensities
    0.1 from health to disease, 0.01 from health to death, 0.1 from disease to health,
    and 0.2 from disease to death.

    If the states represent ordered levels of severity of a disease,
    then this matrix should usually only allow transitions between
    adjacent states.  For example, if someone was observed in state 1
    ("mild") at their first observation, followed by state 3 ("severe")
    at their second observation, they are assumed to have passed through
    state 2 ("moderate") in between, and the 1,3 entry of \code{qmatrix}
    should be zero.

    The initial intensities given here are with any covariates set to
    their means in the data (or set to zero, if \code{center =
    FALSE}). If any intensities are constrained to be equal using
    \code{qconstraint}, then the initial value is taken from the first
    of these (reading across rows).  }

  \item{gen.inits}{If \code{TRUE}, then initial values for the
    transition intensities are generated automatically using the method
    in \code{\link{crudeinits.msm}}. The non-zero
    entries of the supplied \code{qmatrix} are assumed to indicate the
    allowed transitions of the model.   This is not available for hidden
    Markov models, including models with misclassified states.}

  \item{ematrix}{
    If misclassification between states is to be modelled, this should
    be a matrix of initial values for the misclassification probabilities.
    The rows represent underlying states, and the columns represent
    observed states.
    If an observation of state \eqn{s} is not possible when the subject
    occupies underlying state \eqn{r}, then \code{ematrix} should have
    \eqn{(r,s)} entry 0.  Otherwise \code{ematrix} should
    have \eqn{(r,s)} entry corresponding to the probability of
    observing \eqn{s} conditionally on occupying true state \eqn{r}. The
    diagonal of \code{ematrix} is ignored, as rows are constrained to
    sum to 1.  For example, \cr

    \code{
      rbind(
    c( 0, 0.1, 0 ),
    c( 0.1, 0, 0.1 ),
    c( 0, 0.1, 0 )
    )
  }\cr

    represents a model in which misclassifications are only permitted
    between adjacent states.

    If any probabilities are constrained to be equal using
    \code{econstraint}, then the initial value is
    taken from the first of these (reading across rows).

    For an alternative way of specifying misclassification models, see
    \code{hmodel}.
  }

  \item{hmodel}{
    Specification of the hidden Markov model (HMM).  This should be a
    list of return values from HMM constructor functions.  Each element
    of the list corresponds to the outcome model conditionally on the
    corresponding underlying state.  Univariate constructors are
    described in the\code{\link{hmm-dists}} help page.  These may also
    be grouped together to specify a multivariate HMM with a set of
    conditionally independent univariate outcomes at each time, as
    described in \code{\link{hmmMV}}.

    For example, consider a three-state hidden Markov model.  Suppose
    the observations in underlying state 1 are
    generated from a Normal distribution with mean 100 and standard
    deviation 16, while observations in underlying state 2 are Normal
    with mean 54 and standard deviation 18. Observations in state 3,
    representing death, are exactly observed, and coded as 999 in the
    data.  This model is specified as

    \code{hmodel = list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))}

    The mean and standard deviation parameters are estimated starting
    from these initial values. If multiple parameters are constrained
    to be equal using \code{hconstraint}, then the initial value is
    taken from the value given on the first occasion that
    parameter appears in \code{hmodel}.

    See the \code{\link{hmm-dists}} help page
    for details of the constructor functions for each univariate distribution.

    A misclassification model, that is, a hidden Markov model where
    the outcomes are misclassified observations of the underlying
    states, can either be specified using a list of \code{\link{hmmCat}}
    or \code{\link{hmmIdent}} objects, or by using an \code{ematrix}.

    For example, \cr

     \code{
      ematrix = rbind(
    c( 0, 0.1, 0, 0 ),
    c( 0.1, 0, 0.1, 0 ),
    c( 0, 0.1, 0, 0),
    c( 0, 0, 0, 0)
    )
  }\cr

    is equivalent to \cr

    \code{hmodel = list(
    hmmCat(prob=c(0.9, 0.1, 0, 0)),
    hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
    hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent())
    }\cr

  }

  \item{obstype}{A vector specifying the observation scheme for each row
    of the data. This can be included in the data frame \code{data}
    along with the  state, time, subject IDs and covariates.  Its
    elements should be either 1, 2 or 3, meaning as follows:

    \describe{
      \item{1}{An observation of the process at an arbitrary time (a "snapshot"
	of the process, or "panel-observed" data). The states are
	unknown between observation times.}
      \item{2}{An exact transition time, with the state at the previous
	observation retained until the current observation.  An observation may
    represent a transition to a different state or a repeated observation of
    the same state (e.g. at the end of follow-up).  Note that if all
    transition times are known, more flexible models could be fitted
    with packages other than \pkg{msm} - see the note under \code{exacttimes}.

    Note also that if the previous state was censored using
    \code{censor}, for example known only to be state 1 or state 2, then
    \code{obstype} 2 means that either state 1 is retained or state 2 is
    retained until the current observation - this does not allow for a
    change of state in the middle of the observation interval.
    }
      \item{3}{An exact transition time, but the state at the instant before
	entering this state is unknown. A common example is death times in
	studies of chronic diseases.}
    }
    If \code{obstype} is not specified, this defaults to all 1. If
    \code{obstype} is a single number, all observations are assumed to
    be of this type.  The obstype value for the first observation
    from each subject is not used.

    This is a generalisation of the \code{deathexact} and \code{exacttimes}
    arguments to allow different schemes per observation.
    \code{obstype} overrides both \code{deathexact} and \code{exacttimes}.

    \code{exacttimes=TRUE} specifies that all observations are of
    obstype 2.

    \code{deathexact = death.states} specifies that all observations of
    \code{death.states} are of type 3.  \code{deathexact =
      TRUE} specifies that all observations in the final absorbing state
    are of type 3.
  }
  
  \item{obstrue}{In misclassification models specified with
  \code{ematrix}, \code{obstrue} is a vector of logicals (\code{TRUE} or
  \code{FALSE}) or numerics (1 or 0) specifying which observations
  (\code{TRUE}, 1) are observations of the underlying state without
  error, and which (\code{FALSE}, 0) are realisations of a hidden Markov
  model.

  In HMMs specified with \code{hmodel}, where the hidden state is known
  at some times, if \code{obstrue} is supplied it is assumed to contain
  the actual true state data.  Elements of \code{obstrue} at times when
  the hidden state is unknown are set to \code{NA}.  This allows the
  information from HMM outcomes generated conditionally on the known
  state to be included in the model, thus improving the estimation of
  the HMM outcome distributions.
  
  In HMMs where there are also censored states, \code{obstrue} should be
  set to 1 for observed states which are \emph{censored} but not
  \emph{misclassified}.
  }

  \item{covariates}{A formula or a list of formulae representing the
    covariates on the transition intensities via a log-linear model. If
    a single formula is supplied, like

    \code{covariates = ~ age + sex + treatment}

    then these covariates are assumed to apply to all intensities.  If a
    named list is supplied, then this defines a potentially different
    model for each named intensity.  For example,

    \code{covariates = list("1-2" = ~ age, "2-3" = ~ age + treatment)}

    specifies an age effect on the state 1 - state 2 transition,
    additive age and treatment effects on the state 2 - state 3
    transition, but no covariates on any other transitions that are
    allowed by the \code{qmatrix}.

    If covariates are time dependent, they are assumed to be constant in
    between the times they are observed, and the transition probability
    between a pair of times \eqn{(t1, t2)} is assumed to depend on the
    covariate value at \eqn{t1}.   
  }

  \item{covinits}{Initial values for log-linear effects of covariates on the
    transition intensities. This should be a named list with each element
    corresponding to a covariate.  A single element contains the initial
    values for that covariate on each transition intensity, reading
    across the rows in order.  For a pair of effects constrained to be
    equal, the initial value for the first of the two effects is used.

    For example, for a model with the above \code{qmatrix} and age and sex
    covariates, the following initialises all covariate effects to zero
    apart from the age effect on the 2-1 transition, and the sex
    effect on the 1-3 transition.
    \code{ covinits = list(sex=c(0, 0, 0.1, 0), age=c(0, 0.1, 0, 0))}

    For factor covariates, name each level by concatenating the name of
    the covariate with the level name, quoting if necessary. For example,
    for a covariate  \code{agegroup} with three levels \code{0-15,
      15-60, 60-}, use something like

    \code{ covinits = list("agegroup15-60"=c(0, 0.1, 0, 0),
      "agegroup60-"=c(0.1, 0.1, 0, 0))}

    If not specified or wrongly specified, initial values are assumed
    to be zero.
    }

    \item{constraint}{ A list of one numeric vector for each named covariate. The
  vector indicates which covariate effects on intensities are
  constrained to be equal. Take, for example, a model with five
  transition intensities and two covariates. Specifying\cr

    \code{constraint = list (age = c(1,1,1,2,2),  treatment = c(1,2,3,4,5))}\cr

    constrains the effect of age to be equal for the first three
    intensities, and equal for the fourth and fifth. The effect of
    treatment is assumed to be different for each intensity. Any vector of
    increasing numbers can be used as indicators. The intensity parameters are
    assumed to be ordered by reading across the rows of the
    transition matrix, starting at the first row, ignoring the
    diagonals.

    Negative elements of the vector can be used to indicate that particular
    covariate effects are constrained to be equal to minus some other
    effects.  For example:

    \code{constraint = list (age = c(-1,1,1,2,-2), treatment = c(1,2,3,4,5)) }\cr

    constrains the second and third age effects to be equal, the first effect
    to be minus the second, and the fifth age effect to be minus the
    fourth.  For example, it may be realisitic that the effect of a
    covariate on the "reverse" transition rate from state 2 to state 1
    is minus the effect on the "forward" transition rate, state 1 to state
    2.  Note that it is not possible to specify exactly which of the
    covariate effects are constrained to be positive and which negative.
    The maximum likelihood estimation chooses the combination of signs
    which has the higher likelihood.

    For categorical covariates, defined using \code{factor(covname)},
    specify constraints as follows:\cr

    \code{list(..., covnameVALUE1 = c(...), covnameVALUE2 = c(...), ...)}\cr

    where \code{covname} is the name of the factor, and \code{VALUE1},
    \code{VALUE2}, ... are the labels of the factor levels (usually
    excluding the baseline, if using the default contrasts).

    Make sure the \code{contrasts} option is set appropriately, for
    example, the default

    \code{options(contrasts=c(contr.treatment, contr.poly))}

    sets the first (baseline) level of unordered factors to zero,
    then the baseline level is ignored in this specification.

    To assume no covariate effect on a certain transition, use the
    \code{fixedpars} argument to fix it at its initial value (which
    is zero by default) during the optimisation.

  }

  \item{misccovariates}{A formula representing the covariates on the
    misclassification probabilities, analogously to \code{covariates},
    via multinomial logistic regression. Only used if the model is
    specified using \code{ematrix}, rather than \code{hmodel}.

    This must be a single formula - lists are not supported, unlike
    \code{covariates}.  If a different model on each probability is
    required, include all covariates in this formula, and use
    \code{fixedpars} to fix some of their effects (for particular
    probabilities) at their default initial values of zero.
  }

  \item{misccovinits}{Initial values for the covariates on the
    misclassification probabilities, defined in the same way as \code{covinits}.
    Only used if the model is specified using \code{ematrix}.
  }

  \item{miscconstraint}{A list of one vector for each named covariate on
    misclassification probabilities. The vector indicates which
    covariate effects on misclassification probabilities are
    constrained to be equal, analogously to \code{constraint}.
    Only used if the model is specified using \code{ematrix}.
  }

  \item{hcovariates}{
    List of formulae the same length as \code{hmodel}, defining any covariates
    governing the hidden Markov outcome models.  The covariates operate
    on a suitably link-transformed linear scale, for example, log scale
    for a Poisson outcome model. If there are no covariates for a certain
    hidden state, then insert a NULL in the corresponding place in the
    list.  For example,
    \code{hcovariates = list(~acute + age, ~acute, NULL).}
    }

  \item{hcovinits}{
    Initial values for the hidden Markov model covariate effects. A list of the
    same length as \code{hcovariates}. Each element is a vector with
    initial values for the effect of each covariate on that state.  For example, the
    above \code{hcovariates} can be initialised with
    \code{hcovariates = list(c(-8, 0), -8, NULL)}. Initial values must
    be given for all or no covariates, if none are given these are all
    set to zero.  The initial value given in the \code{hmodel}
    constructor function for the corresponding baseline parameter is
    interpreted as the value of that parameter with any covariates fixed
    to their means in the data.  If multiple effects are constrained
    to be equal using \code{hconstraint}, then the
    initial value is taken from the first of the multiple initial values
    supplied.
    }

  \item{hconstraint}{
    A named list. Each element is a vector of constraints on the named hidden
    Markov model parameter. The vector has length equal to the number of
    times that class of parameter appears in the whole model.

    For example consider the three-state hidden Markov model described
    above, with normally-distributed outcomes for states 1 and 2.
    To constrain the outcome variance to be equal for states 1 and 2,
    and to also constrain the effect of \code{acute} on the outcome mean
    to be equal for states 1 and 2, specify

    \code{hconstraint = list(sd = c(1,1), acute=c(1,1))}

    Note this excludes initial state occupancy probabilities and
    covariate effects on those probabilities, which cannot be constrained.
    }

    \item{hranges}{Range constraints for hidden Markov model parameters.
      Supplied as a named list, with each element corresponding to the
      named hidden Markov model parameter.  This element is itself a
      list with two elements, vectors named "lower" and "upper".  These
      vectors each have length equal to the number of times that class
      of parameter appears in the whole model, and give the
      corresponding mininum amd maximum allowable values for that
      parameter.  Maximum likelihood estimation is performed with these
      parameters constrained in these ranges (through a log or
      logit-type  transformation).
      Lower bounds of \code{-Inf} and upper bounds of \code{Inf} can
      be given if the parameter is unbounded above or below.

      For example, in the three-state model above, to constrain the mean
      for state 1 to be between 0 and 6, and the mean of state 2 to be
      between 7 and 12, supply

      \code{hranges=list(mean=list(lower=c(0, 7), upper=c(6, 12)))}

      These default to the natural ranges, e.g. the positive real line
      for variance parameters, and [0,1] for probabilities.  Therefore
      \code{hranges} need not be specified for such parameters unless an
      even stricter constraint is desired.  If only one limit is
      supplied for a parameter, only the first occurrence of that
      parameter is constrained.

      Initial values should be strictly within any ranges, and not on
      the range boundary, otherwise optimisation will fail with a
      "non-finite value" error.
    }

  \item{qconstraint}{A vector of indicators specifying which baseline
    transition intensities are equal. For example,

    \code{qconstraint = c(1,2,3,3)}

    constrains the third and fourth intensities to be equal, in a model
    with four allowed instantaneous transitions.   When there are
    covariates on the intensities and \code{center=TRUE} (the default),
    \code{qconstraint} is applied to the intensities with covariates
    taking the values of the means in the data.  When \code{center=FALSE},
    \code{qconstraint} is applied to the intensities with covariates set
    to zero.
  }

  \item{econstraint}{A similar vector of indicators specifying which
    baseline misclassification probabilities are constrained to be
    equal.      Only used if the model is specified using \code{ematrix}, rather
    than \code{hmodel}.  }

  \item{initprobs}{Only used in hidden Markov models.  Underlying state
    occupancy probabilities at each subject's first observation.  Can
    either be a vector of \eqn{nstates} elements with common probabilities
    to all subjects, or a \eqn{nsubjects} by \eqn{nstates} matrix of
    subject-specific probabilities.  This refers to observations after
    missing data and subjects with only one observation have been excluded.

    If these are estimated (see \code{est.initprobs}), then this
    represents an initial value, and defaults to equal probability for
    each state.  Otherwise this defaults to \code{c(1, rep(0,
    nstates-1))}, that is, in state 1 with a probability of 1.  Scaled
    to sum to 1 if necessary.  The state 1 occupancy probability should
    be non-zero.  }

  \item{est.initprobs}{Only used in hidden Markov models.
    If \code{TRUE}, then the underlying state
    occupancy probabilities at the first observation will be estimated,
    starting from a vector of initial values supplied in the \code{initprobs}
    argument.  Structural zeroes are allowed: if any of these initial
    values are zero they will be fixed at zero during optimisation,
    even if \code{est.initprobs=TRUE}, and no covariate effects on
    them are estimated.  The exception is state 1, which should have
    non-zero occupancy probability.

    Note that the free parameters during this estimation exclude the state 1
    occupancy probability, which is fixed at one minus the sum of the
    other probabilities.
  }

  \item{initcovariates}{ Formula representing covariates on the initial
    state occupancy probabilities, via multinomial logistic regression.  The
    linear effects of these covariates, observed at the individual's first
    observation time, operate on the log ratio of the state \eqn{r}
    occupancy probability to the state 1 occupancy probability, for each
    \eqn{r = 2} to the number of states.  Thus the state 1 occupancy
    probability should be non-zero. If \code{est.initprobs} is
    \code{TRUE}, these effects are estimated starting from their initial
    values.  If \code{est.initprobs} is \code{FALSE}, these effects are fixed at
    theit initial values.}

  \item{initcovinits}{ Initial values for the covariate effects
    \code{initcovariates}.  A named list with each element corresponding
    to a covariate, as in \code{covinits}. Each element is a vector with
    (1 - number of states) elements, containing the initial values for
    the linear effect of that covariate on the log odds of that state
    relative to state 1, from state 2 to the final state.  If
    \code{initcovinits} is not specified, all covariate effects are
    initialised to zero.}

  \item{deathexact}{Vector of indices of absorbing states whose time of entry
    is known exactly, but the individual is assumed to be in an unknown
    transient state ("alive") at the previous instant.  This is the
    usual situation for times of death in chronic disease monitoring
    data.  For example, if you specify \code{deathexact = c(4, 5)} then
    states 4 and 5 are assumed to be exactly-observed death states.

    See the \code{obstype} argument.  States of this kind correspond to
    \code{obstype=3}. \code{deathexact = TRUE} indicates that the final
    absorbing state is of this kind, and \code{deathexact = FALSE} or
    \code{deathexact = NULL} (the default) indicates that there is no
    state of this kind.

    The \code{deathexact} argument is overridden by \code{obstype} or
    \code{exacttimes}.

    Note that you do not always supply a \code{deathexact} argument, even if
    there are states that correspond to deaths, because they do
    not necessarily have
    \code{obstype=3}.  If the state is known between the time of death
    and the previous observation, then you should specify
    \code{obstype=2} for the death times, or \code{exacttimes=TRUE} if
    the state is known at all times, and the \code{deathexact} argument is
    ignored.
  }

  \item{death}{Old name for the \code{deathexact} argument. 
    Overridden by \code{deathexact} if both are supplied. Deprecated.}

  \item{censor}{ A state, or vector of states, which indicates
    censoring.  Censoring means that the observed state is known only to be
    one of a particular set of states. For example, \code{censor=999}
    indicates that all observations of \code{999} in the vector of observed
    states are censored states.  By default, this means that the true
    state could have been any of the transient (non-absorbing) states. To
    specify corresponding true states explicitly, use a \code{censor.states}
    argument.

    Note that in contrast to the usual terminology of survival
    analysis, here it is the \emph{state} which is considered to be censored,
    rather than the \emph{event time}.  If at the end of a study, an individual
    has not died, but their true state is \emph{known}, then \code{censor} is
    unnecessary, since the standard multi-state model likelihood is
    applicable.   Also a "censored" state here can be at any time, not
    just at the end.

    Note in particular that general time-inhomogeneous Markov models with
    piecewise constant transition
    intensities can be constructed using the \code{censor}
    facility. If the true state is unknown on occasions when a
    piecewise constant covariate is known to change, then censored
    states can be inserted in the data on those occasions.
    The covariate may represent time itself, in which case the \code{pci} option to msm can be used
    to perform this trick automatically, or some other time-dependent
    variable.
  }

  \item{censor.states}{
    Specifies the underlying states which censored observations can
    represent. If \code{censor} is a single number (the default) this
    can be a vector, or a list with one element.  If \code{censor} is a
    vector with more than one element, this should be a list, with each
    element a vector corresponding to the equivalent element of
    \code{censor}.  For example

    \code{censor = c(99, 999), censor.states = list(c(2,3), c(3,4))}

    means that observations coded 99 represent either state 2 or state
    3, while observations coded 999 are really either state 3 or
    state 4.
  }

  \item{pci}{Model for piecewise-constant intensities.  Vector of cut
     points defining the times, since the start of the process, at
     which intensities change for all subjects.  For example

     \code{pci = c(5, 10)}

     specifies that the intensity changes at time points 5 and 10.
     This will automatically construct a model with a categorical
     (factor) covariate called \code{timeperiod}, with levels
     \code{"[-Inf,5)"},
     \code{"[5,10)"} and \code{"[10,Inf)"}, where
     the first
     level is the baseline.  This covariate defines the time period in
     which the observation was made.  Initial values and
     constraints on covariate effects are specified the same way as for
     a model with a covariate of this name, for example,

     \code{covinits = list("timeperiod[5,10)"=c(0.1,0.1),
                           "timeperiod[10,Inf)"=c(0.1,0.1))}

     Thus if \code{pci} is supplied, you cannot have a
     previously-existing variable called \code{timeperiod} as a
     covariate in any part of a \code{msm} model.

     Internally, this works by inserting censored observations in the
     data at times when the intensity changes but the state is not
     observed.

     If the supplied times are outside the range of the time variable in
     the data, \code{pci} is ignored and a time-homogeneous model is
     fitted.

     After fitting a time-inhomogeneous model, \code{\link{qmatrix.msm}}
     can be used to obtain the fitted intensity matrices for each time
     period, for example,

     \code{qmatrix.msm(example.msm, covariates=list(timeperiod="[5,Inf)"))}

     This facility does not support interactions between time and
     other covariates.  Such models need to be specified "by hand", using
     a state variable with censored observations inserted.
     Note that the \code{data} component of the \code{msm} object
     returned from a call to \code{msm} with \code{pci} supplied contains the
     states with inserted censored observations and time period
     indicators.  These can be used to construct such models.

     Note that you do not need to use \code{pci} in order to model the
     effect of a time-dependent covariate in the data.  \code{msm} will
     automatically assume that covariates are piecewise-constant and
     change at the times when they are observed.  \code{pci} is for when
     you want all intensities to change at the same pre-specified times
     for all subjects.
   }

   \item{phase.states}{
     Indices of states which have a two-phase sojourn distribution.
     This defines a semi-Markov model, in which the hazard of an onward
     transition depends on the time spent in the state.   

     This uses the technique described by Titman and Sharples (2009).  A
     hidden Markov model is automatically constructed on an expanded
     state space, where the phases correspond to the hidden states.  The
     "tau" proportionality constraint described in this paper is
     currently not supported.

     Covariates, constraints, \code{deathexact} and \code{censor} are
     expressed with respect to the expanded state space.  If not
     supplied by hand, \code{initprobs} is defined automatically so that
     subjects are assumed to begin in the first of the two phases.

     Hidden Markov models can additionally be given phased states.  The
     user supplies an outcome distribution for each original state using
     \code{hmodel}, which is expanded internally so that it is assumed
     to be the same within each of the phased states.  \code{initprobs}
     is interpreted on the expanded state space.  Misclassification
     models defined using \code{ematrix} are not supported, and these
     must be defined using \code{hmmCat} or \code{hmmIdent} constructors,
     as described in the \code{hmodel} section of this help page.  Or
     the HMM on the expanded state space can be defined by hand. 

     Output functions are presented as it were a hidden Markov model on
     the expanded state space, for example, transition probabilities
     between states, covariate effects on transition rates, or
     prevalence counts, are not aggregated over the hidden phases.

     Numerical estimation will be unstable when there is weak evidence
     for a two-phase sojourn distribution, that is, if the model is
     close to Markov. 
     
     See \code{\link{d2phase}} for the definition of the two-phase
     distribution and the interpretation of its parameters.

     This is an experimental feature, and some functions are not
     implemented.  Please report any experiences of using this feature
     to the author!    
     }
   
  \item{phase.inits}{ Initial values for phase-type models.  A list with
    one component for each "two-phased" state.  Each component is itself
    a list of two elements.  The first of these elements is a scalar
    defining the transition intensity from phase 1 to phase 2.  The
    second element is a matrix, with one row for each potential
    destination state from the two-phased state, and two columns.  The
    first column is the transition rate from phase 1 to the destination
    state, and the second column is the transition rate from phase 2 to
    the destination state.  If there is only one destination state, then
    this may be supplied as a vector.

    In phase type models, the initial values for transition rates out of
    non-phased states are taken from the \code{qmatrix} supplied to msm,
    and entries of this matrix corresponding to transitions out of
    phased states are ignored. 
  }

  \item{exacttimes}{By default, the transitions of the Markov process
    are assumed to take place at unknown occasions in between the
    observation times.  If \code{exacttimes} is set to \code{TRUE}, then
    the observation times are assumed to represent the exact times of
    transition of the process.    The subject is assumed
    to be in the same state between these times.  An observation may
    represent a transition to a different state or a repeated observation of
    the same state (e.g. at the end of follow-up).  This is equivalent to
    every row of the data having \code{obstype = 2}.  See the
    \code{obstype} argument.  If both \code{obstype} and
    \code{exacttimes} are specified then \code{exacttimes} is ignored.

    Note that the complete history of the multi-state process is known
    with this type of data.  The models which \pkg{msm} fits have the
    strong assumption of constant (or piecewise-constant) transition
    rates.  Knowing the exact transition times allows
    more realistic models to be fitted with other
    packages.  For example parametric models with sojourn distributions
    more flexible than the exponential can be fitted with the \pkg{flexsurv} package,
    or semi-parametric models can be implemented with \pkg{survival} in
    conjunction with \pkg{mstate}.
  }

  \item{cl}{Width of symmetric confidence intervals for maximum
    likelihood estimates, by default 0.95.}

  \item{fixedpars}{Vector of indices of parameters whose values will be
    fixed at their initial values during the optimisation. These are
    given in the order: transition intensities (reading across rows of
    the transition matrix), covariates on intensities (ordered by
    intensities within covariates), hidden Markov model parameters,
    including misclassification probabilities (ordered by parameters
    within states), hidden Markov model covariate parameters (ordered by
    covariates within parameters within states), initial state occupancy
    probabilities (excluding the first probability, which is fixed at one
    minus the sum of the others).

    If there are equality constraints on certain parameters, then
    \code{fixedpars} indexes the set of unique parameters, excluding
    those which are constrained to be equal to previous parameters.

    To fix all parameters, specify \code{fixedpars = TRUE}.

    This can be useful for profiling likelihoods, and building complex
    models stage by stage. }

  \item{center}{If \code{TRUE} (the default, unless
    \code{fixedpars=TRUE}) then covariates are centered at their means
    during the maximum likelihood estimation. This usually improves
    stability of the numerical optimisation. }

  \item{opt.method}{If "optim", "nlm" or "bobyqa", then the corresponding R
    function will be used for maximum likelihood estimation.
    \code{\link{optim}} is the default.  "bobyqa" requires the package
    \pkg{minqa} to be installed.  See the help of these functions for
    further details.   Advanced users can also add their own
    optimisation methods, see the source for \code{optim.R} in msm for
    some examples.

    If "fisher", then a specialised Fisher scoring method is used
    (Kalbfleisch and Lawless, 1985) which can be faster
    than the generic methods, though less robust.  This is only
    available for Markov models with panel data (\code{obstype=1}), that is, not
    for models with censored states, hidden Markov models, exact
    observation or exact death times (\code{obstype=2,3}).}

  \item{hessian}{If \code{TRUE} then standard errors and confidence
    intervals are obtained from a numerical estimate of the Hessian (the
    observed information matrix).  This is the default when maximum
    likelihood estimation is performed.  If all parameters are fixed at
    their initial values and no optimisation is performed, then this
    defaults to \code{FALSE}.  If requested, the actual Hessian is
    returned in \code{x$paramdata$opt$hessian}, where \code{x} is the
    fitted model object.

    If \code{hessian} is set to \code{FALSE}, then standard
    errors and confidence intervals are obtained from the Fisher
    (expected) information matrix, if this is available.  This may be
    preferable if the numerical estimation of the Hessian is
    computationally intensive, or if the resulting estimate is
    non-invertible or not positive definite.}

  \item{use.deriv}{If \code{TRUE} then analytic first derivatives are
    used in the optimisation of the likelihood, where available and an
    appropriate quasi-Newton optimisation method, such as BFGS, is being
    used.  Analytic derivatives are not available for all models. }

  \item{use.expm}{If \code{TRUE} then any matrix exponentiation needed to
    calculate the likelihood is done using the \pkg{expm} package.
    Otherwise the original routines used in \pkg{msm} 1.2.4 and earlier
    are used.  Set to \code{FALSE} for backward compatibility, and let
    the package maintainer know if this gives any substantive differences.}

  \item{analyticp}{By default, the likelihood for certain simpler 3, 4
    and 5 state models is calculated using an analytic expression for
    the transition probability (P) matrix. For all other models, matrix
    exponentiation is used to obtain P. To revert to the original method
    of using the matrix exponential for all models, specify \code{analyticp=FALSE}. See
    the PDF manual for a list of the models for which analytic P
    matrices are implemented. }

  \item{na.action}{What to do with missing data: either \code{na.omit}
    to drop it and carry on, or \code{na.fail} to stop with an error.
    Missing data includes all NAs in the states, times, \code{subject} or
    \code{obstrue}, all NAs at the first observation for a subject for
    covariates in \code{initcovariates},  all NAs
    in other covariates (excluding the last observation for a subject),
    all NAs in \code{obstype} (excluding the first observation for a
    subject), and any subjects with only one observation (thus no
    observed transitions).
  }

  \item{...}{Optional arguments to the general-purpose \R
    optimisation routine, \code{\link{optim}} by default. For example
    \code{method="Nelder-Mead"} to change the optimisation 
    algorithm from the \code{"BFGS"} method that msm calls by default.

    It is often worthwhile to normalize the optimisation using
    \code{control=list(fnscale = a)}, where \code{a} is the a number of
    the order of magnitude of the -2 log likelihood.
    
    If 'false' convergence is reported and the standard
    errors cannot be calculated due to a non-positive-definite Hessian,
    then consider tightening the tolerance criteria for convergence. If
    the optimisation takes a long time, intermediate steps can be
    printed using the \code{trace} argument of the control list. See
    \code{\link{optim}} for details.

    For the Fisher scoring method, a \code{control} list can be supplied
    in the same way, but the only supported options are \code{reltol},
    \code{trace} and \code{damp}.  The first two are used in the same
    way as for \code{\link{optim}}.  If the algorithm fails with a
    singular information matrix, adjust \code{damp} from the
    default of zero (to, e.g. 1).  This adds a constant identity matrix
    multiplied by \code{damp} to the information matrix during optimisation.
  }
}
\value{
  To obtain summary information from models fitted by the
  \code{\link{msm}} function, it is recommended to use extractor
  functions such as \code{\link{qmatrix.msm}},
  \code{\link{pmatrix.msm}}, \code{\link{sojourn.msm}}, \code{\link{msm.form.qoutput}}.  These provide
  estimates and confidence intervals for quantities such as transition
  probabilities for given covariate values.

  For advanced use, it may be necessary to directly use information
  stored in the object returned by \code{\link{msm}}.  This is
  documented in the help page \code{\link{msm.object}}.

  Printing a \code{msm} object by typing the object's name at the
  command line implicitly invokes \code{\link{print.msm}}.  This
  formats and prints the important information in the model fit,
  and also returns that information in an R object.
  This includes estimates and confidence intervals for the transition intensities and (log) hazard
  ratios for the corresponding covariates. When there is a hidden Markov
  model, the chief information in the
  \code{hmodel} component is also formatted and printed. This includes
  estimates and confidence intervals for each
  parameter.
}
\details{
  For full details about the methodology behind the \pkg{msm} package,
  refer to the PDF manual \file{msm-manual.pdf} in the \file{doc}
  subdirectory of the package. This includes a tutorial in the typical
  use of \pkg{msm}.  The paper by Jackson (2011) in Journal of
  Statistical Software presents the material in this manual in a more
  concise form.

  \pkg{msm} was designed for fitting \emph{continuous-time} Markov
  models, processes where transitions can occur at any time.
  These models are defined by \emph{intensities}, which govern both the
  time spent in the current state and the probabilities of the next
  state.  In \emph{discrete-time models}, transitions are known in
  advance to only occur at multiples of some time unit, and the model is
  purely governed by the probability distributions of the state at the
  next time point, conditionally on the state at the current time. These
  can also be fitted in \pkg{msm}, assuming that there is a
  continuous-time process underlying the data.  Then the fitted
  transition probability matrix over one time period, as
  returned by \code{pmatrix.msm(...,t=1)} is equivalent to the
  matrix that governs the discrete-time model. However, these can be
  fitted more efficiently using multinomial logistic regression, for
  example, using \code{multinom} from the R package \pkg{nnet} (Venables
  and Ripley, 2002).

  For simple continuous-time multi-state Markov models,
  the likelihood is calculated in terms of the transition intensity
  matrix \eqn{Q}. When the data consist of observations of the Markov
  process at arbitrary times, the exact transition times are not known.
  Then the likelihood is calculated using the transition probability
  matrix \eqn{P(t) = \exp(tQ)}{P(t) = exp(tQ)}, where \eqn{\exp}{exp} is
  the matrix exponential.  If state \eqn{i} is observed at time
  \eqn{t} and  state \eqn{j} is observed at time \eqn{u}, then the
  contribution to the likelihood from this pair of observations is
  the \eqn{i,j} element of \eqn{P(u - t)}. See, for example, Kalbfleisch
  and Lawless (1985), Kay (1986), or Gentleman \emph{et al.} (1994).

  For hidden Markov models, the likelihood for an individual
  with \eqn{k} observations is calculated directly by summing over the
  unknown state at each time,  producing a product of \eqn{k} matrices. The
  calculation is a generalisation of the method described by Satten and
  Longini (1996), and also by Jackson and Sharples (2002), and Jackson
  \emph{et al.} (2003).

  There must be enough information in the data on each state to estimate
  each transition rate, otherwise the likelihood will be flat and the
  maximum will not be found.  It may be appropriate to reduce the
  number of states in the model, the number of allowed transitions, or
  the number of covariate effects, to ensure convergence.   Hidden
  Markov models, and situations where the value of the process is only
  known at a series of snapshots, are particularly susceptible to
  non-identifiability, especially when combined with a complex
  transition matrix.
  Choosing an appropriate set of initial values for the optimisation can
  also be important.  For flat likelihoods, 'informative' initial values
  will often be required.   See the PDF manual for other tips.

}
\references{
  Jackson, C.H. (2011). Multi-State Models for Panel Data: The msm Package
  for R., Journal of Statistical Software, 38(8), 1-29. URL
  http://www.jstatsoft.org/v38/i08/.

  Kalbfleisch, J., Lawless, J.F., The analysis of panel data under a
  Markov assumption \emph{Journal of the Americal Statistical
  Association} (1985) 80(392): 863--871.

  Kay, R.  A Markov model for analysing cancer markers and disease
  states in survival studies.  \emph{Biometrics} (1986) 42: 855--865.

  Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P.  Multi-state
  Markov models for analysing incomplete disease history data with
  illustrations for HIV disease.  \emph{Statistics in Medicine} (1994) 13(3):
  805--821.

  Satten, G.A. and Longini, I.M.  Markov chains with measurement error:
  estimating the 'true' course of a marker of the progression of human
  immunodeficiency virus disease (with discussion) \emph{Applied
    Statistics} 45(3): 275-309 (1996)

  Jackson, C.H. and Sharples, L.D. Hidden Markov models for the
  onset and progression of bronchiolitis obliterans syndrome in lung
  transplant recipients \emph{Statistics in Medicine}, 21(1): 113--128
  (2002).

  Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W. and
  Couto, E.  Multi-state Markov models for disease progression with
  classification error. \emph{The Statistician}, 52(2): 193--209 (2003)

  Titman, A.C. and Sharples, L.D.  Semi-Markov models with phase-type
  sojourn distributions.  \emph{Biometrics} 66, 742-752 (2009).

  Venables, W.N. and Ripley, B.D. (2002) \emph{Modern Applied Statistics with
    S}, second edition. Springer.
}

\seealso{
  \code{\link{simmulti.msm}}, \code{\link{plot.msm}},
  \code{\link{summary.msm}}, \code{\link{qmatrix.msm}},
  \code{\link{pmatrix.msm}},  \code{\link{sojourn.msm}}.
}
\examples{
### Heart transplant data
### For further details and background to this example, see
### Jackson (2011) or the PDF manual in the doc directory.
print(cav[1:10,])
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
statetable.msm(state, PTNUM, data=cav)
crudeinits.msm(state ~ years, PTNUM, data=cav, qmatrix=twoway4.q)
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                 qmatrix = twoway4.q, deathexact = 4, 
                 control = list ( trace = 2, REPORT = 1 )  )
cav.msm
qmatrix.msm(cav.msm)
pmatrix.msm(cav.msm, t=10)
sojourn.msm(cav.msm)
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
