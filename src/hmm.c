#include "hmm.h"
#include <Rmath.h>
#include "msm.h"

/* Response (emission) functions for hidden Markov models (probability
   of outcome conditionally on the hidden state).  All of the form
   double f(double x, double *pars) */

/* Categorical distribution on the set (1, 2, 3, ..., pars[0]),
   Baseline category is given by pars[1] (not used any more, to delete?)
   Probabilities are defined by pars[2], ... pars[ncats+1]
   NEW IMPLEMENTATION in v1.3
   pars[2],pars[3]... are absolute probs.  covariates applied in R.
   used to be relative with covariates applied in C.
*/

double hmmCat(double x, double *pars)
{
    int cat = fprec(x, 0);
    int ncats = fprec(pars[0], 0);
    if ((cat > ncats) || (cat < 1)) return 0;
    return pars[1 + cat];
}

double hmmIdent(double x, double *pars)
{
    return all_equal(x, pars[0]);
}

double hmmUnif(double x, double *pars)
{
    double lower = pars[0], upper = pars[1];
    return dunif(x, lower, upper, 0);
}

double hmmNorm(double x, double *pars)
{
    double mean = pars[0], sd = pars[1];
    return dnorm(x, mean, sd, 0);
}

double hmmLNorm(double x, double *pars)
{
    double meanlog = pars[0], sdlog = pars[1];
    return dlnorm(x, meanlog, sdlog, 0);
}

double hmmExp(double x, double *pars)
{
    double mean = 1 / pars[0];
    return dexp(x, mean, 0);
}

double hmmGamma(double x, double *pars)
{
    double shape = pars[0], scale = 1 / pars[1];
    return dgamma(x, shape, scale, 0);
}

double hmmWeibull(double x, double *pars)
{
    double shape = pars[0], scale = pars[1];
    return dweibull(x, shape, scale, 0);
}

double hmmPois(double x, double *pars)
{
    double lambda = pars[0];
    return dpois(x, lambda, 0);
}

double hmmBinom(double x, double *pars)
{
    double size = pars[0], prob = pars[1];
    return dbinom(x, size, prob, 0);
}

/* Truncated normal distribution. Infinite bounds are allowed through
   a parameter with a value of "Inf" or "-Inf" passed from R */

double hmmTNorm(double x, double *pars)
{
    double mean = pars[0], sd = pars[1], lower = pars[2], upper = pars[3];
    double denom = pnorm(upper, mean, sd, 1, 0) - pnorm(lower, mean, sd, 1, 0);
    if (x < lower) return 0;
    if (x > upper) return 0;
    return dnorm(x, mean, sd, 0) / denom;
}

/* Satten and Longini's truncated normal distribution with normal measurement error */

/* To parameterise so covariates go on observation: Put in a dummy
   parameter meanerr = 0 for the measurement error model xobs ~ N(xhid
   + meanerr, sderr), then covs go on meanerr.  */

double hmmMETNorm(double x, double *pars)
{
    double mean = pars[0], sd = pars[1], lower = pars[2], upper = pars[3], sderr = pars[4], meanerr = pars[5];
    double sumsq = sd*sd + sderr*sderr;
    double sigtmp = sd*sderr / sqrt(sumsq);
    double mutmp = ((x - meanerr)*sd*sd + mean*sderr*sderr) / sumsq;
    double nc = 1/(pnorm(upper, mean, sd, 1, 0) - pnorm(lower, mean, sd, 1, 0));
    double nctmp = pnorm(upper, mutmp, sigtmp, 1, 0) - pnorm(lower, mutmp, sigtmp, 1, 0);
    return nc * nctmp * dnorm(x, meanerr + mean, sqrt(sumsq), 0);
}

/* Satten and Longini's uniform distribution with normal measurement error */

double hmmMEUnif(double x, double *pars)
{
    double lower = pars[0], upper = pars[1], sderr = pars[2], meanerr = pars[3];
    return ( pnorm(x, meanerr + lower, sderr, 1, 0) - pnorm(x, meanerr + upper, sderr, 1, 0) ) / (upper - lower) ;
}


double hmmNBinom(double x, double *pars)
{
    double size = pars[0], prob = pars[1];
    return dnbinom(x, size, prob, 0);
}

double hmmBeta(double x, double *pars)
{
    double shape1 = pars[0], shape2 = pars[1];
    return dbeta(x, shape1, shape2, 0);
}

double hmmT(double x, double *pars)
{
    double tmean = pars[0], tscale = pars[1], tdf=pars[2];
    return (1/tscale)*dt((x-tmean)/tscale, tdf, 0);
}
