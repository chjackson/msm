#include "hmm.h"
#include <Rmath.h>

/* Derivatives of PDF w.r.t parameters, for each HMM outcome
 * distribution */

void DhmmCat(double x, double *pars, double *d)
{
/* f(x) = p_y if x=y, and 0 otherwise
   df / dp_y |x   = 1 if x=y, and 0 otherwise */
    int i;
    int cat = fprec(x, 0);
    int ncats = fprec(pars[0], 0);
    for (i=0; i<ncats+1; ++i)
	d[i] = 0;
    if ((cat > ncats) || (cat < 1)) return;
    d[1+cat] = 1;
}

void DhmmIdent(double x, double *pars, double *d){
    d[0] = 0;
}

void DhmmUnif(double x, double *pars, double *d)
{
    d[0]=0; d[1]=0; /* uniform parameters are not estimated */
}

void DhmmNorm(double x, double *pars, double *d)
{
    double mean = pars[0], sd = pars[1], f;
    f = dnorm(x, mean, sd, 0); 
    d[0] = f*(x - mean)/R_pow_di(sd,2);
    d[1] = f*(R_pow_di((x-mean)/sd, 2) - 1)/sd;
}

void DhmmLNorm(double x, double *pars, double *d)
{
    double meanlog = pars[0], sdlog = pars[1], f;
    f = dlnorm(x, meanlog, sdlog, 0);
    d[0] = f*(log(x) - meanlog)/R_pow_di(sdlog,2);
    d[1] = f*(R_pow_di((log(x)-meanlog)/sdlog, 2) - 1)/sdlog;
}

void DhmmExp(double x, double *pars, double *d)
{
    double rate = pars[0];
    d[0] = (1 - rate*x)*exp(-rate*x);
}

void DhmmGamma(double x, double *pars, double *d)
{
    double shape = pars[0], rate=pars[1], scale = 1 / rate, f;
    f = dgamma(x, shape, scale, 0);
    d[0] = f*(log(rate) + log(x) - digamma(shape));
    d[1] = f*(shape/rate - x);
}

void DhmmWeibull(double x, double *pars, double *d)
{
    double shape = pars[0], scale = pars[1], f, rp;
    f = dweibull(x, shape, scale, 0);
    rp = R_pow(x/scale, shape);
    d[0] = f*(1/shape + log(x/scale)*(1 - rp));
    d[1] = f*(shape/scale*(rp - 1));
}

void DhmmPois(double x, double *pars, double *d)
{
    double lambda = pars[0], f;
    f = dpois(x, lambda, 0);
    d[0] = (x/lambda - 1)*f;
}

void DhmmBinom(double x, double *pars, double *d)
{
    double size = pars[0], prob = pars[1], f;
    f = dbinom(x, size, prob, 0);
    d[0] = 0; // fixed
    d[1] = f*(x/prob - (size-x)/(1-prob));
//    printf("f=%.3f, dlf=%.3f, d[1]=%.5f\n", f, (x/prob - (size-x)/(1-prob)), d[1]);
}
 
/* not sure these three are tractable. don't support */
void DhmmTNorm(double x, double *pars, double *d){}
void DhmmMETNorm(double x, double *pars, double *d){}
void DhmmMEUnif(double x, double *pars, double *d){}

void DhmmNBinom(double x, double *pars, double *d)
{
    double size = pars[0], prob = pars[1], f;
    f = dnbinom(x, size, prob, 0);
    d[0] = f*(digamma(x+size) - digamma(size) + log(prob));
    d[1] = f*(size/prob - x/(1-prob));
}

void DhmmBeta(double x, double *pars, double *d)
{
    double shape1 = pars[0], shape2 = pars[1], f;
    f = dbeta(x, shape1, shape2, 0);
    d[0] = f*(digamma(shape1+shape2) - digamma(shape1) + log(x));
    d[1] = f*(digamma(shape1+shape2) - digamma(shape2) + log(1-x));
}

void DhmmT(double x, double *pars, double *d)
{
    double tmean = pars[0], tscale = pars[1], tdf=pars[2], f, xmsq;
    f = (1/tscale)*dt((x-tmean)/tscale, tdf, 0);
    xmsq = (x - tmean)*(x-tmean);
    d[0] = f * (x - tmean)*(tdf + 1) / (tdf*tscale*tscale + xmsq);
    d[1] = f * (-1/tscale + (tdf+1)*xmsq / (tdf*R_pow(tscale,3) + tscale*xmsq));
    d[2] = 0.5 * f * (digamma((tdf + 1)/2) - digamma(tdf/2) - 1/tdf - 
		      log(1 + xmsq / (tdf*tscale*tscale)) + (tdf+1)*xmsq / (R_pow(tdf*tscale,2) + tdf*xmsq));
}
