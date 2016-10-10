/* *****************************************************************
   PROGRAM: lik.c
   AUTHOR:  Chris Jackson
   DATE:    July 2004

   Routines for calculating likelihoods for multi-state Markov and
   hidden Markov models.

   ******************************************************************  */

#include "msm.h"
#include "hmm.h"
#include <Rmath.h>
#include <Rdefines.h>
#define NODEBUG
#define NOVITDEBUG0
#define NOVITDEBUG
#define NODERIVDEBUG

/* MUST KEEP THIS IN SAME ORDER AS .msm.HMODELPARS IN R/constants.R */
hmmfn HMODELS[] = {
    hmmCat,
    hmmIdent,
    hmmUnif,
    hmmNorm,
    hmmLNorm,
    hmmExp,
    hmmGamma,
    hmmWeibull,
    hmmPois,
    hmmBinom,
    hmmTNorm,
    hmmMETNorm,
    hmmMEUnif,
    hmmNBinom,
    hmmBeta,
    hmmT
};

/* MUST KEEP THIS IN SAME ORDER AS .msm.HMODELPARS IN R/constants.R */
dhmmfn DHMODELS[] = {
    DhmmCat,
    DhmmIdent,
    DhmmUnif,
    DhmmNorm,
    DhmmLNorm,
    DhmmExp,
    DhmmGamma,
    DhmmWeibull,
    DhmmPois,
    DhmmBinom,
    DhmmTNorm,
    DhmmMETNorm,
    DhmmMEUnif,
    DhmmNBinom,
    DhmmBeta,
    DhmmT
};

/* MUST MATCH order of .msm.CTASKS in R/constants.R */
#define DO_LIK 0
#define DO_DERIV 1
#define DO_INFO 2
#define DO_VITERBI 3
#define DO_LIK_SUBJ 4
#define DO_DERIV_SUBJ 5
#define DO_DPMAT 6

#define OBS_SNAPSHOT 1
#define OBS_PANEL 1 /* preferred term now */
#define OBS_EXACT 2
#define OBS_DEATH 3

double logit(double x)
{
    return log(x / (1 - x));
}

double expit(double x)
{
    return exp(x) / ( 1 + exp(x) );
}

double identity(double x)
{
    return x;
}

/* Good-enough floating point equality comparison */

int all_equal(double x, double y)
{
    return fabs (x - y) <= DBL_EPSILON * fabs(x);
}

/* For models with censoring: */
/* Return a vector of the nc possible true states that a censored state could represent */
/* These will be summed over when calculating the likelihood */
/* Compare one-indexed obs against one-indexed cm->censor. Return one-indexed current (*states) */

void GetCensored (double obs, cmodel *cm, int *nc, double **states)
{
    int j, k=0, n, cens=0;
    if (cm->ncens == 0)
	n = 1;
    else {
	while ((k < cm->ncens) && !all_equal(obs, cm->censor[k])){ 
	    ++k;
	}
	if (k < cm->ncens) {
	    cens = 1;
	    n =  cm->index[k+1] - cm->index[k];
	}
	else n = 1;
    }
    if (cm->ncens == 0 || !cens)
	(*states)[0] = obs;
    else { for (j = cm->index[k]; j < cm->index[k+1]; ++j)
	(*states)[j - cm->index[k]] = cm->states[j]; }
    *nc = n;
}

/*
  Calculate p (obs curr | true i) for hidden Markov or censoring models
  If observation is not necessarily of the true state (!obstrue),
  then this is just the HMM outcome probability (summed over censor set if necessary)
  If obstrue (observation not misclassified) or censored state.
  e.g. censor set 1,2,3,    state set 1,2,3,4,
  pout =   if i in curr 1, else 0
*/

/* TODO does find_exactdeath_hmm need updating for multivariate observations with different models ? */

/* New obstrue facility 
On entry, obstrue will contain 0 if state unknown, and state if state known
But how do we know if there are any extra outcome data in the outcome variable?
If this is NA, we can ignore it but still use obstrue (TODO put in na.find.msmdata)
If this is a state (eg in misc models) prob of observing it cond on true state is 1. 
If this is a general outcome,  get prob of observing it from HMODELS
*/

void GetOutcomeProb(double *pout, double *outcome, int nc, int nout, double *hpars, hmodel *hm, qmodel *qm, int obstrue)
{
    int i, j, k, ind;
    for (i=0; i<qm->nst; ++i) {
	if (hm->hidden && (obstrue==0)) { /* HMMs with true state not known */
	    if (nout > 1) {  /* multivariate outcomes. Censored states not supported */ 
		pout[i] = 1;
		for (k=0; k<nout; ++k) {
		    ind = (hm->mv ? 
			   MI(k,i,nout) :  /* different models for different variables */
			   i);             /* same model for all */
		    if (!ISNA(outcome[k]) && !ISNA(hm->models[ind])){
			pout[i] *= ((HMODELS[hm->models[ind]])(outcome[k], &(hpars[hm->firstpar[ind]])));
		    }		    
		}
	    } else {  /* Standard univariate HMM (with or without censored state) */
		pout[i] = 0;
		for (j=0; j<nc; ++j){		    
		    pout[i] += (HMODELS[hm->models[i]])(outcome[j], &(hpars[hm->firstpar[i]]));
		}
	    }
	}
	else {      /* True state is known at this time, and appears here as "obstrue" */
	    if (nout > 1){
		pout[i] = 0;
		if (obstrue == i+1){   /* "state" data contain an actual observation. 
					  get its distribution here conditional on the supplied true state */
		    pout[i] = 1;
		    for (k=0; k<nout; ++k) {
			ind = (hm->mv ? MI(k,i,nout) : i);
			if (!ISNA(outcome[k]) && !ISNA(hm->models[ind])){
			    pout[i] *= ((HMODELS[hm->models[ind]])(outcome[k], &(hpars[hm->firstpar[ind]])));
			}
		    }
		}
	    } else {
		pout[i] = 0;
		if (hm->hidden && nc == 1 && !hm->ematrix){ /* "state" data contain an actual observation. 
				 get its distribution here conditional on the supplied true state */
		    if (obstrue == i+1){
			pout[i] = (HMODELS[hm->models[i]])(outcome[0], &(hpars[hm->firstpar[i]]));
		    } 
		} else {  /* "state" data contain the true state (for hm->ematrix with obstrue) or censor indicator */
		    for (j=0; j<nc; ++j){
			if ((int) outcome[j] == i+1)
			    pout[i] = 1;
		    }
		}
	    }
	}
    }
}

/* Get derivative of each state i outcome probability w.r.t. HMM pars
   HMM outcome dist is f(x | p1(a1,b11,b12..), p2(a1,b11,b12),....)
 Need derivs wrt subparameters a1,b11,etc.
 First get derivs wrt p1,p2... from DHMODELS (dptmp).

 Then search for single p1,p2.. that depends on subparameters, and
 multiply by dp1/da, dp1/dpb1 etc, (hm->dpars, calculated in R msm.form.dh)

 We assume there's only one of these (true for standard dists with one
 location parameter, and for categorical outcome where f(x) depends on only
 one pr for a given x.  */

void GetDOutcomeProb(double *dpout, /* qm->nst x hm->nopt */
		     double *outcome, int nc, int nout, double *hpars, hmodel *hm, qmodel *qm, int obsno, int obstrue)
{
    int i, j, k, l, r, s, ind;
    int p=0; /* indexes parameters up to totpars */
    double *pout, *dptmp = Calloc(hm->totpars, double); /* will only use hm->npars[i] slots for each i */
#ifdef DERIVDEBUG
    printf("GetDOutcomeProb:\n");
#endif
    for (i=0; i<qm->nst; ++i) {
	for (l=0; l<hm->nopt; ++l)
	    dpout[MI(i,l,qm->nst)] = 0;
	if (hm->hidden && (!obstrue || (obstrue==(i+1) && !hm->ematrix))) {
	    if (nout > 1) {  /* multivariate outcomes. Censored states not supported			
				TODO.  This is fiddlier than first thought.
				not considered what hm->dpars should be, particularly with constraints  */ 
		pout = Calloc(nout, double);
		for (r=0; r<nout; ++r) {
		    pout[r] = 0; 
		    ind = (hm->mv ? MI(r,i,nout) : i);
		    if (!ISNA(outcome[r]) && !ISNA(hm->models[ind])){
			pout[r] = ((HMODELS[hm->models[ind]])(outcome[r], &(hpars[hm->firstpar[ind]])));
		    }
		}
		for (r=0; r<nout; ++r){
		    ind = (hm->mv ? MI(r,i,nout) : i);
		    if (!ISNA(outcome[r]) && !ISNA(hm->models[ind])){
			(DHMODELS[hm->models[ind]])(outcome[r], &(hpars[hm->firstpar[ind]]), dptmp);
			for (k=0; k<hm->npars[ind]; ++k){
			    for (s=0; s<nout; ++s)
				if (r != s && !ISNA(outcome[s])) dptmp[k] *= pout[s];
			    for (l=0; l<hm->nopt; ++l){
				dpout[MI(i,l,qm->nst)] += dptmp[k] * hm->dpars[MI3(p+k,l,obsno,hm->totpars,hm->nopt)];
#ifdef DERIVDEBUG
				    printf("dpars[%d,%d]=%.2f,", p+k, l, hm->dpars[MI3(p+k,l,obsno,hm->totpars,hm->nopt)]);
#endif
			    }
#ifdef DERIVDEBUG
			    printf("\n");
#endif
#ifdef DERIVDEBUG
			    printf("MI=%d,hm=%d,fp=%d,hp=%f,outcome=%2.0f,dptmp=%f\n",
				   ind,hm->models[ind], hm->firstpar[ind], hpars[hm->firstpar[ind]+1], outcome[r],dptmp[k]);
#endif
			}
#ifdef DERIVDEBUG
			for (l=0; l<hm->nopt; ++l)
			    printf("dpout[%d,%d]=%f,",i,l,dpout[MI(i,l,qm->nst)]);
			printf("\n");
#endif
		    }
		    if (hm->mv) p += hm->npars[ind];
		}
		if (!hm->mv) p += hm->npars[i];
		Free(pout);
	    }
	    else { 
		for (j=0; j<nc; ++j){
		    (DHMODELS[hm->models[i]])(outcome[j], &(hpars[hm->firstpar[i]]), dptmp);
		    for (k=0; k<hm->npars[i]; ++k){
			for (l=0; l<hm->nopt; ++l){
			    dpout[MI(i,l,qm->nst)] += dptmp[k] *
				hm->dpars[MI3(p+k,l,obsno,hm->totpars,hm->nopt)];
			}
		    }
		}
		p += hm->npars[i];
	    }
	}
	else {
	    for (l=0; l<hm->nopt; ++l)
		dpout[MI(i,l,qm->nst)] = 0;	    
	    if (nout > 1 && hm->mv)
		for (r=0; r<nout; ++r)
		    p += hm->npars[MI(r,i,nout)];
	    else p += hm->npars[i];
	}
    }
    Free(dptmp);
}

void normalize(double *in, double *out, int n, double *lweight)
{
    int i; double ave;
    for (i=0, ave=0; i<n; ++i)
	ave += in[i];
    ave /= n;
    if (ave == 0) ave = 1;
    for (i=0; i<n; ++i)
	out[i] = in[i] / ave;
    *lweight -= log(ave);
}

/* Calculate and store all distinct P matrices in the data, given qmat
   stored in regular per-observation format and given "pcomb" (dt,
   covariates, obstype, comb indicator for each obs) */

void calc_p(msmdata *d, qmodel *qm, double *pmat)
{
    double *qmat;
    int pt, i, c, *comb_done = Calloc(d->npcombs, int);
    for (i=0; i<d->npcombs; ++i)
	comb_done[i] = 0;
    for (pt = 0;  pt < d->npts; ++pt){
	for (i = d->firstobs[pt]+1; i <= d->firstobs[pt+1] - 1; ++i) {
	    c = d->pcomb[i];
	    if (!comb_done[c]) {
		qmat = &(qm->intens[MI3(0, 0, i-1, qm->nst, qm->nst)]);
		Pmat(&pmat[MI3(0, 0, c, qm->nst, qm->nst)], d->time[i] - d->time[i-1], qmat, qm->nst,
		     (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm, qm->qperm, qm->expm);
		comb_done[c] = 1;
	    }
	}
    }
    Free(comb_done);
}

void calc_dp(msmdata *d, qmodel *qm, double *dpmat)
{
    double *qmat, *dqmat;
    int pt, i, c, np=qm->nopt;
    //    int r,s,p;
    int *comb_done = Calloc(d->npcombs, int);
    for (i=0; i<d->npcombs; ++i)
	comb_done[i] = 0;
    for (pt = 0;  pt < d->npts; ++pt){
	for (i = d->firstobs[pt]+1; i <= d->firstobs[pt+1] - 1; ++i) {
	    c = d->pcomb[i];
	    if (!comb_done[c]) {
		qmat = &(qm->intens[MI3(0, 0, i-1, qm->nst, qm->nst)]);
		dqmat = &(qm->dintens[MI4(0, 0, 0, i-1, qm->nst, qm->nst, np)]);
		/* dpmat is nst*nst*np*npcombs, leftmost index varies fastest and
		   DPmat returns nst*nst*np */
		DPmat(&dpmat[MI4(0, 0, 0, c, qm->nst, qm->nst, np)], d->time[i] - d->time[i-1], dqmat, qmat, qm->nst, np, (d->obstype[i] == OBS_EXACT));

		//		printf("calc_dp: i=%d,c=%d,r=%d,s=%d,p=%d,dp=%f\n",i,1,0,1,0,dpmat[MI4(0, 1, 0, 1, qm->nst, qm->nst, np)]);
		comb_done[c] = 1;
	    }
	}
    }
    Free(comb_done);
}

/* Find the true state that observation with an exact death time
 * represents in a HMM.  This should be the state with outcome model
 * hmmIdent(obs).  This function also works for non-HMM censoring
 * models, just returning the observed state. */

int find_exactdeath_hmm(double *outcome, int obsno, msmdata *d, qmodel *qm, hmodel *hm){
    int ideath;
    double *hpars = &(hm->pars[MI(0, obsno, hm->totpars)]);
    if (!hm->hidden || d->obstrue[obsno])
	ideath = outcome[0] - 1;
    else
	for (ideath=0; ideath < qm->nst; ++ideath)
	    if (hm->models[ideath] == 1 && hmmIdent(outcome[0], &(hpars[hm->firstpar[ideath]])))
		break;
    return ideath;
}

/* Post-multiply the row-vector cump by matrix T to accumulate the likelihood */

void update_likhidden(double *outcome, int nc, int obsno, msmdata *d, qmodel *qm,
		      hmodel *hm, double *cump, double *newp, double *lweight, Array3 pmat)
{
    int i, j, ideath=0;
    double T, *pout = Calloc(qm->nst, double);
    double *qmat = &(qm->intens[MI3(0, 0, obsno-1, qm->nst, qm->nst)]);
    double *hpars = &(hm->pars[MI(0, obsno, hm->totpars)]);

    GetOutcomeProb(pout, outcome, nc, d->nout, hpars, hm, qm, d->obstrue[obsno]);
    if (d->obstype[obsno] == OBS_DEATH)
	ideath = find_exactdeath_hmm(outcome, obsno, d, qm, hm);

    for(j = 0; j < qm->nst; ++j)
	{
	    newp[j] = 0.0;
	    for(i = 0; i < qm->nst; ++i)
		{
		    if (d->obstype[obsno] == OBS_DEATH)
			T = pmat[MI(i,j,qm->nst)] * qmat[MI(j,ideath,qm->nst)];
		    else {
			T = pmat[MI(i, j, qm->nst)] * pout[j];
//			printf("pmat[%d,%d]=%16.12lf,pout[%d]=%16.12lf\n,", i, j, pmat[MI(i, j, qm->nst)], j, pout[j]);

		    }
		    if (T < 0) T = 0;
		    newp[j] = newp[j] + cump[i]*T;
		}
	}
    /* re-scale the likelihood at each step to prevent it getting too small and underflowing */
    /*  while cumulatively recording the log scale factor   */
    normalize (newp, cump, qm->nst, lweight);
    Free(pout);
}

/* Likelihood for the hidden Markov model for one individual */

double likhidden(int pt, /* ordinal subject ID */
		 msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, Array3 pmat)
{
    double *curr = Calloc (qm->nst, double);
    double *cump     = Calloc(qm->nst, double);
    double *newp     = Calloc(qm->nst, double);
    double *pout = Calloc(qm->nst, double);
    double lweight, lik, *hpars, *outcome;
    int i, obsno, nc=1, allzero=1;
    if (d->firstobs[pt] + 1 == d->firstobs[pt+1])
	return 0; /* individual has only one observation. Shouldn't happen since 1.3.2 */
    /* Likelihood for individual's first observation */
    hpars = &(hm->pars[MI(0, d->firstobs[pt], hm->totpars)]);
    if (d->nout > 1) outcome = &d->obs[MI(0, d->firstobs[pt], d->nout)];
    else {  /* TODO these four lines or similar are pasted a few times */
	GetCensored((double)d->obs[d->firstobs[pt]], cm, &nc, &curr);
	outcome = curr;
    }
    GetOutcomeProb(pout, outcome, nc, d->nout, hpars, hm, qm, d->obstrue[d->firstobs[pt]]);
    /* Likelihood contribution for initial observation */
//    printf("\nlikhidden:\n");
    for (i = 0; i < qm->nst; ++i) {
	cump[i] = pout[i];
//	printf("pout[%d]=%.4f\n",i,pout[i]);
	/* Ignore initprobs if observation is known to be the true state
	   or TODO, can we set it in R to one for obs state, zero for others? */
	if (!d->obstrue[d->firstobs[pt]]) cump[i] = cump[i]*hm->initp[MI(pt,i,d->npts)];
	if (!all_equal(cump[i], 0)) allzero = 0;
    }
    if (allzero && (qm->nliks==1)) {
	warning("First observation of %f for subject number %d out of %d is impossible for given initial state probabilities and outcome model\n", curr[0], pt+1, d->npts);
    }
    lweight=0;
    /* Matrix product loop to accumulate the likelihood for subsequent observations */

    for (obsno = d->firstobs[pt]+1; obsno <= d->firstobs[pt+1] - 1; ++obsno)
    {
	R_CheckUserInterrupt();
	if (d->nout > 1) outcome = &d->obs[MI(0, obsno, d->nout)];
	else {
	    GetCensored((double)d->obs[obsno], cm, &nc, &curr);
	    outcome = curr;
	}
	update_likhidden(outcome, nc, obsno, d, qm, hm, cump, newp, &lweight,
			 &pmat[MI3(0,0,d->pcomb[obsno],qm->nst,qm->nst)]);
    }

    for (i = 0, lik = 0; i < qm->nst; ++i) {
	lik = lik + cump[i];
    }
    Free(curr); Free(cump);  Free(newp); Free(pout);
    /* Transform the likelihood back to the proper scale */
    return -2*(log(lik) - lweight);
}

/* DERIVATIVES FOR HMM using notation from Titman (Lifetime Data Analysis 2009) */

/* Lik and deriv at first obs */
/* Don't support random initprobs for the moment */

void init_hmm_deriv(double *curr, int nc, int pt, int obsno, double *hpars,
		    double *a, double *phi, double *xi, double *dxi,
		    msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm,
		    double *pok, double *dpok)
{
    int i, j, p, n=qm->nst, nqp=qm->nopt, nhp = hm->nopt, np = nqp + nhp;
    double suma, sumphi;
    double *pout = Calloc(n, double);
    double *dpout = Calloc(n * nhp, double);
    int cens_not_hmm = (cm->ncens > 0 && !hm->hidden);
    GetOutcomeProb(pout, curr, nc, d->nout, hpars, hm, qm, d->obstrue[obsno]);
    GetDOutcomeProb(dpout, curr, nc, d->nout, hpars, hm, qm, obsno, d->obstrue[obsno]);
    for (p=0; p<nqp; ++p){
	dpok[p] = 0;
	for (i = 0; i < n; ++i)
	    phi[MI(i,p,n)] = 0;
    }
    suma = 0;
#ifdef DERIVDEBUG
    printf("init_hmm_deriv:\n");
#endif
    for (i = 0; i < n; ++i) {
	/* printf("pout[%d]=%f, ", i, pout[i]); */
	a[i] = (cens_not_hmm ? pout[i] : hm->initp[MI(pt,i,d->npts)] * pout[i]);
#ifdef DERIVDEBUG
	printf("i=%d,initp=%f,pout=%f,a=%f\n", i, hm->initp[MI(pt,i,d->npts)], pout[i], a[i]);
#endif
	suma += a[i];
    }
    /* printf("\n"); */
    *pok = (cens_not_hmm ? 1 : suma);
    for (i = 0; i < n; ++i)
	xi[i] = a[i] / (*pok);
    for (p=0; p<nhp; ++p){
	dpok[nqp+p] = 0;
	for (i = 0; i < n; ++i){
	    phi[MI(i,nqp+p,n)] = (cens_not_hmm ? 0 : hm->initp[MI(pt,i,d->npts)] * dpout[MI(i,p,n)]);
#ifdef DERIVDEBUG
	    printf("p=%d,i=%d,initp=%f,dpout=%f,phi=%f\n", p, i, hm->initp[MI(pt,i,d->npts)], dpout[MI(i,p,n)], phi[MI(i,nqp+p,n)]);
#endif
	    dpok[nqp+p] += phi[MI(i,nqp+p,n)];
	}
    }
    for (p=0; p<np; ++p){
	sumphi = 0;
	for (j=0; j<n; ++j) sumphi += phi[MI(j,p,n)];
	for (i = 0; i < n; ++i){
	    /* Note the denominator should be squared here, error in Titman (2009) */
	    dxi[MI(i,p,n)] = ((*pok)*phi[MI(i,p,n)] - a[i]*sumphi) / ((*pok)*(*pok));
	}
    }
    Free(pout); Free(dpout);
}

void update_hmm_deriv(double *curr, int nc, int obsno,
		      double *pmat, double *dpmat, double *qmat, double *dqmat, double *hpars,
		      double *aold, double *phiold, double *xiold, double *dxiold,
		      double *anew, double *phinew, double *xinew, double *dxinew,
		      msmdata *d, qmodel *qm, hmodel *hm,
		      double *pok, double *dpok)
{
    int i, j, p, s, n=qm->nst, nqp=qm->nopt, nhp = hm->nopt, np = nqp + nhp, ideath=0;
    double qs=0, suma, sumphi, ptrans, dptrans, dqs, dhp;
    double *pout = Calloc(n, double);
    double *dpout = Calloc(n * nhp, double);
    GetOutcomeProb(pout, curr, nc, d->nout, hpars, hm, qm, d->obstrue[obsno]);
    GetDOutcomeProb(dpout, curr, nc, d->nout, hpars, hm, qm, obsno, d->obstrue[obsno]);
    if (d->obstype[obsno] == OBS_DEATH)
	ideath = find_exactdeath_hmm(curr, obsno, d, qm, hm);
#ifdef DERIVDEBUG
    printf("update_hmm_deriv:\n");
#endif
    for (i=0; i<n; ++i){
	anew[i] = 0;
	if (d->obstype[obsno] == OBS_DEATH)
	    qs = qmat[MI(i,ideath,n)];
	for (p=0; p<np; ++p)
	    phinew[MI(i,p,n)] = 0;
	for (j=0; j<n; ++j){
	    ptrans = pmat[MI3(j, i, d->pcomb[obsno], n,n)];
	    if (d->obstype[obsno] == OBS_DEATH)
		anew[i] += aold[j] * ptrans * qs;
	    else
		anew[i] += aold[j] * ptrans * pout[i]; /* typo here in Titman (2009) */
	    for (p=0; p<np; ++p){
		dptrans = (p<nqp ? dpmat[MI4(j, i, p, d->pcomb[obsno], n, n, nqp)] : 0);
		dhp = (p<nqp ? 0 : dpout[MI(i, p-nqp, n)]);
		if (d->obstype[obsno] == OBS_DEATH){
		    dqs = (p<nqp ? dqmat[MI3(i,ideath,p,n,n)] : 0);
		    phinew[MI(i,p,n)] += ptrans * (phiold[MI(j,p,n)]*qs  + aold[j]*dqs) +  aold[j]*qs*dptrans;
		}
		else
		    phinew[MI(i,p,n)] += ptrans * (phiold[MI(j,p,n)]*pout[i] + aold[j]*dhp) + aold[j]*pout[i]*dptrans; /* error in Titman 2009 */
	    }
#ifdef DERIVDEBUG
	    printf("i=%d,j=%d,anew=%f,ideath=%d,ptrans=%f\n",i,j,anew[0],ideath,ptrans);
#endif
	}
    }
    suma = 0;
    for (j=0; j<n; ++j) suma += anew[j];
    for (i=0; i<n; ++i){
	/* def of xi (and dxi) has k off by one in Titman (2009) */
	xinew[i] = anew[i] / suma;
	for (p=0; p<np; ++p){
	    sumphi = 0;
	    for (j=0; j<n; ++j)
		sumphi += phinew[MI(j,p,n)];
	    /* another error here in Titman (2009), denominator not squared */
	    dxinew[MI(i,p,n)] = (suma*phinew[MI(i,p,n)] - anew[i]*sumphi) / (suma*suma);
	}
    }
    *pok = 0;
    for (p=0; p<np; ++p)
	dpok[p] = 0;
    for (j=0; j<n; ++j) {
	for (s=0; s<n; ++s) {
	    if (d->obstype[obsno] == OBS_DEATH){
		qs = qmat[MI(s,ideath,n)];
#ifdef DERIVDEBUG
		printf("s=%d,ideath=%d,MI=%d,qs=%f\n",s,ideath,MI(s,ideath,n),qs);
#endif
	    }
	    ptrans = pmat[MI3(j, s, d->pcomb[obsno], n,n)];
	    if (d->obstype[obsno] == OBS_DEATH)
		*pok += xiold[j] * ptrans * qs;
	    else
		*pok += xiold[j] * ptrans * pout[s];
	    for (p=0; p<np; ++p){
		dptrans = (p<nqp ? dpmat[MI4(j, s, p, d->pcomb[obsno], n, n, nqp)] : 0);
		dhp = (p<nqp ? 0 : dpout[MI(s,p-nqp,n)]);
		if (d->obstype[obsno] == OBS_DEATH){
		    dqs = (p<nqp ? dqmat[MI3(s,ideath,p,n,n)] : 0);
		    dpok[p] += dxiold[MI(j,p,n)] * ptrans * qs  +  xiold[j] * dptrans * qs  + xiold[j] * ptrans * dqs;
		}
		else
		    dpok[p] += dxiold[MI(j,p,n)] * ptrans * pout[s]  +  xiold[j] * dptrans * pout[s]  + xiold[j] * ptrans * dhp;
	    }
#ifdef DERIVDEBUG
	    printf("j=%d,s=%d,qs=%f,pok=%f,xiold=%f,dpok[0]=%f,dpok[1]=%f\n",j,s,qs,*pok,xiold[j],dpok[0],dpok[1]);
#endif
	}
    }
    Free(pout); Free(dpout);
}

void hmm_deriv(int pt, msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, Array3 pmat, Array4 dpmat, double *dlp)
{
    int i, k, p, obsno, nc=1, n=qm->nst, nqp=qm->nopt, nhp = hm->nopt, np = nqp + nhp;
    double lp, pok;
    double *curr = Calloc (n, double);
    int nobspt = d->firstobs[pt+1] - d->firstobs[pt];
    double *anew = Calloc(n, double); /* alpha_k(i) = P(X(t_k) = i, o_1, ..., o_k),    X is true, O is obs */
    double *aold = Calloc(n, double);
    double *phinew = Calloc(n * np, double);   /* phi_k(theta_m, i) = deriv of alpha wrt theta_m. */
    double *phiold = Calloc(n * np, double);
    double *xinew = Calloc(n, double); /* xi(k,i) = P(x(k-1) = i | o_1, ..., o_(k-1)) */
    double *xiold = Calloc(n, double);
    double *dxinew = Calloc(n * np, double);
    double *dxiold = Calloc(n * np, double);
    double *dpok = Calloc(np, double);
    double *qmat, *dqmat, *hpars=NULL, *outcome=NULL;
    if (hm->hidden) hpars = &(hm->pars[MI(0, d->firstobs[pt], hm->totpars)]);
    if (d->nout > 1) outcome = &d->obs[MI(0, d->firstobs[pt], d->nout)];
    else { 
	GetCensored((double)d->obs[d->firstobs[pt]], cm, &nc, &curr);
	outcome = curr;
    }
    // Get lik and deriv at first obs
    init_hmm_deriv(outcome, nc, pt, d->firstobs[pt], hpars,
		   aold, phiold, xiold, dxiold,
		   d, qm, cm, hm, &pok, dpok);
    lp = log(pok);
    /* for (i=0;i<d->nout;++i) printf("outcome[%d]=%f,",i,outcome[i]); printf("\n"); */
#ifdef DERIVDEBUG
    printf("hmm_deriv:\n");
#endif
    for (p=0; p<np; ++p){
	dlp[p] = dpok[p] / pok;
#ifdef DERIVDEBUG
	printf("k=0,dpok[%d]=%f,pok=%f,dlp[%d]=%f\n",p,dpok[p],pok,p,dlp[p]);
#endif
    }
    /* printf("\n"); */
    // Subsequent observations, using forward algorithm
    
    for (k=1; k<nobspt; ++k) {
	obsno = d->firstobs[pt] + k;
	qmat = &(qm->intens[MI3(0, 0, obsno-1, n, n)]);
	dqmat = &(qm->dintens[MI4(0, 0, 0, obsno - 1, n, n, nqp)]);
	hpars = &(hm->pars[MI(0, obsno, hm->totpars)]);
	if (d->nout > 1) outcome = &d->obs[MI(0, obsno, d->nout)];
	else { 
	    GetCensored((double)d->obs[obsno], cm, &nc, &curr);
	    outcome = curr;
	}
	update_hmm_deriv(outcome, nc, obsno,
			 pmat, dpmat, qmat, dqmat, hpars,
			 aold, phiold, xiold, dxiold,
			 anew, phinew, xinew, dxinew,
			 d, qm, hm,
			 &pok, dpok);
	for (i=0; i<n; ++i){
	    // Set a to xi to avoid underflow
	    aold[i] = xinew[i];  // works since pok,dpok only depend on a through xi
	    xiold[i] = xinew[i];
	    for (p=0; p<np; ++p){
		phiold[MI(i,p,n)] = dxinew[MI(i,p,n)];
		dxiold[MI(i,p,n)] = dxinew[MI(i,p,n)];
	    }
	}
	lp += log(pok);
	for (p=0; p<np; ++p){
	    dlp[p] += dpok[p] / pok;
#ifdef DERIVDEBUG
	    printf("k=%d,dlp[%d]=%f,",k,p,dlp[p]); 
#endif
	
#ifdef DERIVDEBUG
	printf("\n"); 
#endif
	}
    }
	
    lp *= -2;
    Free(curr); Free(aold); Free(anew); Free(phiold); Free(phinew); Free(xiold); Free(xinew); Free(dxiold); Free(dxinew); Free(dpok);
}


void hmm_info(int pt, msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, Array3 pmat, Array4 dpmat, double *info)
{
    int i, j, k, p, q, obsno, nc=1, n=qm->nst, nqp=qm->nopt, nhp = hm->nopt, np = nqp + nhp;
    double pok;
    double *curr = Calloc (n, double);
    double *potential = Calloc (n, double);
    int nobspt = d->firstobs[pt+1] - d->firstobs[pt];
    double *anew = Calloc(n, double); /* alpha_k(i) = P(X(t_k) = i, o_1, ..., o_k),    X is true, O is obs */
    double *aold = Calloc(n, double);
    double *phinew = Calloc(n * np, double);   /* phi_k(theta_m, i) = deriv of alpha wrt theta_m. */
    double *phiold = Calloc(n * np, double);
    double *xinew = Calloc(n, double); /* xi(k,i) = P(x(k-1) = i | o_1, ..., o_(k-1)) */
    double *xiold = Calloc(n, double);
    double *dxinew = Calloc(n * np, double);
    double *dxiold = Calloc(n * np, double);
    double *dpok = Calloc(np, double);
    double *qmat, *dqmat, *hpars=NULL, *outcome;
    if (hm->hidden) hpars = &(hm->pars[MI(0, d->firstobs[pt], hm->totpars)]);
    for (p=0; p<np; ++p)
	for (q=0; q<np; ++q)
	    info[MI(q,p,np)] = 0;
    for (j=0; j<n; ++j) {  // calc lik and deriv for potential obs at each time, given previous actual ones
	*potential = j+1; nc=1;
	init_hmm_deriv(potential, nc, pt, d->firstobs[pt], hpars,
		       anew, phinew, xinew, dxinew, // discard, not used for updating
		       d, qm, cm, hm, &pok, dpok);
	for (p=0; p<np; ++p)
	    for (q=0; q<np; ++q){
		if (pok > 0)
		    info[MI(q,p,np)] += dpok[p]*dpok[q] / pok;
//		printf("k=0,j=%d,dpok[%d]=%f,dpok[%d]=%f,pok=%f,info=%f\n",j,p,dpok[p],q,dpok[q],pok,info[MI(q,p,np)]);
	    }
    }
    if (d->nout > 1) outcome = &d->obs[MI(0, d->firstobs[pt], d->nout)];
    else { 
	GetCensored((double)d->obs[d->firstobs[pt]], cm, &nc, &curr);
	outcome = curr;
    }
    init_hmm_deriv(outcome, nc, pt, d->firstobs[pt], hpars,
		   aold, phiold, xiold, dxiold, // use actual observation to update these
		   d, qm, cm, hm, &pok, dpok);
    // Subsequent observations, using forward algorithm
    for (k=1; k<nobspt; ++k) {
	obsno = d->firstobs[pt] + k;
	if (d->obstype[obsno] != OBS_PANEL)
	    error("Fisher information only available for panel data\n");
	qmat = &(qm->intens[MI3(0, 0, obsno-1, n, n)]);
	dqmat = &(qm->dintens[MI4(0, 0, 0, obsno - 1, n, n, nqp)]);
	hpars = &(hm->pars[MI(0, obsno, hm->totpars)]);
	for (j=0; j<n; ++j) {
	    *potential = j+1; nc=1;
	    update_hmm_deriv(potential, nc, obsno, pmat, dpmat, qmat, dqmat, hpars,
			     aold, phiold, xiold, dxiold,
			     anew, phinew, xinew, dxinew, // these are discarded
			     d, qm, hm,  &pok, dpok);
	    for (p=0; p<np; ++p)
		for (q=0; q<np; ++q){
		    if (pok > 0)
			info[MI(q,p,np)] += dpok[p]*dpok[q] / pok;
//		    printf("k=%d,j=%d,dpok[%d]=%f,dpok[%d]=%f,pok=%f,info=%f\n",k,j,p,dpok[p],q,dpok[q],pok,info[MI(q,p,np)]);
		}
	}
	if (d->nout > 1) outcome = &d->obs[MI(0, obsno, d->nout)];
	else {
	    GetCensored((double)d->obs[obsno], cm, &nc, &curr); // update using observed data
	    outcome = curr;
	}
	update_hmm_deriv(outcome, nc, obsno, pmat, dpmat, qmat, dqmat, hpars,
			 aold, phiold, xiold, dxiold,
			 anew, phinew, xinew, dxinew,
			 d, qm, hm, &pok, dpok);
	for (i=0; i<n; ++i){
	    aold[i] = xinew[i];
	    xiold[i] = xinew[i];
	    for (p=0; p<np; ++p){
		phiold[MI(i,p,n)] = dxinew[MI(i,p,n)];
		dxiold[MI(i,p,n)] = dxinew[MI(i,p,n)];
	    }
	}
    }
    Free(curr); Free(potential); Free(anew); Free(aold); Free(phiold); Free(phinew); Free(xinew); Free(xiold); Free(dxiold); Free(dxinew); Free(dpok);
}



void update_likcensor(int obsno, double *prev, double *curr, int np, int nc,
			msmdata *d, qmodel *qm,  hmodel *hm,
		      double *cump, double *newp, double *lweight, Array3 pmat)
{
    double *qmat = &(qm->intens[MI3(0, 0, obsno-1, qm->nst, qm->nst)]);
    double contrib;
    int i, j, k;

    for(i = 0; i < nc; ++i)
	{
	    newp[i] = 0.0;
	    for(j = 0; j < np; ++j) {
		if (d->obstype[obsno] == OBS_DEATH) {
		    contrib = 0;
		    for (k = 0; k < qm->nst; ++k)
			if (k != curr[i]-1)
			    contrib += pmat[MI((int) prev[j]-1, k, qm->nst)] *
				qmat[MI(k, (int) curr[i]-1, qm->nst)];
		    newp[i] += cump[j] * contrib;

		}
		else {
		    newp[i] += cump[j] * pmat[MI((int) prev[j]-1, (int) curr[i]-1, qm->nst)];
		}
	    }
	}
    normalize(newp, cump, nc, lweight);
}

double likcensor(int pt, /* ordinal subject ID */
		 msmdata *d, qmodel *qm,
		 cmodel *cm, hmodel *hm, Array3 pmat
		 )
{
    double *cump     = Calloc(qm->nst, double);
    double *newp     = Calloc(qm->nst, double);
    double *prev     = Calloc(qm->nst, double);
    double *curr     = Calloc(qm->nst, double);
    double lweight = 0, lik;
    int i, obs, np=0, nc=0;
    if (d->firstobs[pt] + 1 == d->firstobs[pt+1])
      return 0; /* individual has only one observation */
    for (i = 0; i < qm->nst; ++i)
	cump[i] = 1;
    GetCensored((double)d->obs[d->firstobs[pt]], cm, &np, &prev);
    for (obs = d->firstobs[pt]+1; obs <= d->firstobs[pt+1] - 1; ++obs)
	{
	    /* post-multiply by sub-matrix of P at each obs */
	    GetCensored((double)d->obs[obs], cm, &nc, &curr);
	    update_likcensor(obs, prev, curr, np, nc, d, qm, hm,
			     cump, newp, &lweight,
			     &pmat[MI3(0,0,d->pcomb[obs],qm->nst,qm->nst)]);
	    np = nc;
	    for (i=0; i<nc; ++i) prev[i] = curr[i];
	}
    for (i = 0, lik = 0; i < nc; ++i)
	lik = lik + cump[i];
    Free(cump);  Free(newp);  Free(prev); Free(curr);
    return -2*(log(lik) - lweight);
}

/* Likelihood for the non-hidden multi-state Markov model. Data of
   form "time-difference, covariates, from-state, to-state, number of
   occurrences" */

double liksimple(msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm)
{
    int i;
    double lik=0, contrib=0;
    double *pmat = Calloc((qm->nst)*(qm->nst), double), *qmat;
    qmat = &(qm->intens[MI3(0, 0, 0, qm->nst, qm->nst)]);
    for (i=0; i < d->nagg; ++i)
	{
	    R_CheckUserInterrupt();
	    if ((i==0) || (d->whicha[i] != d->whicha[i-1]) || (d->obstypea[i] != d->obstypea[i-1])) {
		/* we have a new timelag/covariates/obstype combination. Recalculate the
		   P matrix for this */
		/* pointer to Q matrix for ith datapoint */
		qmat = &(qm->intens[MI3(0, 0, i, qm->nst, qm->nst)]);
		Pmat(pmat, d->timelag[i], qmat, qm->nst, (d->obstypea[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
	    }
	    if (d->obstypea[i] == OBS_DEATH) {
		contrib = pijdeath(d->fromstate[i], d->tostate[i], pmat, qmat, qm->nst);
	    }
	    else
		contrib = pmat[MI(d->fromstate[i], d->tostate[i], qm->nst)];
	    lik += d->nocc[i] * log(contrib);
#ifdef DEBUG
/*	    printf("obs %d, from %d, to %d, time %lf, obstypea %d, ", i, d->fromstate[i], d->tostate[i], d->timelag[i], d->obstypea[i]);
	    printf("nocc %d, con %lf, lik %lf\n", d->nocc[i], log(contrib), lik);*/
//	    printf("%d-%d in %lf, q=%lf,%lf, lik=%20.20lf, ll=%lf\n",d->fromstate[i], d->tostate[i], d->timelag[i],qmat[0],qmat[1], contrib, d->nocc[i] * log(contrib));
#endif
	}
    Free(pmat);
    return (-2*lik);
}

/* Likelihood for the non-hidden multi-state Markov model, by subject */

double liksimple_subj(int pt, /* ordinal subject ID */
		      msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm)
{
    int i, from, to;
    double lik=0, pm=0, dt;
    double *pmat = Calloc((qm->nst)*(qm->nst), double), *qmat;
    for (i = d->firstobs[pt]+1; i < d->firstobs[pt+1]; ++i) {
	R_CheckUserInterrupt();
	dt = d->time[i] - d->time[i-1];
	from = fprec(d->obs[i-1] - 1, 0); /* convert state outcome to integer */
	to = fprec(d->obs[i] - 1, 0);
	qmat = &(qm->intens[MI3(0, 0, i-1, qm->nst, qm->nst)]); /* use covariate at start of transition */
	Pmat(pmat, dt, qmat, qm->nst, (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
	if (d->obstype[i] == OBS_DEATH)
	    pm = pijdeath(from, to, pmat, qmat, qm->nst);
	else
	    pm = pmat[MI(from, to, qm->nst)];
#ifdef DEBUG
	printf("i=%d, %d-%d in %lf, q=%lf,%lf, ll=%lf\n",i,from,to,dt,qmat[0],qmat[1],log(pm));
#endif
	lik += log(pm);
    }
    Free(pmat);
    return (-2*lik);
}

void msmLikelihood (msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *returned)
{
    int pt;
    double likone;
    double *pmat = Calloc((qm->nst)*(qm->nst)*(d->npcombs), double);
    *returned = 0;
    /* Likelihood for hidden Markov model */
    if (hm->hidden)
	{
 	    calc_p(d, qm, pmat);
	    for (pt = 0;  pt < d->npts; ++pt){
		likone = likhidden (pt, d, qm, cm, hm, pmat);
#ifdef DEBUG
		printf("pt %d, lik %lf\n", pt, likone);
#endif
		*returned += likone;
	    }
	}
    /* Likelihood for Markov model with censored outcomes */
    else if (cm->ncens > 0)
	{
 	    calc_p(d, qm, pmat);
	    for (pt = 0;  pt < d->npts; ++pt){
		likone = likcensor (pt, d, qm, cm, hm, pmat);
		*returned += likone;
	    }
	}
    /* Likelihood for simple non-hidden, non-censored Markov model */
    else {
	*returned = liksimple (d, qm, cm, hm);
    }
    Free(pmat);
}

/* First derivatives of the log-likelihood for the non-hidden
   multi-state Markov model. */

void derivsimple(msmdata *d, qmodel *qm,  cmodel *cm, hmodel *hm, double *deriv)
{
    int i, p, np=qm->nopt;
    double *qmat, *dqmat;
    double *pmat = Calloc(qm->nst * qm->nst, double);
    double *dpmat = Calloc(qm->nst * qm->nst * np, double);
    double *dp = Calloc(np, double);
    double pm;
    qmat = &(qm->intens[MI3(0, 0, 0, qm->nst, qm->nst)]);
    dqmat = &(qm->dintens[MI4(0, 0, 0, 0, qm->nst, qm->nst, np)]);
    for (p = 0; p < np; ++p) deriv[p] = 0;
    for (i=0; i < d->nagg; ++i)
    {
	    R_CheckUserInterrupt();
	    if ((i==0) || (d->whicha[i] != d->whicha[i-1]) || (d->obstypea[i] != d->obstypea[i-1])) {
		/* we have a new timelag/covariates/obstype combination. Recalculate the
		   P matrix and its derivatives for this */
		qmat = &(qm->intens[MI3(0, 0, i, qm->nst, qm->nst)]);
		Pmat(pmat, d->timelag[i], qmat, qm->nst, (d->obstypea[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
		dqmat = &(qm->dintens[MI4(0, 0, 0, i, qm->nst, qm->nst, np)]);
 		DPmat(dpmat, d->timelag[i], dqmat, qmat, qm->nst, np, (d->obstypea[i] == OBS_EXACT));
	    }
	    if (d->obstypea[i] == OBS_DEATH) {
		pm = pijdeath(d->fromstate[i], d->tostate[i], pmat, qmat, qm->nst);
		dpijdeath(d->fromstate[i], d->tostate[i], dpmat, pmat, dqmat, qmat, qm->nst, np, dp);
	    }
	    else {
		pm = pmat[MI(d->fromstate[i], d->tostate[i], qm->nst)];
		for (p = 0; p < np; ++p)
		    dp[p] = dpmat[MI3(d->fromstate[i], d->tostate[i], p, qm->nst, qm->nst)];
	    }
	    for (p = 0; p < np; ++p) {
		if (pm > 0)
		    deriv[p] += d->nocc[i] * dp[p] / pm;
	    }
    }
    for (p = 0; p < np; ++p) {
	deriv[p] *= -2;   /* above is dlogL/dtu as in Kalb+Law, we want
			     deriv of -2 loglik  */
    }
    Free(pmat); Free(dpmat); Free(dp);
}

/* First derivatives of the likelihood for the non-hidden multi-state
   Markov model. Uses data of form "subject, time, state, covariates".
   Returns derivatives by individual and parameter for use in the
   score residual diagnostic.
*/

void derivsimple_subj(msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *deriv)
{
    int pt, i, p, np=qm->nopt;
    double *qmat, *dqmat;
    double *pmat = Calloc(qm->nst * qm->nst, double);
    double *dpmat = Calloc(qm->nst * qm->nst * np, double);
    double *dp = Calloc(np, double);
    double pm=0, dt;
    int from, to;
    for (pt = 0;  pt < d->npts; ++pt){
	{
	    R_CheckUserInterrupt();
	    if (d->firstobs[pt+1] > d->firstobs[pt] + 1) { /* individual has more than one observation? */
	      for (p = 0; p < np; ++p) {
		deriv[MI(pt,p,d->npts)] = 0;
	      }
	      for (i = d->firstobs[pt]+1; i < d->firstobs[pt+1]; ++i) {
		    dt = d->time[i] - d->time[i-1];
		    from = fprec(d->obs[i-1] - 1, 0); /* convert state outcome to integer */
		    to = fprec(d->obs[i] - 1, 0);
		    qmat = &(qm->intens[MI3(0, 0, i-1, qm->nst, qm->nst)]); /* use covariate at start of transition */
		    Pmat(pmat, dt, qmat, qm->nst, (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
		    dqmat = &(qm->dintens[MI4(0, 0, 0, i-1, qm->nst, qm->nst, np)]);
		    DPmat(dpmat, dt, dqmat, qmat, qm->nst, np, (d->obstype[i] == OBS_EXACT));
		    if (d->obstype[i] == OBS_DEATH) {
			pm = pijdeath(from, to, pmat, qmat, qm->nst);
			dpijdeath(from, to, dpmat, pmat, dqmat, qmat, qm->nst, np, dp);
		    }
		    else {
			pm = pmat[MI(from, to, qm->nst)];
			for (p = 0; p < np; ++p)
			    dp[p] = dpmat[MI3(from, to, p, qm->nst, qm->nst)];
		    }
		    for (p = 0; p < np; ++p) {
		      deriv[MI(pt,p,d->npts)] += dp[p] / pm; /* on loglik scale not -2*loglik */
#ifdef DEBUG
		      printf("i=%d, p=%d, pt=%d, from=%d, to=%d, dt=%6.4f, %f, %f, %lf, %lf\n", i, p, pt, from, to, dt, dp[p], pm, -2 * dp[p] / pm, -2*deriv[MI(pt,p,d->npts)]);
#endif
		    }
		}
		for (p = 0; p < np; ++p)
		    deriv[MI(pt,p,d->npts)] *= -2;
	    }
	    else for (p = 0; p < np; ++p)
		deriv[MI(pt,p,d->npts)] = 0;
	}
    }
    Free(pmat); Free(dpmat); Free(dp);
}

void infosimple(msmdata *d, qmodel *qm,  cmodel *cm, hmodel *hm, double *info)
{
    int i, j, p, q, np=qm->nopt;
    double *qmat, *dqmat;
    double *pmat = Calloc(qm->nst * qm->nst, double);
    double *dpmat = Calloc(qm->nst * qm->nst * np, double);
    double *dpm = Calloc(qm->nst* np, double);
    double *pm = Calloc(qm->nst, double);
    for (p = 0; p < np; ++p)
	for (q = 0; q < np; ++q)
	    info[MI(p,q,np)] = 0;
    for (i=0; i < d->nagg; ++i)
	{
	    R_CheckUserInterrupt();
	    if ((i==0) || (d->whicha[i] != d->whicha[i-1]) || (d->obstypea[i] != d->obstypea[i-1])) {
		/* we have a new timelag/covariates/obstype combination. Recalculate the
		   P matrix and its derivatives for this */
		qmat = &(qm->intens[MI3(0, 0, i, qm->nst, qm->nst)]);
		Pmat(pmat, d->timelag[i], qmat, qm->nst, (d->obstypea[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
		dqmat = &(qm->dintens[MI4(0, 0, 0, i, qm->nst, qm->nst, np)]);
 		DPmat(dpmat, d->timelag[i], dqmat, qmat, qm->nst, np, (d->obstypea[i] == OBS_EXACT));
	    }
	    if (d->obstypea[i] != OBS_PANEL)
		error("Fisher information only available for panel data\n");
	    for (j=0; j < qm->nst; ++j)  {
		pm[j] = pmat[MI(d->fromstate[i], j, qm->nst)];
		for (p = 0; p < np; ++p){
		    dpm[MI(j,p,qm->nst)] = dpmat[MI3(d->fromstate[i], j, p, qm->nst, qm->nst)];
		}
	    }
	    if ((i==0) || (d->whicha[i] != d->whicha[i-1]) || (d->obstypea[i] != d->obstypea[i-1]) ||
		(d->fromstate[i] != d->fromstate[i-1])) {
		for (p = 0; p < np; ++p) {
		    for (q = 0; q < np; ++q) {
			/* for expected information, sum over all possible destination states for this fromstate/time/cov combination */
			for(j = 0; j<qm->nst; ++j)  {
			    if (pm[j] > 0)
				info[MI(p,q,np)] +=  d->noccsum[i] * dpm[MI(j,p,qm->nst)] * dpm[MI(j,q,qm->nst)] /  pm[j];
			}
		    }
		}
	    }
	}
    for (p = 0; p < np; ++p) {
	for (q = 0; q < np; ++q) {
	    info[MI(p,q,np)] *= 2; /* above is E(-d2logL/dtudtv) as in Kalb+Law.
					   we want second deriv of of -2 loglik */
	}
    }
    Free(pm); Free(dpm); Free(dpmat); Free(pmat);
}

/* Derivatives of the P matrix, used for the Pearson test p-value */
/* Returns a ntrans * ntostates * npars matrix */
/* Panel data only (obstype 1) */

void dpmat_obs(msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *deriv)
{
    int pt, i, j, k, p, from, np = qm->nopt;
    double *dpmat = Calloc(qm->nst * qm->nst * np, double), *qmat, *dqmat;
    double dt;

    j=0;
    for (pt = 0;  pt < d->npts; ++pt)
	{
	    R_CheckUserInterrupt();
	    if (d->firstobs[pt+1] > d->firstobs[pt] + 1) { /* individual has more than one observation? */
		for (i = d->firstobs[pt]+1; i < d->firstobs[pt+1]; ++i) {
		    ++j;
		    qmat = &(qm->intens[MI3(0, 0, i, qm->nst, qm->nst)]);
		    dqmat = &(qm->dintens[MI4(0, 0, 0, i, qm->nst, qm->nst, np)]);
		    dt = d->time[i] - d->time[i-1];
		    from = fprec(d->obs[i-1] - 1, 0); /* convert state outcome to integer */
		    DPmat(dpmat, dt, dqmat, qmat, qm->nst, np, (d->obstype[i] == OBS_EXACT));
		    for (p = 0; p < np; ++p) {
			for (k=0; k < qm->nst; ++k) {
			    deriv[MI3(j-1,k,p,d->ntrans,qm->nst)] = dpmat[MI3(from, k, p, qm->nst, qm->nst)];
			}
		    }
		}
	    }
	}
    Free(dpmat);
}

void derivhidden(msmdata *d, qmodel *qm,  cmodel *cm, hmodel *hm, double *deriv, int by_subject){
    int pt, p, np = qm->nopt + hm->nopt;
    double *pmat = Calloc((qm->nst)*(qm->nst)*(d->npcombs), double);
    double *dpmat = Calloc((qm->nst)*(qm->nst)*qm->nopt*(d->npcombs), double);
    double *dlp = Calloc(np, double);
    calc_p(d, qm, pmat);
    calc_dp(d, qm, dpmat);
    if (!by_subject)
	for (p=0; p<np; ++p)
	    deriv[p] = 0;
    for (pt = 0;  pt < d->npts; ++pt){
	hmm_deriv(pt, d, qm, cm, hm, pmat, dpmat, dlp);
	for (p=0; p<np; ++p){
	    if (by_subject)
		deriv[MI(pt,p,d->npts)] = -2*dlp[p];
	    else
		deriv[p] += -2*dlp[p];
	}
    }
    Free(pmat); Free(dpmat); Free(dlp);
}

void msmDeriv_subj (msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *returned)
{
    if (hm->hidden || (cm->ncens > 0))
	derivhidden(d, qm, cm, hm, returned, 1);
    else
	derivsimple_subj(d, qm, cm, hm, returned);
}

// Information matrix for HMMs.  Only support panel data and misclassification models

void infohidden(msmdata *d, qmodel *qm,  cmodel *cm, hmodel *hm, double *info){
    int pt, p, q, np = qm->nopt + hm->nopt;
    double *pmat = Calloc((qm->nst)*(qm->nst)*(d->npcombs), double);
    double *dpmat = Calloc((qm->nst)*(qm->nst)*qm->nopt*(d->npcombs), double);
    double *itmp = Calloc(np*np, double);
    calc_p(d, qm, pmat);
    calc_dp(d, qm, dpmat);
    for (p=0; p<np; ++p)
	for (q=0; q<np; ++q)
	    info[MI(q,p,np)] = 0;
    for (pt = 0;  pt < d->npts; ++pt){
	hmm_info(pt, d, qm, cm, hm, pmat, dpmat, itmp);
	for (p=0; p<np; ++p)
	    for (q=0; q<np; ++q){
		info[MI(q,p,np)] += 2 * itmp[MI(q,p,np)];
//		printf("itmp[%d,%d]=%f,info[%d,%d]=%f\n",q,p,itmp[MI(q,p,np)],q,p,info[MI(q,p,np)]);
	    }
    }
    Free(pmat); Free(dpmat); Free(itmp);
}

void msmDeriv (msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *returned)
{
    if (hm->hidden || cm->ncens > 0)
	derivhidden(d, qm, cm, hm, returned, 0);
    else
	derivsimple(d, qm, cm, hm, returned);
}

void msmInfo (msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *returned)
{
    if (hm->hidden || (cm->ncens > 0))
	infohidden(d, qm, cm, hm, returned);
    else
	infosimple(d, qm, cm, hm, returned);
}

/* Return vector of subject-specific log likelihoods */

void msmLikelihood_subj (msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *returned)
{
    int pt;
    double *pmat = Calloc((d->npcombs)*(qm->nst)*(qm->nst), double);
    if (hm->hidden || (cm->ncens > 0))
	calc_p(d, qm, pmat);
    for (pt = 0;  pt < d->npts; ++pt){
	if (hm->hidden)
	    returned[pt] = likhidden (pt, d, qm, cm, hm, pmat);
	else if (cm->ncens > 0)
	    returned[pt] = likcensor (pt, d, qm, cm, hm, pmat);
	else
	    returned[pt] = liksimple_subj (pt, d, qm, cm, hm);
    }
    Free(pmat);
}

/* Find zero-based index of maximum element of a vector x */

void pmax(double *x, int n, int *maxi)
{
    int i=0;
    *maxi = i;
    for (i=1; i<n; ++i) {
	if (x[i] > x[*maxi]) {
	    *maxi = i;
	}
    }
}

/* Calculates the most likely path through underlying states
   and "posterior" probabilities of underlying states 
 */

void Viterbi(msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *fitted, double *pstate)
{
    int i, j, tru, k, kmax, obs, nc = 1;

    double *pmat = Calloc((qm->nst)*(qm->nst), double);
    int *ptr = Calloc((d->n)*(qm->nst), int);
    double *lvold = Calloc(qm->nst, double);
    double *lvnew = Calloc(qm->nst, double);
    double *lvp = Calloc(qm->nst, double);
    double *curr = Calloc (qm->nst, double), *outcome;
    double *pout = Calloc(qm->nst, double);
    double *pfwd = Calloc((d->n)*(qm->nst), double);
    double *pbwd = Calloc((d->n)*(qm->nst), double);

    double dt, logpall, psum;
    double *qmat, *hpars;
    double *ucfwd = Calloc(d->n, double);
    double *ucbwd = Calloc(d->n, double);

    i = 0;
    if (d->obstrue[i]) {
      for (k = 0; k < qm->nst; ++k)
	lvold[k] = (k+1 == d->obstrue[i] ? 0 : R_NegInf);
    }
    else {
      if (d->nout > 1) outcome = &d->obs[MI(0, i, d->nout)];
      else {
        GetCensored(d->obs[i], cm, &nc, &curr);
        outcome = curr; 
      }
      /* initial observation is a censored state. No HMM here, so initprobs not needed */
      if (nc > 1) {
	for (k = 0, j = 0; k < qm->nst; ++k) {
	  if (k+1 == outcome[j]) {
	    lvold[k] = 0;
	    ++j;
	  }
	  else lvold[k] = R_NegInf;
	}
      }
      else { /* use initprobs */
	for (k = 0; k < qm->nst; ++k)
	    lvold[k] = log(hm->initp[MI(0, k, d->npts)]);
      }
    }
    for (k = 0; k < qm->nst; ++k)
    	pfwd[MI(i,k,d->n)] = exp(lvold[k]);
    ucfwd[0] = 0;
		  
    for (i = 1; i <= d->n; ++i)
	{
	    R_CheckUserInterrupt();
#ifdef VITDEBUG0
	    printf("obs %d\n", i);
#endif
	    if ((i < d->n) && (d->subject[i] == d->subject[i-1]))
		{
#ifdef VITDEBUG0
		    printf("subject %d\n", d->subject[i]);
#endif
		    dt = d->time[i] - d->time[i-1];
		    qmat = &(qm->intens[MI3(0, 0, i-1, qm->nst, qm->nst)]);
		    hpars = &(hm->pars[MI(0, i, hm->totpars)]); /* not i-1 as pre 1.2.3 */
		    if (d->nout > 1) outcome = &d->obs[MI(0, i, d->nout)];
		    else {
			GetCensored(d->obs[i], cm, &nc, &curr);
			outcome = curr;
		    }
		    GetOutcomeProb(pout, outcome, nc, d->nout, hpars, hm, qm, d->obstrue[i]);
#ifdef VITDEBUG0
		    for (tru=0;tru<nc;++tru) printf("outcome[%d] = %1.0lf, ",tru, outcome[tru]); printf("\n");
#endif
		    Pmat(pmat, dt, qmat, qm->nst,
			 (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);

		    psum = 0;
		    for (tru = 0; tru < qm->nst; ++tru)
			{
/* lvnew =  log prob of most likely path ending in tru at current obs.
   kmax  = most likely state at the previous obs  
   pfwd = p(data up to and including current obs, and hidden state at current obs), then scaled at each step to sum to 1 to avoid underflow.
*/
			    pfwd[MI(i,tru,d->n)] = 0;
			    for (k = 0; k < qm->nst; ++k) {
				lvp[k] = lvold[k] + log(pmat[MI(k, tru, qm->nst)]);
				pfwd[MI(i,tru,d->n)] += pfwd[MI(i-1,k,d->n)] * pmat[MI(k,tru,qm->nst)];
			    }
			    if (d->obstrue[i-1])
//				kmax = d->obs[MI(0, i-1, d->nout)] - 1;
				kmax = d->obstrue[i-1] - 1;
			    else pmax(lvp, qm->nst, &kmax);
			    lvnew[tru] = log ( pout[tru] )  +  lvp[kmax];
			    ptr[MI(i, tru, d->n)] = kmax;
			    pfwd[MI(i,tru,d->n)] *= pout[tru];
			    psum += pfwd[MI(i,tru,d->n)];
#ifdef VITDEBUG0
			    printf("true %d, pout[%d] = %lf, lvold = %lf, pmat = %lf, lvnew = %lf, ptr[%d,%d]=%d\n",
				   tru, tru, pout[tru], lvold[tru], pmat[MI(kmax, tru, qm->nst)], lvnew[tru], i, tru, ptr[MI(i, tru, d->n)]);
#endif
			}
		    ucfwd[i] = ucfwd[i-1] + log(psum);
		    for (k = 0; k < qm->nst; ++k){
			pfwd[MI(i,k,d->n)] /= psum;
			lvold[k] = lvnew[k];
		    }
		}
	    else
		{
		    /* Traceback and backward algorithm for current individual
		       pall = p(data at all times)
		       pbwd = p(data at future times | hidden state at current obs)
		       so that p(hidden state at current obs) = pfwd * pbwd / pall
		     */
		    pmax(lvold, qm->nst, &kmax);
		    obs = i-1;
//		    fitted[obs] = (d->obstrue[obs] ? d->obs[MI(0,obs,d->nout)]-1 : kmax);
		    fitted[obs] = (d->obstrue[obs] ? d->obstrue[obs]-1 : kmax);

		    logpall = 0;  // compute full likelihood.  
		    ucbwd[obs] = 0;
		    for (k = 0; k < qm->nst; ++k){
			pbwd[MI(obs,k,d->n)] = 1;
			logpall += pfwd[MI(obs,k,d->n)];
		    }
		    logpall = log(logpall) + ucfwd[obs];
		    for (k = 0; k < qm->nst; ++k){
			pstate[MI(obs,k,d->n)] = exp(log(pfwd[MI(obs,k,d->n)]) + log(pbwd[MI(obs,k,d->n)]) - logpall + ucfwd[obs]);
		    }
#ifdef VITDEBUG
		    printf("traceback for subject %d\n", d->subject[i-1]);
		    printf("obs=%d,fitted[%d]=%1.0lf\n",obs,obs,fitted[obs]);
		    for (tru = 0; tru < qm->nst; ++tru){
			printf("pfwd[%d,%d]=%f, ", obs, tru, pfwd[MI(obs,tru,d->n)]);
			printf("pbwd[%d,%d]=%f, ", obs, tru, pbwd[MI(obs,tru,d->n)]);
			printf("ucfwd[%d]=%f, ", obs, ucfwd[obs]);
			printf("ucbwd[%d]=%f, ", obs, ucfwd[obs]);
			printf("logpall=%f,",logpall);
			printf("pstate[%d,%d]=%f, ", obs, tru, pstate[MI(obs,tru,d->n)]);
			printf("\n");

		    }
#endif
		    while   ( (obs > 0) && (d->subject[obs] == d->subject[obs-1]) )
			{
			    fitted[obs-1] = ptr[MI(obs, fitted[obs], d->n)];
#ifdef VITDEBUG
			    printf("fitted[%d] = ptr[%d,%1.0lf] = %1.0lf\n", obs-1, obs, fitted[obs], fitted[obs-1]);
#endif

			    dt = d->time[obs] - d->time[obs-1];
			    qmat = &(qm->intens[MI3(0, 0, obs-1, qm->nst, qm->nst)]);
			    hpars = &(hm->pars[MI(0, obs, hm->totpars)]);
			    if (d->nout > 1) outcome = &d->obs[MI(0, obs, d->nout)];
			    else {
				GetCensored(d->obs[obs], cm, &nc, &curr);
				outcome = curr;
			    }
			    GetOutcomeProb(pout, outcome, nc, d->nout, hpars, hm, qm, d->obstrue[obs]);
			    Pmat(pmat, dt, qmat, qm->nst,
				 (d->obstype[obs] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);

			    psum=0;
			    for (tru = 0; tru < qm->nst; ++tru){
				pbwd[MI(obs-1,tru,d->n)] = 0;
				for (k = 0; k < qm->nst; ++k)
				    pbwd[MI(obs-1,tru,d->n)] += pmat[MI(tru,k,qm->nst)] * pout[k] * pbwd[MI(obs,k,d->n)];
				psum += pbwd[MI(obs-1,tru,d->n)];				
			    }
			    ucbwd[obs-1] = ucbwd[obs] + log(psum);
			    for (tru = 0; tru < qm->nst; ++tru){
				pbwd[MI(obs-1,tru,d->n)] /= psum;
				pstate[MI(obs-1,tru,d->n)] = exp(log(pfwd[MI(obs-1,tru,d->n)]) + log(pbwd[MI(obs-1,tru,d->n)]) - logpall + ucfwd[obs-1] + ucbwd[obs-1]);
#ifdef VITDEBUG
				printf("pfwd[%d,%d]=%f, ", obs-1, tru, pfwd[MI(obs-1,tru,d->n)]);
				printf("pbwd[%d,%d]=%f, ", obs-1, tru, pbwd[MI(obs-1,tru,d->n)]);
				printf("ucfwd[%d]=%f, ", obs, ucfwd[obs]);
				printf("ucbwd[%d]=%f, ", obs, ucfwd[obs]);
				printf("logpall=%f,",logpall);
				printf("pstate[%d,%d]=%f, ", obs-1, tru, pstate[MI(obs-1,tru,d->n)]);
				printf("\n");
#endif
			    }
#ifdef VITDEBUG
			    printf("logpall=%f\n",logpall);
#endif
			    --obs;
			}
#ifdef VITDEBUG
		    printf("\n");
#endif
		    if (i < d->n) {
			if (d->obstrue[i]) {
			    for (k = 0; k < qm->nst; ++k)
//				lvold[k] = (k+1 == d->obs[MI(0,i,d->nout)] ? 0 : R_NegInf);
				lvold[k] = (k+1 == d->obstrue[i] ? 0 : R_NegInf);
			}
			else {
			    if (d->nout > 1) outcome = &d->obs[MI(0, i, d->nout)];
			    else {
				GetCensored(d->obs[i], cm, &nc, &curr);
				outcome = curr;
			    }
			    /* initial observation is a censored state. No HMM here, so initprobs not needed */
			    if (nc > 1) {
				for (k = 0, j = 0; k < qm->nst; ++k) {
				    if (k+1 == outcome[j]) {
					lvold[k] = 0;
					++j;
				    }
				    else lvold[k] = R_NegInf;
				}
			    }
			    else { /* use initprobs */
				for (k = 0; k < qm->nst; ++k)
				    lvold[k] = log(hm->initp[MI(d->subject[i]-1, k, d->npts)]);
			    }
			}
			for (k = 0; k < qm->nst; ++k)
			    pfwd[MI(i,k,d->n)] = exp(lvold[k]); 
			ucfwd[i] = 0;
		    }
		}
	}
    Free(pmat); Free(ptr); Free(lvold); Free(lvnew); Free(lvp); Free(curr); Free(pout);
    Free(pfwd); Free(pbwd); Free(ucfwd); Free(ucbwd);
}

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

    for (int i = 0; i < length(list); i++)
	if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
    return elmt;
}

double *list_double_vec(SEXP list, const char *str) {
    SEXP elt = getListElement(list, str);
    return REAL(elt);
}

int *list_int_vec(SEXP list, const char *str) {
    SEXP elt = getListElement(list, str);
    return INTEGER(elt);
}

double list_double(SEXP list, const char *str) {
    SEXP elt = getListElement(list, str);
    return REAL(elt)[0];
}

int list_int(SEXP list, const char *str) {
    SEXP elt = getListElement(list, str);
    return INTEGER(elt)[0];
}

SEXP msmCEntry(SEXP do_what_s, SEXP mf_agg_s, SEXP mf_s, SEXP auxdata_s, SEXP qmodel_s, SEXP cmodel_s, SEXP hmodel_s, SEXP pars_s)
{
    msmdata d; qmodel qm; cmodel cm; hmodel hm;
    int do_what = INTEGER(do_what_s)[0], nopt;
    double lik, *ret, *fitted, *pfwd;
    SEXP ret_s;

/* type coercion for all these is done in R */
    d.fromstate = list_int_vec(mf_agg_s, "(fromstate)");
    d.tostate = list_int_vec(mf_agg_s, "(tostate)");
    d.timelag = list_double_vec(mf_agg_s, "(timelag)");
    d.nocc = list_int_vec(mf_agg_s, "(nocc)");
    d.noccsum = list_int_vec(mf_agg_s, "(noccsum)");
    d.whicha = list_int_vec(mf_agg_s, "(whicha)");
    d.obstypea = list_int_vec(mf_agg_s, "(obstype)");

    d.subject = list_int_vec(mf_s, "(subject)");
    d.time = list_double_vec(mf_s, "(time)");
    d.obs = list_double_vec(mf_s, "(state)");
    d.obstype = list_int_vec(mf_s, "(obstype)");
    d.obstrue = list_int_vec(mf_s, "(obstrue)");
    d.pcomb = list_int_vec(mf_s, "(pcomb)");

    d.nagg = list_int(auxdata_s, "nagg");
    d.n = list_int(auxdata_s, "n");
    d.npts = list_int(auxdata_s, "npts");
    d.ntrans = list_int(auxdata_s, "ntrans");
    d.npcombs = list_int(auxdata_s, "npcombs");
    d.firstobs = list_int_vec(auxdata_s, "firstobs");
    d.nout = list_int(auxdata_s, "nout");

    qm.nst = list_int(qmodel_s,"nstates");
    qm.npars = list_int(qmodel_s,"npars");
    qm.nopt = list_int(qmodel_s,"nopt");
    qm.iso = list_int(qmodel_s,"iso");
    qm.perm = list_int_vec(qmodel_s,"perm");
    qm.qperm = list_int_vec(qmodel_s,"qperm");
    qm.expm = list_int(qmodel_s,"expm");
    qm.nliks = list_int(auxdata_s,"nliks");
    qm.intens = list_double_vec(pars_s,"Q"); // TODO name
    qm.dintens = list_double_vec(pars_s,"DQ"); // TODO name

    cm.ncens = list_int(cmodel_s,"ncens");
    cm.censor = list_int_vec(cmodel_s,"censor");
    cm.states = list_int_vec(cmodel_s,"states");
    cm.index = list_int_vec(cmodel_s,"index");

    hm.hidden = list_int(hmodel_s,"hidden");
    hm.mv = list_int(hmodel_s,"mv");
    hm.ematrix = list_int(hmodel_s,"ematrix");
    hm.models = list_int_vec(hmodel_s,"models");
    hm.totpars = list_int(hmodel_s,"totpars");
    hm.npars = list_int_vec(hmodel_s,"npars");
    hm.firstpar = list_int_vec(hmodel_s,"firstpar");
    hm.pars = list_double_vec(pars_s,"H");
    hm.dpars = list_double_vec(pars_s,"DH");
    hm.nopt = list_int(hmodel_s,"nopt");
    hm.initp = list_double_vec(pars_s,"initprobs");
    nopt = list_int(pars_s, "nopt");

    if (do_what == DO_LIK) {
	msmLikelihood(&d, &qm, &cm, &hm, &lik);
	ret_s = ScalarReal(lik);
    }

    else if (do_what == DO_DERIV) {
	ret_s = PROTECT(NEW_NUMERIC(nopt));
	ret = REAL(ret_s);
	msmDeriv(&d, &qm, &cm, &hm, ret);
	UNPROTECT(1);
    }

    else if (do_what == DO_INFO) {
	ret_s = PROTECT(allocMatrix(REALSXP, nopt, nopt));
	ret = REAL(ret_s);
	msmInfo(&d, &qm, &cm, &hm, ret);
	UNPROTECT(1);
    }

    else if (do_what == DO_LIK_SUBJ) {
	ret_s = PROTECT(NEW_NUMERIC(d.npts));
	ret = REAL(ret_s);
	msmLikelihood_subj(&d, &qm, &cm, &hm, ret);
	UNPROTECT(1);
    }

    else if (do_what == DO_DERIV_SUBJ) {
	ret_s = PROTECT(allocMatrix(REALSXP, d.npts, nopt));
	ret = REAL(ret_s);
	msmDeriv_subj(&d, &qm, &cm, &hm, ret);
	UNPROTECT(1);
    }

    else if (do_what == DO_VITERBI) {
	/* Return a list of a vector (fitted states) and a matrix (hidden state probabilities) */
        /* see https://stat.ethz.ch/pipermail/r-devel/2013-April/066246.html */
	ret_s = PROTECT(allocVector(VECSXP, 2));
	SEXP fitted_s = SET_VECTOR_ELT(ret_s, 0, NEW_NUMERIC(d.n));
	SEXP pfwd_s = SET_VECTOR_ELT(ret_s, 1, allocMatrix(REALSXP, d.n, qm.nst));
	fitted = REAL(fitted_s);
	pfwd = REAL(pfwd_s);
	Viterbi(&d, &qm, &cm, &hm, fitted, pfwd);
	UNPROTECT(1);
    }

    else if (do_what == DO_DPMAT) {
	ret_s = PROTECT(alloc3DArray(REALSXP, d.ntrans, qm.nst, nopt));
	ret = REAL(ret_s);
	dpmat_obs(&d, &qm, &cm, &hm, ret);
	UNPROTECT(1);
    }

    else error("Unknown C task.\n");
    return ret_s;
}
