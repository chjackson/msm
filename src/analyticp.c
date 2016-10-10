/* Analytic formulae for P(t) in terms of transition intensities, for
selected 2, 3, 4 and 5-state models. These were derived using symbolic
algebra software (Mathematica).

Increases speed and stability by avoiding the numeric calculation of
the matrix exponential.

### The numbered label gives the indices into the matrix of rates (vectorised by reading across rows)
### e.g. the 3-state model with qmatrix of the form
### well-disease, well-death, disease-death transitions allowed.
### *,1,1
### 0,*,1
### 0,0,*    uses the function p3q124()

### Some models are isomorphic and use the same p?q? function.
### See .msm.graphs in R/constants.R for permutations.
*/

#include "msm.h"
#include <stdio.h>

void p2q1(Matrix pmat, double t, Matrix qmat, int *degen);
void p2q12(Matrix pmat, double t, Matrix qmat, int *degen);

void p3q12(Matrix pmat, double t, Matrix qmat, int *degen);
void p3q14(Matrix pmat, double t, Matrix qmat, int *degen);
void p3q16(Matrix pmat, double t, Matrix qmat, int *degen);
void p3q124(Matrix pmat, double t, Matrix qmat, int *degen);
void p3q135(Matrix pmat, double t, Matrix qmat, int *degen);
void p3q1246(Matrix pmat, double t, Matrix qmat, int *degen);

void p4q159(Matrix pmat, double t, Matrix qmat, int *degen);
void p4q13569(Matrix pmat, double t, Matrix qmat, int *degen);

void p5q1_6_11_16(Matrix pmat, double t, Matrix qmat, int *degen);
void p5q1_4_6_8_11_12_16(Matrix pmat, double t, Matrix qmat, int *degen);
void p5q1_6_7_11_12(Matrix pmat, double t, Matrix qmat, int *degen);

typedef void (*pfn)(Matrix pmat, double t, Matrix qmat, int *degen);

pfn P2FNS[] = {
    p2q1, p2q12
};

pfn P3FNS[] = {
    p3q12,p3q14,p3q16,p3q124,p3q135,p3q1246
};

pfn P4FNS[] = {
  p4q159,p4q13569
};

pfn P5FNS[] = {
  p5q1_6_11_16,p5q1_4_6_8_11_12_16,p5q1_6_7_11_12
};

void AnalyticP(Matrix pmat, double t, int nstates, int iso, int *perm, int *qperm, Matrix qmat, int *degen)
{
    int i, j;
    Matrix qmat_base = (Matrix) Calloc( (nstates)*(nstates), double);
    Matrix pmat_base = (Matrix) Calloc( (nstates)*(nstates), double);
    for (i=0; i<nstates; ++i) {
	for (j=0; j<nstates; ++j) {
	    qmat_base[MI(i,j,nstates)] = qmat[MI(qperm[i]-1,qperm[j]-1,nstates)];
	}
    }
    if (nstates==2)
	(P2FNS[iso-1])(pmat_base, t, qmat_base, degen);
    else if (nstates==3)
	(P3FNS[iso-1])(pmat_base, t, qmat_base, degen);
    else if (nstates==4)
	(P4FNS[iso-1])(pmat_base, t, qmat_base, degen);
    else if (nstates==5)
	(P5FNS[iso-1])(pmat_base, t, qmat_base, degen);
    else error("internal error in GetAnalyticP. Send a bug report to the package maintainer.");
    if (*degen) return;
    for (i=0; i<nstates; ++i)
	for (j=0; j<nstates; ++j) {
	    pmat[MI(i,j,nstates)] = pmat_base[MI(perm[i]-1,perm[j]-1,nstates)];
	}
    Free(pmat_base); Free(qmat_base);
}



/* TWO STATES */
/* 1 transition */

void p2q1(Matrix pmat, double t, Matrix qmat, int *degen)
{
    double a = qmat[MI(0,1,2)];
    double e1 = exp(-(a*t));
    pmat[MI(0,0,2)] = e1;
    pmat[MI(0,1,2)] = 1 - e1;
    pmat[MI(1,0,2)] = 0;
    pmat[MI(1,1,2)] = 1;
}

/* 2 transitions */

void p2q12(Matrix pmat, double t, Matrix qmat, int *degen)
{
    double a = qmat[MI(0,1,2)], b=qmat[MI(1,0,2)];
    double e1 = exp(-((a + b)*t));
    if (all_equal(a+b, 0)) {
	pmat[MI(0,0,2)] = 1; pmat[MI(1,1,2)] = 1;
	pmat[MI(0,1,2)] = 0; pmat[MI(1,0,2)] = 0;
    }
    else {
	pmat[MI(0,0,2)] = (b + a*e1)/(a + b);
	pmat[MI(0,1,2)] = (a - a*e1)/(a + b);
	pmat[MI(1,0,2)] = (b - b*e1)/(a + b);
	pmat[MI(1,1,2)] = (a + b*e1)/(a + b);
    }
}


/* THREE STATES */
/* 2 transitions */

void p3q12(Matrix pmat, double t, Matrix qmat, int *degen)
{
    double a = qmat[MI(0,1,3)], b=qmat[MI(0,2,3)];
    double e1 = exp(-(a + b)*t);
    pmat[MI(0,0,3)] = e1;
    if (all_equal(a+b, 0)) {
	pmat[MI(0,1,3)] = 0; pmat[MI(0,2,3)] = 0;
    }
    else {
	pmat[MI(0,1,3)] = (a - a*e1)/(a + b);
	pmat[MI(0,2,3)] = (b - b*e1)/(a + b);
    }
    pmat[MI(1,0,3)] = 0;
    pmat[MI(1,1,3)] = 1;
    pmat[MI(1,2,3)] = 0;
    pmat[MI(2,0,3)] = 0;
    pmat[MI(2,1,3)] = 0;
    pmat[MI(2,2,3)] = 1;
}

void p3q14(Matrix pmat, double t, Matrix qmat, int *degen)
{
    double a = qmat[MI(0,1,3)], b=qmat[MI(1,2,3)];
    double e1 = exp(-a*t);
    double e2 = exp(-b*t);
    pmat[MI(0,0,3)] = e1;
    pmat[MI(0,1,3)] = ( all_equal(a,b) ?
			a*t*e1 :
			(a*(e1 - e2)/(b - a)) );
    pmat[MI(0,2,3)] = ( all_equal(a,b) ?
			1 - e1 - a*t*e1 :
			1 - e1 - pmat[MI(0,1,3)] );
    /*			(a - a*e2 + b*(-1 + e1))/(a - b) ); */
    /*    printf("t=%f,a=%f,b=%f,p01=%f,p02=%f\n",t,a,b,pmat[MI(0,1,3)],pmat[MI(0,2,3)]); */
    pmat[MI(1,0,3)] = 0;
    pmat[MI(1,1,3)] = e2;
    pmat[MI(1,2,3)] = 1 - e2;
    pmat[MI(2,0,3)] = 0;
    pmat[MI(2,1,3)] = 0;
    pmat[MI(2,2,3)] = 1;
}

void p3q16(Matrix pmat, double t, Matrix qmat, int *degen)
{
    double a = qmat[MI(0,1,3)], b=qmat[MI(2,1,3)];
    double e1 = exp(-(a*t));
    double e2 = exp(-(b*t));
    pmat[MI(0,0,3)] = e1;
    pmat[MI(0,1,3)] = 1 - e1;
    pmat[MI(0,2,3)] = 0;
    pmat[MI(1,0,3)] = 0;
    pmat[MI(1,1,3)] = 1;
    pmat[MI(1,2,3)] = 0;
    pmat[MI(2,0,3)] = 0;
    pmat[MI(2,1,3)] = 1 - e2;
    pmat[MI(2,2,3)] = e2;
}

/* 3 transitions */

void p3q124(Matrix pmat, double t, Matrix qmat, int *degen)
{
    double a = qmat[MI(0,1,3)], b=qmat[MI(0,2,3)], c=qmat[MI(1,2,3)];
    double e1 = exp(-((a + b)*t));
    double e2 = exp(-c*t);
    pmat[MI(0,0,3)] = e1;
    pmat[MI(0,1,3)] = ( all_equal(a + b, c) ?
			(a*t)*e1 :
			(a*(-e1 + e2))/(a + b - c) );
    pmat[MI(0,2,3)] = ( all_equal(a + b, c) ?
			(-e1 + 1 - a*t*e1) :
			1 + (-b + c)*e1/(a + b - c) - a*e2/(a + b - c) );
    pmat[MI(1,0,3)] = 0;
    pmat[MI(1,1,3)] = e2;
    pmat[MI(1,2,3)] = 1 - e2;
    pmat[MI(2,0,3)] = 0;
    pmat[MI(2,1,3)] = 0;
    pmat[MI(2,2,3)] = 1;
}

void p3q135(Matrix pmat, double t, Matrix qmat, int *degen)
{
    double a = qmat[MI(0,1,3)], b=qmat[MI(1,0,3)], c=qmat[MI(2,0,3)];
    double e1 = exp(-(a + b)*t);
    double e2 = exp(-(c*t));
    double e3 = exp((a + b - c)*t);
    if (all_equal(a+b, 0)) {
	pmat[MI(0,0,3)] = 1; pmat[MI(1,1,3)] = 1;
	pmat[MI(0,1,3)] = 0; pmat[MI(1,0,3)] = 0;
    }
    else {
	pmat[MI(0,0,3)] = (b + a*e1)/(a + b);
	pmat[MI(0,1,3)] = (a - a*e1)/(a + b);
	pmat[MI(1,0,3)] = (b - b*e1)/(a + b);
	pmat[MI(1,1,3)] = (a + b*e1)/(a + b);
    }
    pmat[MI(0,2,3)] = 0;
    pmat[MI(1,2,3)] = 0;
    pmat[MI(2,0,3)] = (all_equal(a + b, c) ?
		       (pow(a,2)*t*e1 + b*(-e1 + 1 + a*t*e1))/
		       (a + b) :
		       (b*(b - c)*(-e2 + 1) +
			a*(c*e2 - c*e2/e3 + b*(-e2 + 1)))/
		       ((a + b)*(a + b - c)) );
    pmat[MI(2,1,3)] = (all_equal(a + b, c) ?
		       (a*(-e1 + 1 - (a+b)*e1*t))/(a + b) :
		       (a*(e1*c - c + (a + b)* (1 - e1*e3)))/
		       ((a + b)*(a + b - c)) );
    pmat[MI(2,2,3)] = e2;
}

/* 4 transitions */

void p3q1246(Matrix pmat, double t, Matrix qmat, int *degen)
{
    double a = qmat[MI(0,1,3)], b=qmat[MI(0,2,3)], c=qmat[MI(1,2,3)], d=qmat[MI(2,1,3)];
    double e1 = exp(-((a + b)*t));
    double e2 = exp(-((c + d)*t));
    pmat[MI(0,0,3)] = e1;
    pmat[MI(0,1,3)] = (all_equal( a + b, c + d) ?
		       (a + b - c)/(a + b) - (a + b - c)*e1/(a + b) +
		       ((-b + c)*t)*e1 :
		       (a*(d*(-1 + e1) +
			   c*(e1 - e2)) +
			d*(((c + d)*(-e1 + 1)) +
			   b*(-1 + e2)))/((c + d)*(-a - b + c + d)) );
    pmat[MI(0,2,3)] = (all_equal(a + b, c + d) ?
		       (b*(a + b)*t*e1 + c*(-e1 + 1 - a*e1*t - b*e1*t))/
		       ((a + b)) :
		       ((c*(c + d)*(-e1 + 1)) +
			a*c*(-1 + e2) +
			b*(c*(-1 + e1) +
			   d*(e1 - e2)))/
		       ((c + d)*(-a - b + c + d)) );
    pmat[MI(1,0,3)] = 0;
    pmat[MI(1,1,3)] = (d + c*e2)/(c + d);
    pmat[MI(1,2,3)] = (c - c*e2)/(c + d);
    pmat[MI(2,0,3)] = 0;
    pmat[MI(2,1,3)] = (d - d*e2)/(c + d);
    pmat[MI(2,2,3)] = (c + d*e2)/(c + d);
}

/* TODO: p3q1234, reversible illness-death model. No more complex than the five-state models below... */

/* FOUR STATES */

/* Progression only */

void p4q159(Matrix pmat, double t, Matrix qmat, int *degen){
    double a = qmat[MI(0,1,4)], b = qmat[MI(1,2,4)], c=qmat[MI(2,3,4)];
    double e1 = exp(-(a*t));
    double e2 = exp(-(b*t));
    double e3 = exp(-(c*t));
    pmat[MI(0,0,4)] = e1;
    pmat[MI(1,0,4)] = 0;
    pmat[MI(1,1,4)] = e2;
    pmat[MI(2,0,4)] = 0;
    pmat[MI(2,1,4)] = 0;
    pmat[MI(2,2,4)] = e3;
    pmat[MI(2,3,4)] = 1 - e3;
    pmat[MI(3,0,4)] = 0;
    pmat[MI(3,1,4)] = 0;
    pmat[MI(3,2,4)] = 0;
    pmat[MI(3,3,4)] = 1;

    if (all_equal(a,b) && !all_equal(b,c)){
      pmat[MI(0,1,4)] = (a*t)*e1;
      pmat[MI(0,2,4)] = -((pow(a,2)*(-e3 + e1*(1 + a*t - c*t)))/
			  (pow(a - c,2)));
      pmat[MI(0,3,4)] = 1 + ((2*a - c)*c*e1)/pow(a - c,2) -
	  pow(a,2)*e3/pow(a - c,2) + (a*c*t*e1)/(a - c);
      pmat[MI(1,2,4)] = -((a*(e1 - e3))/(a - c));
      pmat[MI(1,3,4)] = (a - a*e3 + c*(-1 + e1))/(a - c);
    }
    else if (all_equal(a,c) && !all_equal(b,c)){
      pmat[MI(0,1,4)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,4)] = -((a*b*(-e2 + e1*(1 + a*t - b*t)))/
			  (pow(a - b,2)));
      pmat[MI(0,3,4)] = 1 + ((2*a - b)*b*e1)/(pow(a - b,2)) -
	pow(a,2)*e2/(pow(a - b,2)) + (a*b*t*e1)/((a - b));
      pmat[MI(1,2,4)] = -((b*(e1 - e2))/(a - b));
      pmat[MI(1,3,4)] = (a - a*e2 + b*(-1 + e1))/(a - b);
    }
    else if (!all_equal(a,b) && all_equal(b,c)){
      pmat[MI(0,1,4)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,4)] = (a*b*(e1 + e2*(-1 + a*t - b*t)))/
	  (pow(a - b,2));
      pmat[MI(0,3,4)] = 1 - pow(b,2)*e1/(pow(a - b,2)) +
	(a*b*e2)/(pow(a - b,2)) - (a*(1 + b*t)*e2)/((a - b));
      pmat[MI(1,2,4)] = (b*t)*e2;
      pmat[MI(1,3,4)] = (-1 + 1/e2 - b*t)*e2;
    }
    else if (all_equal(a,b) && all_equal(b,c)){
      pmat[MI(0,1,4)] = (a*t)*e1;
      pmat[MI(0,2,4)] = (pow(a,2)*pow(t,2)*e1)/(2);
      pmat[MI(0,3,4)] = (-2*e1 + 2 - 2*e1*a*t - pow(a,2)*pow(t,2)*e1)/(2);
      pmat[MI(1,2,4)] = (a*t)*e1;
      pmat[MI(1,3,4)] = (-1 + 1/e1 - a*t)*e1;
    }
    else {
      pmat[MI(0,1,4)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,4)] = (a*b*(a*(e3 - e2) +
			      c*(e2 - e1) +
			      b*(-e3 + e1)))/
	  ((a - b)*(a - c)*(b - c));
      pmat[MI(0,3,4)] = 1 + (a*c*e2)/((a - b)*(b - c)) +
	  (b*(-(c*e1/(a - b)) + a*e3/(-b + c)))/(a - c);
      pmat[MI(1,2,4)] = -((b*(e2 - e3))/(b - c));
      pmat[MI(1,3,4)] = (b - b*e3 + c*(-1 + e2))/(b - c);
    }
}

/* Progression and death from every state */
/* TODO more rearrangements to avoid division by zero */

void p4q13569(Matrix pmat, double t, Matrix qmat, int *degen){
  double a = qmat[MI(0,1,4)], b = qmat[MI(0,3,4)], c=qmat[MI(1,2,4)], d=qmat[MI(1,3,4)], e=qmat[MI(2,3,4)];
  double e1 = exp(-((a + b)*t));
  double e2 = exp(-((c + d)*t));
  double e3 = exp(-(e*t));
  pmat[MI(0,0,4)] = e1;
  pmat[MI(1,0,4)] = 0;
  pmat[MI(1,1,4)] = e2;
  pmat[MI(2,0,4)] = 0;
  pmat[MI(2,1,4)] = 0;
  pmat[MI(2,2,4)] = e3;
  pmat[MI(2,3,4)] = 1 - e3;
  pmat[MI(3,0,4)] = 0;
  pmat[MI(3,1,4)] = 0;
  pmat[MI(3,2,4)] = 0;
  pmat[MI(3,3,4)] = 1;
  if (all_equal(a+b, c+d) && !all_equal(a+b, e)) {
    pmat[MI(0,1,4)] = (a*t)*e1;
    pmat[MI(0,2,4)] = (a*c*(-e1 + e3 + e1*(-a*t - b*t + e*t)))/
      (pow(a + b - e,2));
    pmat[MI(0,3,4)] = 1 - (a*(1/a - c/pow(a + b - e,2)))*e1 -
	(a*c*e3)/(pow(a + b - e,2)) -
	(a*(a + b - c - e)*t*e1)/(a + b - e);
    pmat[MI(1,2,4)] = (c*(-e1 + e3))/(a + b - e);
    pmat[MI(1,3,4)] = 1 + (-a - b + c + e)*e1/(a + b - e) -
      c*e3/(a + b - e);
  }
  else if (!all_equal(a+b, c+d) && all_equal(a+b, e)) {
    pmat[MI(0,1,4)] = (a*(e2 - e1))/(a + b - c - d);
    pmat[MI(0,2,4)] = -((a*c*(-1/e1 + 1/e2*
			      (1 + a*t + b*t - c*t - d*t)))/
			(pow(a + b - c - d,2) / (e1*e2) ));
    pmat[MI(0,3,4)] = 1 - (a*(b - 2*c - d) + pow(-b + c + d,2))/
      (pow(a + b - c - d,2)/e1) -
      (a*(a + b - d))/(pow(a + b - c - d,2)/e2) +
      (a*c*t)/((a + b - c - d)/e1);
    pmat[MI(1,2,4)] = (c*(e1 - e2))/(-a - b + c + d);
    pmat[MI(1,3,4)] = (a + b - c - d + c*e1 - a*e2 -
		       b*e2 + d*e2)/(a + b - c - d);
  }
  else if (!all_equal(a+b, c+d) && all_equal(c+d, e)) {
    pmat[MI(0,1,4)] = (a*(e2 - e1))/(a + b - c - d);
    pmat[MI(0,2,4)] = (a*c*(1/e2 + 1/e1*(-1 + a*t + b*t - c*t - d*t)))/
      (pow(a + b - c - d,2) / (e1*e2));
    pmat[MI(0,3,4)] = 1 - (a*(b - d) + pow(-b + c + d,2))/
      (pow(a + b - c - d,2)/e1) +
      (a*c)/(pow(a + b - c - d,2)/e2) -
      (a*(1 + c*t))/((a + b - c - d)/e2);
    pmat[MI(1,2,4)] = (c*t)*e2;
    pmat[MI(1,3,4)] = (-1 + 1/e2 - c*t)*e2;
  }
  else if (all_equal(a+b, c+d) && all_equal(a+b, e)) {
    pmat[MI(0,1,4)] = (a*t)*e1;
    pmat[MI(0,2,4)] = (a*c*pow(t,2)*e1)/(2);
    pmat[MI(0,3,4)] = (-2*e1 + 2 - a*t*e1*(2 + c*t))/(2);
    pmat[MI(1,2,4)] = (c*t)*e1;
    pmat[MI(1,3,4)] = (-e1 + 1 - c*t*e1);
  }
  else {
    pmat[MI(0,1,4)] = (a*(e2 - e1))/(a + b - c - d);
    pmat[MI(0,2,4)] = a*c*(e1/((a + b - c - d)*(a + b - e)) -
			   e2/((a + b - c - d)*(c + d - e)) -
			   e3/((a + b - e)*(-c - d + e)));
    pmat[MI(0,3,4)] = 1 - (a*(b - d) + (b - c - d)*(b - e))*e1/
	((a + b - c - d)*(a + b - e)) +
	(a*(-d + e)*e2)/((a + b - c - d)*(c + d - e)) -
	(a*c*e3)/((a + b - e)*(c + d - e));
    pmat[MI(1,2,4)] = (c*(-e2 + e3))/((c + d - e));
    pmat[MI(1,3,4)] = 1 + (-d + e)*e2/(c + d - e) - c*e3/(c + d - e);
  }
}


/* FIVE STATES */

/* Five states with progression only */

void p5q1_6_11_16(Matrix pmat, double t, Matrix qmat, int *degen){
  double a = qmat[MI(0,1,5)], b=qmat[MI(1,2,5)], c=qmat[MI(2,3,5)], d=qmat[MI(3,4,5)];
  double e1 = exp(-(a*t));
  double e2 = exp(-(b*t));
  double e3 = exp(-(c*t));
  double e4 = exp(-(d*t));
  pmat[MI(0,0,5)] = e1;
  pmat[MI(1,0,5)] = 0;
  pmat[MI(1,1,5)] = e2;
  pmat[MI(2,0,5)] = 0;
  pmat[MI(2,1,5)] = 0;
  pmat[MI(2,2,5)] = e3;
  pmat[MI(3,0,5)] = 0;
  pmat[MI(3,1,5)] = 0;
  pmat[MI(3,2,5)] = 0;
  pmat[MI(3,3,5)] = e4;
  pmat[MI(3,4,5)] = 1 - e4;
  pmat[MI(4,0,5)] = 0;
  pmat[MI(4,1,5)] = 0;
  pmat[MI(4,2,5)] = 0;
  pmat[MI(4,3,5)] = 0;
  pmat[MI(4,4,5)] = 1;
  if (all_equal(a,b) && !all_equal(a,c) && !all_equal(a,d) && !all_equal(c,d)){
      pmat[MI(0,1,5)] = (a*t)*e1;
      pmat[MI(0,2,5)] = -((pow(a,2)*(-1/e1 + 1/e3*(1 + a*t - c*t)))/
			  (pow(a - c,2)/(e1*e3)));
      pmat[MI(0,3,5)] = pow(a,2)*c*(-((-2*a + c + d)/
				      (pow(a - c,2)*pow(a - d,2)/e1)) -
				    1/(pow(a - c,2)*(c - d)/e3) +
				    1/(pow(a - d,2)*(c - d)/e4) +
				    t/((a - c)*(a - d)/e1));
      pmat[MI(0,4,5)] = 1 - (c*d*(3*pow(a,2) + c*d - 2*a*(c + d)))/
	  (pow(a - c,2)*pow(a - d,2)/e1) +
	  (pow(a,2)*d)/(pow(a - c,2)*(c - d)/e3) +
	  (pow(a,2)*c)/(pow(a - d,2)*(-c + d)/e4) -
	  (a*c*d*t)/((a - c)*(a - d)/e1);
      pmat[MI(1,2,5)] = -((a*(e1 - e3))/(a - c));
      pmat[MI(1,3,5)] = (a*c*(a*(1/(e1*e3) - 1/(e1*e4)) +
			      d*(1/(e1*e4) - 1/(e3*e4)) +
			      c*(-1/(e1*e3) + 1/(e3*e4))))/
	  ((a - c)*(a - d)*(c - d)/(e1*e3*e4));
      pmat[MI(1,4,5)] = 1 + (a*d)/((a - c)*(c - d)/e3) +
	  (c*(-(d/((a - c)/e1)) + a/((-c + d)/e4)))/(a - d);
      pmat[MI(2,3,5)] = -((c*(e3 - e4))/(c - d));
      pmat[MI(2,4,5)] = (c - c*e4 + d*(-1 + e3))/(c - d);
  }
  else if (all_equal(a,c) && !all_equal(a,b) && !all_equal(a,d) && !all_equal(b,d)){
      pmat[MI(0,1,5)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,5)] = -((a*b*(-1/e1 + 1/e2*(1 + a*t - b*t)))/
			  (pow(a - b,2)/(e1*e2)));
      pmat[MI(0,3,5)] = (pow(a,2)*b*(-(pow(d,2)*(1/(e1*e4) - 1/(e2*e4))) -
				     b*pow(d,2)*1/(e2*e4)*t +
				     pow(a,2)*(1/(e1*e2) - 1/(e1*e4) +
					       (b - d)*1/(e2*e4)*t) +
				     pow(b,2)/e2*(1/e1 + 1/e4*(-1 + d*t)) +
				     a*(-2*b*(1/(e1*e2) - 1/(e2*e4)) -
					pow(b,2)*1/(e2*e4)*t +
					d/e4*(2/e1 + 1/e2*(-2 + d*t)))))/
	  (pow(a - b,2)*pow(a - d,2)*(b - d)/(e1*e2*e4));
      pmat[MI(0,4,5)] = 1 - (b*d*(3*pow(a,2) + b*d - 2*a*(b + d)))/
	  (pow(a - b,2)*pow(a - d,2)/e1) +
	  (pow(a,2)*d)/(pow(a - b,2)*(b - d)/e2) +
	  (pow(a,2)*b)/(pow(a - d,2)*(-b + d)/e4) -
	  (a*b*d*t)/((a - b)*(a - d)/e1);
      pmat[MI(1,2,5)] = -((b*(e1 - e2))/(a - b));
      pmat[MI(1,3,5)] = (a*b*(a*(1/(e1*e2) - 1/(e1*e4)) +
			      d*(1/(e1*e4) - 1/(e2*e4)) +
			      b*(-1/(e1*e2) + 1/(e2*e4))))/
	  ((a - b)*(a - d)*(b - d)/(e1*e2*e4));
      pmat[MI(1,4,5)] = 1 + (a*d)/((a - b)*(b - d)/e2) +
	  (b*(-(d/((a - b)/e1)) + a/((-b + d)/e4)))/(a - d);
      pmat[MI(2,3,5)] = -((a*(e1 - e4))/(a - d));
      pmat[MI(2,4,5)] = (a - a*e4 + d*(-1 + e1))/(a - d);
  }
  else if (all_equal(a,d) && !all_equal(a,b) && !all_equal(a,c) && !all_equal(b,c)){
      pmat[MI(0,1,5)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,5)] = (a*b*(a*(1/(e1*e2) - 1/(e1*e3)) +
			      c*(1/(e1*e3) - 1/(e2*e3)) +
			      b*(-1/(e1*e2) + 1/(e2*e3))))/
	  ((a - b)*(a - c)*(b - c)/(e1*e2*e3));
      pmat[MI(0,3,5)] = (a*b*c*(-(pow(c,2)*(1/(e1*e3) - 1/(e2*e3))) -
				b*pow(c,2)*1/(e2*e3)*t +
				pow(a,2)*(1/(e1*e2) - 1/(e1*e3) +
					  (b - c)*1/(e2*e3)*t) +
				pow(b,2)/e2*(1/e1 + 1/e3*(-1 + c*t)) +
				a*(-2*b*(1/(e1*e2) - 1/(e2*e3)) -
				   pow(b,2)*1/(e2*e3)*t +
				   c/e3*(2/e1 + 1/e2*(-2 + c*t)))))/
	  (pow(a - b,2)*pow(a - c,2)*(b - c)/(e1*e2*e3));
      pmat[MI(0,4,5)] = 1 - (b*c*(3*pow(a,2) + b*c - 2*a*(b + c)))/
	  (pow(a - b,2)*pow(a - c,2)/e1) +
	  (pow(a,2)*c)/(pow(a - b,2)*(b - c)/e2) -
	  (pow(a,2)*b)/(pow(a - c,2)*(b - c)/e3) -
	  (a*b*c*t)/((a - b)*(a - c)/e1);
      pmat[MI(1,2,5)] = -((b*(e2 - e3))/(b - c));
      pmat[MI(1,3,5)] = (b*c*(a*(1/(e1*e2) - 1/(e1*e3)) +
			      c*(1/(e1*e3) - 1/(e2*e3)) +
			      b*(-1/(e1*e2) + 1/(e2*e3))))/
	  ((a - b)*(a - c)*(b - c)/(e1*e2*e3));
      pmat[MI(1,4,5)] = 1 + (a*c)/((a - b)*(b - c)/e2) +
	  (b*(-(c/((a - b)/e1)) + a/((-b + c)/e3)))/(a - c);
      pmat[MI(2,3,5)] = -((c*(e1 - e3))/(a - c));
      pmat[MI(2,4,5)] = (a - a*e3 + c*(-1 + e1))/(a - c);
  }
  else if (all_equal(b,c) && !all_equal(b,a) && !all_equal(b,d) && !all_equal(a,d)){
      pmat[MI(0,1,5)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,5)] = (a*b*(1/e2 + 1/e1*(-1 + a*t - b*t)))/
	  (pow(a - b,2)/(e1*e2));
      pmat[MI(0,3,5)] = a*pow(b,2)*(-(1/(pow(a - b,2)*(a - d)/e1)) +
				    1/((a - d)*pow(b - d,2)/e4) -
				    (a - 2*b + d + a*b*t - pow(b,2)*t - a*d*t + b*d*t)/
				    (pow(a - b,2)*pow(b - d,2)/e2));
      pmat[MI(0,4,5)] = 1 + (pow(b,2)*d)/(pow(a - b,2)*(a - d)/e1) -
	  (a*b*d)/(pow(a - b,2)*(b - d)/e2) -
	  (a*pow(b,2))/((a - d)*pow(b - d,2)/e4) +
	  (a*d*(-d + pow(b,2)*t + b*(2 - d*t)))/
	  ((a - b)*pow(b - d,2)/e2);
      pmat[MI(1,2,5)] = (b*t)*e2;
      pmat[MI(1,3,5)] = -((pow(b,2)*(-1/e2 + 1/e4*(1 + b*t - d*t)))/
			  (pow(b - d,2)*1/(e2*e4)));
      pmat[MI(1,4,5)] = 1 + ((2*b - d)*d)/(pow(b - d,2)/e2) -
	  pow(b,2)/(pow(b - d,2)/e4) + (b*d*t)/((b - d)/e2);
      pmat[MI(2,3,5)] = -((b*(e2 - e4))/(b - d));
      pmat[MI(2,4,5)] = (b - b*e4 + d*(-1 + e2))/(b - d);
  }
  else if (all_equal(b,d) && !all_equal(b,a) && !all_equal(b,c) && !all_equal(a,c)){
      pmat[MI(0,1,5)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,5)] = (a*b*(a*(1/(e1*e2) - 1/(e1*e3)) +
			      c*(1/(e1*e3) - 1/(e2*e3)) +
			      b*(-1/(e1*e2) + 1/(e2*e3))))/
	  ((a - b)*(a - c)*(b - c)/(e1*e2*e3));
      pmat[MI(0,3,5)] = a*b*c*(-(1/(pow(a - b,2)*(a - c)/e1)) +
			       1/((a - c)*pow(b - c,2)/e3) -
			       (a - 2*b + c + a*b*t - pow(b,2)*t - a*c*t + b*c*t)/
			       (pow(a - b,2)*pow(b - c,2)/e2));
      pmat[MI(0,4,5)] = 1 + (pow(b,2)*c)/(pow(a - b,2)*(a - c)/e1) -
	  (a*b*c)/(pow(a - b,2)*(b - c)/e2) -
	  (a*pow(b,2))/((a - c)*pow(b - c,2)/e3) +
	  (a*c*(-c + pow(b,2)*t + b*(2 - c*t)))/
	  ((a - b)*pow(b - c,2)/e2);
      pmat[MI(1,2,5)] = -((b*(e2 - e3))/(b - c));
      pmat[MI(1,3,5)] = -((b*c*(-1/e2 + 1/e3*(1 + b*t - c*t)))/
			  (pow(b - c,2)*1/(e2*e3)));
      pmat[MI(1,4,5)] = 1 + ((2*b - c)*c)/(pow(b - c,2)/e2) -
	  pow(b,2)/(pow(b - c,2)/e3) + (b*c*t)/((b - c)/e2);
      pmat[MI(2,3,5)] = -((c*(e2 - e3))/(b - c));
      pmat[MI(2,4,5)] = (b - b*e3 + c*(-1 + e2))/(b - c);
  }
  else if (all_equal(c,d) && !all_equal(c,a) && !all_equal(c,b) && !all_equal(a,b)){
      pmat[MI(0,1,5)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,5)] = (a*b*(a*(1/(e1*e2) - 1/(e1*e3)) +
			      c*(1/(e1*e3) - 1/(e2*e3)) +
			      b*(-1/(e1*e2) + 1/(e2*e3))))/
	  ((a - b)*(a - c)*(b - c)/(e1*e2*e3));
      pmat[MI(0,3,5)] = (a*b*c*(pow(c,2)*(1/(e1*e3) - 1/(e2*e3)) +
				a/e1*(2*c*(1/e2 - 1/e3) -
					    pow(b,2)/e2*t + pow(c,2)/e2*t) +
				pow(a,2)/e1*
				(1/e3 + 1/e2*(-1 + b*t - c*t)) +
				pow(b,2)/e2*(-1/e3 + 1/e1*(1 + c*t)) -
				b*c/e2*(-2/e3 + 1/e1*(2 + c*t))))/
	  ((a - b)*pow(a - c,2)*pow(b - c,2)/(e1*e2*e3));
      pmat[MI(0,4,5)] = (b*pow(b - c,2)*pow(c,2)*(1/(e2*e3) - 1/(e1*e2*e3)) +
			 pow(a,3)/e1*(pow(c,2)/e3*(-1 + 1/e2) +
					    pow(b,2)/e2*(-1 + 1/e3 - c*t) +
					    b*c/e2*(2 - 2/e3 + c*t)) -
			 a*c*(pow(c,3)*(1/(e1*e3) - 1/(e1*e2*e3)) -
			      pow(b,2)*c/(e1*e2)*(3 - 3/e3 + c*t) +
			      pow(b,3)/(e1*e2)*(2 - 2/e3 + c*t)) +
			 pow(a,2)*(-2*pow(c,3)/(e1*e3)*(-1 + 1/e2) -
				   b*pow(c,2)/(e1*e2)*(3 - 3/e3 + c*t) +
				   pow(b,3)/(e1*e2)*(1 - 1/e3 + c*t)))/
	  ((a - b)*pow(a - c,2)*pow(b - c,2)/(e1*e2*e3));
      pmat[MI(1,2,5)] = -((b*(e2 - e3))/(b - c));
      pmat[MI(1,3,5)] = (b*c*(1/e3 + 1/e2*(-1 + b*t - c*t)))/
	  (pow(b - c,2)*1/(e2*e3));
      pmat[MI(1,4,5)] = 1 - pow(c,2)/(pow(b - c,2)/e2) +
	  (b*c)/(pow(b - c,2)/e3) - (b*(1 + c*t))/((b - c)/e3);
      pmat[MI(2,3,5)] = (c*t)*e3;
      pmat[MI(2,4,5)] = (-1 + 1/e3 - c*t)*e3;
  }
  else if (all_equal(a,b) && all_equal(c,d) && !all_equal(a,c)){
    pmat[MI(0,1,5)] = (a*t)/exp(a*t);
    pmat[MI(0,2,5)] = -((pow(a,2)*(-exp(a*t) + exp(c*t)*(1 + a*t - c*t)))/
			(pow(a - c,2)*exp((a + c)*t)));
    pmat[MI(0,3,5)] = (pow(a,2)*c*(exp(a*t)*(-2 + a*t - c*t) + exp(c*t)*(2 + a*t - c*t)))/
      (pow(a - c,3)*exp((a + c)*t));
    pmat[MI(0,4,5)] = 1 + pow(c,3)/(pow(a - c,3)*exp(a*t)) -
      (a*pow(c,2)*(3 + a*t - c*t))/(pow(a - c,3)*exp(a*t)) -
      (pow(a,2)*(a + a*c*t - c*(3 + c*t)))/(pow(a - c,3)*exp(c*t));
    pmat[MI(1,2,5)] = -((a*(exp(-(a*t)) - exp(-(c*t))))/(a - c));
    pmat[MI(1,3,5)] = (a*c*(exp(c*t) + exp(a*t)*(-1 + a*t - c*t)))/
      (pow(a - c,2)*exp((a + c)*t));
    pmat[MI(1,4,5)] = 1 - pow(c,2)/(pow(a - c,2)*exp(a*t)) +
      (a*c)/(pow(a - c,2)*exp(c*t)) - (a*(1 + c*t))/((a - c)*exp(c*t));
    pmat[MI(2,3,5)] = (c*t)/exp(c*t);
    pmat[MI(2,4,5)] = (-1 + exp(c*t) - c*t)/exp(c*t);
  }
  else if (all_equal(a,c) && all_equal(b,d) && !all_equal(a,b)){
    pmat[MI(0,1,5)] = -((a*(exp(-(a*t)) - exp(-(b*t))))/(a - b));
    pmat[MI(0,2,5)] = -((a*b*(-exp(a*t) + exp(b*t)*(1 + a*t - b*t)))/
			(pow(a - b,2)*exp((a + b)*t)));
    pmat[MI(0,3,5)] = (pow(a,2)*b*(exp(a*t)*(-2 + a*t - b*t) + exp(b*t)*(2 + a*t - b*t)))/
      (pow(a - b,3)*exp((a + b)*t));
    pmat[MI(0,4,5)] = (-(pow(b,3)*exp(b*t)*(-1 + exp(a*t))) +
		       pow(a,3)*exp(a*t)*(-1 + exp(b*t) - b*t) +
		       a*pow(b,2)*exp(b*t)*(-3 + 3*exp(a*t) + b*t) +
		       pow(a,2)*b*(-3*exp((a + b)*t) - b*exp(b*t)*t +
				   exp(a*t)*(3 + b*t)))/(pow(a - b,3)*exp((a + b)*t));
    pmat[MI(1,2,5)] = -((b*(exp(-(a*t)) - exp(-(b*t))))/(a - b));
    pmat[MI(1,3,5)] = (a*b*(exp(b*t) + exp(a*t)*(-1 + a*t - b*t)))/
      (pow(a - b,2)*exp((a + b)*t));
    pmat[MI(1,4,5)] = 1 - pow(b,2)/(pow(a - b,2)*exp(a*t)) -
      (a*(a - 2*b))/(pow(a - b,2)*exp(b*t)) -
      (a*b*t)/((a - b)*exp(b*t));
    pmat[MI(2,3,5)] = -((a*(exp(-(a*t)) - exp(-(b*t))))/(a - b));
    pmat[MI(2,4,5)] = (a - a/exp(b*t) + b*(-1 + exp(-(a*t))))/(a - b);
  }
  else if (all_equal(a,d) && all_equal(b,c) && !all_equal(a,b)){
    pmat[MI(0,1,5)] = -((a*(exp(-(a*t)) - exp(-(b*t))))/(a - b));
    pmat[MI(0,2,5)] = (a*b*(exp(b*t) + exp(a*t)*(-1 + a*t - b*t)))/
      (pow(a - b,2)*exp((a + b)*t));
    pmat[MI(0,3,5)] = (a*pow(b,2)*(exp(a*t)*(-2 + a*t - b*t) + exp(b*t)*(2 + a*t - b*t)))/
      (pow(a - b,3)*exp((a + b)*t));
    pmat[MI(0,4,5)] = (-(pow(b,3)*exp(b*t)*(-1 + exp(a*t))) +
		       pow(a,3)*exp(a*t)*(-1 + exp(b*t) - b*t) +
		       a*pow(b,2)*exp(b*t)*(-3 + 3*exp(a*t) + b*t) +
		       pow(a,2)*b*(-3*exp((a + b)*t) - b*exp(b*t)*t +
				   exp(a*t)*(3 + b*t)))/(pow(a - b,3)*exp((a + b)*t));
    pmat[MI(1,2,5)] = (b*t)/exp(b*t);
    pmat[MI(1,3,5)] = (pow(b,2)*(exp(b*t) + exp(a*t)*(-1 + a*t - b*t)))/
      (pow(a - b,2)*exp((a + b)*t));
    pmat[MI(1,4,5)] = 1 - pow(b,2)/(pow(a - b,2)*exp(a*t)) -
      (a*(a - 2*b))/(pow(a - b,2)*exp(b*t)) -
      (a*b*t)/((a - b)*exp(b*t));
    pmat[MI(2,3,5)] = -((b*(exp(-(a*t)) - exp(-(b*t))))/(a - b));
    pmat[MI(2,4,5)] = (a - a/exp(b*t) + b*(-1 + exp(-(a*t))))/(a - b);
  }
  else if (all_equal(a,b) && all_equal(a,c) && !all_equal(a,d)){
      pmat[MI(0,1,5)] = (a*t)*e1;
      pmat[MI(0,2,5)] = (pow(a,2)*pow(t,2))/(2/e1);
      pmat[MI(0,3,5)] = -(pow(a,3)*(-2/e1 +
				    1/e4*(2 - 2*d*t + pow(a,2)*pow(t,2) +
					      pow(d,2)*pow(t,2) - 2*a*t*(-1 + d*t))))/
	  (2*pow(a - d,3)*1/(e1*e4));
      pmat[MI(0,4,5)] = 1 + (d*(3*pow(a,2) - 3*a*d + pow(d,2)))/(pow(a - d,3)/e1) -
	  pow(a,3)/(pow(a - d,3)/e4) +
	  (a*(2*a - d)*d*t)/(pow(a - d,2)/e1) +
	  (pow(a,2)*d*pow(t,2))/((2*a - 2*d)/e1);
      pmat[MI(1,2,5)] = (a*t)*e1;
      pmat[MI(1,3,5)] = -((pow(a,2)*(-1/e1 + 1/e4*(1 + a*t - d*t)))/
			  (pow(a - d,2)*1/(e1*e4)));
      pmat[MI(1,4,5)] = (pow(d,2)/e4*(-1 + 1/e1) -
			 a*d/e4*(-2 + 2/e1 + d*t) +
			 pow(a,2)*(-1/e1 + 1/(e1*e4) + d/e4*t))/
	  (pow(a - d,2)*1/(e1*e4));
      pmat[MI(2,3,5)] = -((a*(e1 - e4))/(a - d));
      pmat[MI(2,4,5)] = (a - a*e4 + d*(-1 + e1))/(a - d);
  }
  else if (all_equal(a,b) && all_equal(a,d) && !all_equal(a,c)){
      pmat[MI(0,1,5)] = (a*t)*e1;
      pmat[MI(0,2,5)] = -((pow(a,2)*(-1/e1 + 1/e3*(1 + a*t - c*t)))/
			  (pow(a - c,2)/(e1*e3)));
      pmat[MI(0,3,5)] = -(pow(a,2)*c*(-2/e1 +
				      1/e3*(2 - 2*c*t + pow(a,2)*pow(t,2) +
						pow(c,2)*pow(t,2) - 2*a*t*(-1 + c*t))))/
	  (2*pow(a - c,3)/(e1*e3));
      pmat[MI(0,4,5)] = 1 + (c*(3*pow(a,2) - 3*a*c + pow(c,2)))/(pow(a - c,3)/e1) -
	  pow(a,3)/(pow(a - c,3)/e3) +
	  (a*(2*a - c)*c*t)/(pow(a - c,2)/e1) +
	  (pow(a,2)*c*pow(t,2))/((2*a - 2*c)/e1);
      pmat[MI(1,2,5)] = -((a*(e1 - e3))/(a - c));
      pmat[MI(1,3,5)] = -((a*c*(-1/e1 + 1/e3*(1 + a*t - c*t)))/
			  (pow(a - c,2)/(e1*e3)));
      pmat[MI(1,4,5)] = (pow(c,2)/e3*(-1 + 1/e1) -
			 a*c/e3*(-2 + 2/e1 + c*t) +
			 pow(a,2)*(-1/e1 + 1/(e1*e3) + c/e3*t))/
	  (pow(a - c,2)/(e1*e3));
      pmat[MI(2,3,5)] = -((c*(e1 - e3))/(a - c));
      pmat[MI(2,4,5)] = (a - a*e3 + c*(-1 + e1))/(a - c);
  }
  else if (all_equal(a,c) && all_equal(a,d) && !all_equal(a,b)){
      pmat[MI(0,1,5)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,5)] = -((a*b*(-1/e1 + 1/e2*(1 + a*t - b*t)))/
			  (pow(a - b,2)/(e1*e2)));
      pmat[MI(0,3,5)] = -(pow(a,2)*b*(-2/e1 +
				      1/e2*(2 - 2*b*t + pow(a,2)*pow(t,2) +
						pow(b,2)*pow(t,2) - 2*a*t*(-1 + b*t))))/
	  (2*pow(a - b,3)/(e1*e2));
      pmat[MI(0,4,5)] = 1 + (b*(3*pow(a,2) - 3*a*b + pow(b,2)))/(pow(a - b,3)/e1) -
	  pow(a,3)/(pow(a - b,3)/e2) +
	  (a*(2*a - b)*b*t)/(pow(a - b,2)/e1) +
	  (pow(a,2)*b*pow(t,2))/((2*a - 2*b)/e1);
      pmat[MI(1,2,5)] = -((b*(e1 - e2))/(a - b));
      pmat[MI(1,3,5)] = -((a*b*(-1/e1 + 1/e2*(1 + a*t - b*t)))/
			  (pow(a - b,2)/(e1*e2)));
      pmat[MI(1,4,5)] = (pow(b,2)/e2*(-1 + 1/e1) -
			 a*b/e2*(-2 + 2/e1 + b*t) +
			 pow(a,2)*(-1/e1 + 1/(e1*e2) + b/e2*t))/
	  (pow(a - b,2)/(e1*e2));
      pmat[MI(2,3,5)] = (a*t)*e1;
      pmat[MI(2,4,5)] = (-1 + 1/e1 - a*t)*e1;
  }
  else if (all_equal(b,c) && all_equal(b,d) && !all_equal(b,a)){
      pmat[MI(0,1,5)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,5)] = (a*b*(1/e2 + 1/e1*(-1 + a*t - b*t)))/
	  (pow(a - b,2)/(e1*e2));
      pmat[MI(0,3,5)] = (a*pow(b,2)*(-2/e2 +
				     1/e1*(2 + 2*b*t + pow(a,2)*pow(t,2) +
					       pow(b,2)*pow(t,2) - 2*a*t*(1 + b*t))))/
	  (2*pow(a - b,3)/(e1*e2));
      pmat[MI(0,4,5)] = 1 + pow(b,3)/(pow(a - b,3)/e1) -
	  (a*pow(b,2))/(pow(a - b,3)/e2) +
	  (a*b*(1 + b*t))/(pow(a - b,2)/e2) -
	  (a*(2 + 2*b*t + pow(b,2)*pow(t,2)))/(2*(a - b)/e2);
      pmat[MI(1,2,5)] = (b*t)*e2;
      pmat[MI(1,3,5)] = (pow(b,2)*pow(t,2))/(2/e2);
      pmat[MI(1,4,5)] = (-2 + 2/e2 - 2*b*t - pow(b,2)*pow(t,2))/(2/e2);
      pmat[MI(2,3,5)] = (b*t)*e2;
      pmat[MI(2,4,5)] = (-1 + 1/e2 - b*t)*e2;
  }
  else if (all_equal(a,b) && all_equal(a,c) && all_equal(a,d)){
      pmat[MI(0,1,5)] = (a*t)*e1;
      pmat[MI(0,2,5)] = (pow(a,2)*pow(t,2))/(2/e1);
      pmat[MI(0,3,5)] = (pow(a,3)*pow(t,3))/(6./e1);
      pmat[MI(0,4,5)] = (-6 + 6/e1 - 6*a*t - 3*pow(a,2)*pow(t,2) -
			 pow(a,3)*pow(t,3))/(6./e1);
      pmat[MI(1,2,5)] = (a*t)*e1;
      pmat[MI(1,3,5)] = (pow(a,2)*pow(t,2))/(2/e1);
      pmat[MI(1,4,5)] = (-2 + 2/e1 - 2*a*t - pow(a,2)*pow(t,2))/(2/e1);
      pmat[MI(2,3,5)] = (a*t)*e1;
      pmat[MI(2,4,5)] = (-1 + 1/e1 - a*t)*e1;
  }
  else {
      pmat[MI(0,1,5)] = -((a*(e1 - e2))/(a - b));
      pmat[MI(0,2,5)] = (a*b*(a*(1/(e1*e2) - 1/(e1*e3)) +
			      c*(1/(e1*e3) - 1/(e2*e3)) +
			      b*(-1/(e1*e2) + 1/(e2*e3))))/
	  ((a - b)*(a - c)*(b - c)/(e1*e2*e3));
      pmat[MI(0,3,5)] = a*b*c*(-(1/((a - b)*(a - c)*(a - d)/e1)) +
			       1/((a - b)*(b - c)*(b - d)/e2) +
			       1/((a - c)*(-b + c)*(c - d)/e3) +
			       1/((a - d)*(-b + d)*(-c + d)/e4));
      pmat[MI(0,4,5)] = 1 - (a*c*d)/((a - b)*(b - c)*(b - d)/e2) +
	  b*((a*d)/((a - c)*(b - c)*(c - d)/e3) +
	     (c*(d/((a - b)*(a - c)/e1) +
		 a/((c - d)*(-b + d)/e4)))/(a - d));
      pmat[MI(1,2,5)] = -((b*(e2 - e3))/(b - c));
      pmat[MI(1,3,5)] = (b*c*(b*(1/(e2*e3) - 1/(e2*e4)) +
			      d*(1/(e2*e4) - 1/(e3*e4)) +
			      c*(-1/(e2*e3) + 1/(e3*e4))))/
	  ((b - c)*(b - d)*(c - d)/(e2*e3*e4));
      pmat[MI(1,4,5)] = 1 + (b*d)/((b - c)*(c - d)/e3) +
	  (c*(-(d/((b - c)/e2)) + b/((-c + d)/e4)))/(b - d);
      pmat[MI(2,3,5)] = -((c*(e3 - e4))/(c - d));
      pmat[MI(2,4,5)] = (c - c*e4 + d*(-1 + e3))/(c - d);
  }
}

/* Five states with progression and death from every state */

void p5q1_4_6_8_11_12_16(Matrix pmat, double t, Matrix qmat, int *degen){
    double a = qmat[MI(0,1,5)], b=qmat[MI(0,4,5)], c=qmat[MI(1,2,5)], d=qmat[MI(1,4,5)], e=qmat[MI(2,3,5)], f=qmat[MI(2,4,5)], g=qmat[MI(3,4,5)];
    double e1 = exp(-((a + b)*t));
    double e2 = exp(-((c + d)*t));
    double e3 = exp(-((e + f)*t));
    double e4 = exp(-(g*t));
    pmat[MI(0,0,5)] = e1;
    pmat[MI(1,0,5)] = 0;
    pmat[MI(1,1,5)] = e2;
    pmat[MI(2,0,5)] = 0;
    pmat[MI(2,1,5)] = 0;
    pmat[MI(2,2,5)] = e3;
    pmat[MI(3,0,5)] = 0;
    pmat[MI(3,1,5)] = 0;
    pmat[MI(3,2,5)] = 0;
    pmat[MI(3,3,5)] = e4;
    pmat[MI(3,4,5)] = 1 - e4;
    pmat[MI(4,0,5)] = 0;
    pmat[MI(4,1,5)] = 0;
    pmat[MI(4,2,5)] = 0;
    pmat[MI(4,3,5)] = 0;
    pmat[MI(4,4,5)] = 1;
    if (all_equal(a+b,c+d) && !all_equal(a+b,e+f) && !all_equal(a+b,g)){
	pmat[MI(0,1,5)] = (a*t)*e1;
	pmat[MI(0,2,5)] = -((a*c*(-1/e1 + 1/e3*
				  (1 + a*t + b*t - e*t - f*t)))/
			    (1/(e1*e3)*pow(a + b - e - f,2)));
	pmat[MI(0,3,5)] = a*c*e*(-(1/(1/e3*pow(a + b - e - f,2)*(e + f - g))) +
				 1/(1/e4*pow(a + b - g,2)*(e + f - g)) -
				 (-2*a - 2*b + e + f + g)/
				 (1/e1*pow(a + b - e - f,2)*pow(a + b - g,2)) +
				 t/(1/e1*(a + b - e - f)*(a + b - g)));
	pmat[MI(0,4,5)] = 1 - (a*c*e)/(1/e4*pow(a + b - g,2)*(e + f - g)) -
	    (a*c*(f - g))/(1/e3*pow(a + b - e - f,2)*(e + f - g)) +
	    (a*(-(1/a) + (c*(pow(a,2) + pow(b,2) + pow(e,2) + e*f +
			     2*a*(b - e - g) + e*g + pow(g,2) - 2*b*(e + g)))/
		(pow(a + b - e - f,2)*pow(a + b - g,2))))*e1 -
	    (a*(pow(a,2) + pow(b,2) + c*e + a*(2*b - c - e - f - g) + c*g + e*g +
		f*g - b*(c + e + f + g))*t)/
	    (1/e1*(a + b - e - f)*(a + b - g));
	pmat[MI(1,2,5)] = (c*(-e1 + e3))/(a + b - e - f);
	pmat[MI(1,3,5)] = c*e*(1/(1/e1*(a + b - e - f)*(a + b - g)) +
			       1/(1/e3*(-a - b + e + f)*(e + f - g)) -
			       1/(1/e4*(e + f - g)*(-a - b + g)));
	pmat[MI(1,4,5)] = 1 + (c*(f - g))/(1/e3*(-a - b + e + f)*(e + f - g)) +
	    (c*e)/(1/e4*(e + f - g)*(-a - b + g)) +
	    (-pow(a,2) - pow(b,2) - c*e - c*g - e*g - f*g + b*(c + e + f + g) +
	     a*(-2*b + c + e + f + g))/
	    (1/e1*(a + b - e - f)*(a + b - g));
	pmat[MI(2,3,5)] = (e*(-1 + e4/e3))/(1/e3*(e + f - g));
	pmat[MI(2,4,5)] = 1 - e/(1/e4*(e + f - g)) + (-f + g)/(1/e3*(e + f - g));
    }
    else if (all_equal(a+b,e+f) && !all_equal(a+b,c+d) && !all_equal(a+b,g)){
	pmat[MI(0,1,5)] = (a*(e2 - e1))/(a + b - c - d);
	pmat[MI(0,2,5)] = -((a*c*(-1/e1 + 1/e2*
				  (1 + a*t + b*t - c*t - d*t)))/
			    (pow(a + b - c - d,2)/(e1*e2)));
	pmat[MI(0,3,5)] = a*c*e*((2*a + 2*b - c - d - g)/
				 (pow(a + b - c - d,2)/e1*pow(a + b - g,2)) -
				 1/(pow(a + b - c - d,2)/e2*(c + d - g)) +
				 1/(1/e4*pow(a + b - g,2)*(c + d - g)) +
				 t/((a + b - c - d)/e1*(a + b - g)));
	pmat[MI(0,4,5)] = 1 - (a*c*e)/(1/e4*pow(a + b - g,2)*(c + d - g)) -
	    (a*(-(c*d) - pow(d,2) - c*e + a*(c + d - g) + b*(c + d - g) + d*g))/
	    (pow(a + b - c - d,2)/e2*(c + d - g)) -
	    (pow(a,3)*(b - 2*c - d) + pow(-b + c + d,2)*pow(b - g,2) +
	     pow(a,2)*(3*pow(b,2) + pow(c,2) - 2*b*(3*c + 2*d + g) +
		       d*(d + 2*g) + 2*c*(d + e + 2*g)) +
	     a*(3*pow(b,3) - d*g*(2*d + g) - pow(c,2)*(e + 2*g) -
		pow(b,2)*(6*c + 5*d + 4*g) - c*(g*(e + 2*g) + d*(e + 4*g)) +
		b*(2*pow(c,2) + 2*pow(d,2) + 6*d*g + pow(g,2) +
		   2*c*(2*d + e + 4*g))))/
	    (pow(a + b - c - d,2)/e1*pow(a + b - g,2)) +
	    (a*c*(a + b - e - g)*t)/((a + b - c - d)/e1*(a + b - g));
	pmat[MI(1,2,5)] = (c*(e1 - e2))/(-a - b + c + d);
	pmat[MI(1,3,5)] = c*e*(1/((a + b - c - d)/e1*(a + b - g)) -
			       1/((a + b - c - d)/e2*(c + d - g)) -
			       1/(1/e4*(a + b - g)*(-c - d + g)));
	pmat[MI(1,4,5)] = 1 - (c*(a + b - e - g))/((-a - b + c + d)/e1*(a + b - g)) +
	    (c*e)/(1/e4*(c + d - g)*(-a - b + g)) +
	    (c*d + pow(d,2) + c*e - a*(c + d - g) - b*(c + d - g) - d*g)/
	    ((a + b - c - d)/e2*(c + d - g));
	pmat[MI(2,3,5)] = (e*(-e1 + e4))/(a + b - g);
	pmat[MI(2,4,5)] = (e - e*(e1/e4) -
			   ((-1 + 1/e1)*(a + b - g))*(e1/e4))/
	    (1/e4*(-a - b + g));
    }
    else if (all_equal(a+b,g) && !all_equal(a+b,c+d) && !all_equal(a+b,e+f)){
	pmat[MI(0,1,5)] = (a*(e2 - e1))/(a + b - c - d);
	pmat[MI(0,2,5)] = a*c*(1/((a + b - c - d)/e1*(a + b - e - f)) -
			       1/((a + b - c - d)/e2*(c + d - e - f)) +
			       1/(1/e3*(-a - b + e + f)*(-c - d + e + f)));
	pmat[MI(0,3,5)] = a*c*e*((2*a + 2*b - c - d - e - f)/
				 (pow(a + b - c - d,2)/e1*pow(a + b - e - f,2)) -
				 1/(pow(a + b - c - d,2)/e2*(c + d - e - f)) -
				 1/(1/e3*pow(a + b - e - f,2)*(-c - d + e + f)) +
				 t/((a + b - c - d)/e1*(a + b - e - f)));
	pmat[MI(0,4,5)] = 1 - (a*c*(a + b - f))/
	    (1/e3*pow(a + b - e - f,2)*(c + d - e - f)) -
	    (a*(-(c*d) - pow(d,2) + d*e + a*(d - e - f) + b*(d - e - f) + c*f +
		d*f))/(pow(a + b - c - d,2)/e2*(c + d - e - f)) +
	    (a*c*e*(1/(-a - b + c + d) +
		    (a*(-b + d) - (b - c - d)*(b - e - f))/(a*c*e) + 1/(-a - b + e + f)))/
	    ((a + b - c - d)/e1*(a + b - e - f)) -
	    (a*c*e*t)/((a + b - c - d)/e1*(a + b - e - f));
	pmat[MI(1,2,5)] = (c*(-e2 + e3))/(c + d - e - f);
	pmat[MI(1,3,5)] = c*e*(1/((a + b - c - d)/e1*(a + b - e - f)) -
			       1/((a + b - c - d)/e2*(c + d - e - f)) +
			       1/(1/e3*(-a - b + e + f)*(-c - d + e + f)));
	pmat[MI(1,4,5)] = 1 + (c*e)/((-a - b + c + d)/e1*(a + b - e - f)) +
	    (c*(a + b - f))/(1/e3*(c + d - e - f)*(-a - b + e + f)) +
	    (c*d + pow(d,2) - d*e - c*f - d*f + a*(-d + e + f) + b*(-d + e + f))/
	    ((a + b - c - d)/e2*(c + d - e - f));
	pmat[MI(2,3,5)] = (e*(e1 - e3))/(-a - b + e + f);
	pmat[MI(2,4,5)] = (a + b - e + e*e1 - a*e3 -
			   b*e3 - f + f*e3)/(a + b - e - f);
    }
    else if (all_equal(c+d,e+f) && !all_equal(c+d,a+b) && !all_equal(c+d,g)){
	pmat[MI(0,1,5)] = (a*(e2 - e1))/(a + b - c - d);
	pmat[MI(0,2,5)] = (a*c*(1/e2 + 1/e1*(-1 + a*t + b*t - c*t - d*t)))/
	    (pow(a + b - c - d,2)/(e1*e2));
	pmat[MI(0,3,5)] = a*c*e*(-(1/(pow(a + b - c - d,2)/e1*(a + b - g))) +
				 1/(1/e4*(a + b - g)*pow(c + d - g,2)) +
				 1/(pow(a + b - c - d,2)/e2*(c + d - g)) -
				 (1 + c*t + d*t - g*t)/
				 ((a + b - c - d)/e2*pow(c + d - g,2)));
	pmat[MI(0,4,5)] = 1 - (a*c*e)/(1/e4*(a + b - g)*pow(c + d - g,2)) +
	    (a*c*(c + d - e - g))/
	    (pow(a + b - c - d,2)/e2*(c + d - g)) -
	    (pow(a,2)*(b - d) + pow(-b + c + d,2)*(b - g) +
	     a*(2*pow(b,2) + pow(c,2) + c*(2*d - e) + d*(d + g) -
		b*(2*c + 3*d + g)))/
	    (pow(a + b - c - d,2)/e1*(a + b - g)) +
	    (a*c*(-(1/c) + e/pow(c + d - g,2) - ((c + d - e - g)*t)/(c + d - g)))/
	    ((a + b - c - d)/e2);
	pmat[MI(1,2,5)] = (c*t)*e2;
	pmat[MI(1,3,5)] = (c*e*(-1 + e4/e2 - c*t - d*t + g*t))/
	    (1/e2*pow(c + d - g,2));
	pmat[MI(1,4,5)] = 1 + (c*(-(1/c) + e/pow(c + d - g,2)))*e2 -
	    (c*e)/(1/e4*pow(c + d - g,2)) -
	    (c*(c + d - e - g)*t)/(1/e2*(c + d - g));
	pmat[MI(2,3,5)] = (e*(-e2 + e4))/(c + d - g);
	pmat[MI(2,4,5)] = 1 - e/(1/e4*(c + d - g)) +
	    (-c - d + e + g)/(1/e2*(c + d - g));
    }
    else if (all_equal(c+d,g) && !all_equal(c+d,a+b) && !all_equal(c+d,e+f)){
	pmat[MI(0,1,5)] = (a*(e2 - e1))/(a + b - c - d);
	pmat[MI(0,2,5)] = a*c*(1/((a + b - c - d)/e1*(a + b - e - f)) -
			       1/((a + b - c - d)/e2*(c + d - e - f)) +
			       1/(1/e3*(-a - b + e + f)*(-c - d + e + f)));
	pmat[MI(0,3,5)] = a*c*e*(-(1/(pow(a + b - c - d,2)/e1*(a + b - e - f))) +
				 1/(1/e3*(a + b - e - f)*pow(c + d - e - f,2)) +
				 1/(pow(a + b - c - d,2)/e2*(c + d - e - f)) +
				 (-1 - c*t - d*t + e*t + f*t)/
				 ((a + b - c - d)/e2*pow(c + d - e - f,2)));
	pmat[MI(0,4,5)] = 1 - (a*c*(c + d - f))/
	    (1/e3*(a + b - e - f)*pow(c + d - e - f,2)) -
	    (a*c*e)/(pow(a + b - c - d,2)/e2*(c + d - e - f)) -
	    (pow(a,2)*(b - d) + pow(-b + c + d,2)*(b - e - f) +
	     a*(2*pow(b,2) + c*(d + f) + d*(d + e + f) - b*(2*c + 3*d + e + f)))/
	    (pow(a + b - c - d,2)/e1*(a + b - e - f)) +
	    (a*(-pow(-d + e + f,2) + pow(c,2)*e*t +
		c*(2*e + f - pow(e,2)*t - e*f*t + d*(-1 + e*t))))/
	    ((a + b - c - d)/e2*pow(c + d - e - f,2));
	pmat[MI(1,2,5)] = (c*(-e2 + e3))/(c + d - e - f);
	pmat[MI(1,3,5)] = -((c*e*(-1/e2 + 1/e3*
				  (1 + c*t + d*t - e*t - f*t)))/
			    (1/(e2*e3)*pow(c + d - e - f,2)));
	pmat[MI(1,4,5)] = 1 - (c*(c + d - f))/(1/e3*pow(c + d - e - f,2)) -
	    (c*(d - 2*e - f) + pow(-d + e + f,2))/
	    (1/e2*pow(c + d - e - f,2)) +
	    (c*e*t)/(1/e2*(c + d - e - f));
	pmat[MI(2,3,5)] = (e*(e2 - e3))/(-c - d + e + f);
	pmat[MI(2,4,5)] = (c + d - e + e*e2 - c*e3 -
			   d*e3 - f + f*e3)/(c + d - e - f);
    }
    else if (all_equal(e+f,g) && !all_equal(e+f,a+b) && !all_equal(e+f,c+d)){
	pmat[MI(0,1,5)] = (a*(e2 - e1))/(a + b - c - d);
	pmat[MI(0,2,5)] = a*c*(1/((a + b - c - d)/e1*(a + b - e - f)) -
			       1/((a + b - c - d)/e2*(c + d - e - f)) +
			       1/(1/e3*(-a - b + e + f)*(-c - d + e + f)));
	pmat[MI(0,3,5)] = a*c*e*(-(1/((a + b - c - d)/e1*pow(a + b - e - f,2))) +
				 1/((a + b - c - d)/e2*pow(c + d - e - f,2)) -
				 (a + b + c + d - 2*e - 2*f)/
				 (1/e3*pow(a + b - e - f,2)*pow(c + d - e - f,2)) +
				 t/(1/e3*(-a - b + e + f)*(-c - d + e + f)));
	pmat[MI(0,4,5)] = 1 + (a*c*e*(a + b + c + d - 2*e - 2*f))/
	    (1/e3*pow(a + b - e - f,2)*pow(c + d - e - f,2)) -
	    (a*(c*(d - f) + pow(-d + e + f,2)))/
	    ((a + b - c - d)/e2*pow(c + d - e - f,2)) -
	    (pow(a,2)*(b - d) + (b - c - d)*pow(-b + e + f,2) +
	     a*(2*pow(b,2) + c*f + 2*d*(e + f) - b*(c + 2*(d + e + f))))/
	    ((a + b - c - d)/e1*pow(a + b - e - f,2)) -
	    (a*c*(1 + e*t))/(1/e3*(a + b - e - f)*(c + d - e - f));
	pmat[MI(1,2,5)] = (c*(-e2 + e3))/(c + d - e - f);
	pmat[MI(1,3,5)] = (c*e*(1/e3 + 1/e2*(-1 + c*t + d*t - e*t - f*t)))/
	    (1/(e2*e3)*pow(c + d - e - f,2));
	pmat[MI(1,4,5)] = 1 + (c*e)/(1/e3*pow(c + d - e - f,2)) -
	    (c*(d - f) + pow(-d + e + f,2))/
	    (1/e2*pow(c + d - e - f,2)) -
	    (c*(1 + e*t))/(1/e3*(c + d - e - f));
	pmat[MI(2,3,5)] = (e*t)*e3;
	pmat[MI(2,4,5)] = (-1 + 1/e3 - e*t)*e3;
    }
    else if (all_equal(a+b,c+d) && all_equal(e+f,g) && !all_equal(a+b,e+f)){
      pmat[MI(0,1,5)] = (a*t)/exp((a + b)*t);
      pmat[MI(0,2,5)] = -((a*c*(-exp((a + b)*t) + exp((e + f)*t)*
				(1 + a*t + b*t - e*t - f*t)))/
			  (exp((a + b + e + f)*t)*pow(a + b - e - f,2)));
      pmat[MI(0,3,5)] = (a*c*e*(exp((a + b)*t)*(-2 + a*t + b*t - e*t - f*t) +
				exp((e + f)*t)*(2 + a*t + b*t - e*t - f*t)))/
	(exp((a + b + e + f)*t)*pow(a + b - e - f,3));

      pmat[MI(0,4,5)] = 1 + (a*(-(1/a) + (c*(a + b - 3*e - f))/pow(a + b - e - f,3)))/
	exp((a + b)*t) + (2*a*c*e)/
	(exp((e + f)*t)*pow(a + b - e - f,3)) -
	(a*(pow(a,2) + pow(b,2) + 2*c*e + pow(e,2) + c*f + 2*e*f +
	    pow(f,2) + a*(2*b - c - 2*(e + f)) - b*(c + 2*(e + f)))*t)/
	(exp((a + b)*t)*pow(a + b - e - f,2)) -
	(a*c*(1 + e*t))/(exp((e + f)*t)*pow(a + b - e - f,2));
      pmat[MI(1,2,5)] = (c*(-exp(-((a + b)*t)) + exp(-((e + f)*t))))/(a + b - e - f);
      pmat[MI(1,3,5)] = -((c*e*(-exp((e + f)*t) + exp((a + b)*t)*
				(1 - a*t - b*t + e*t + f*t)))/
			  (exp((a + b + e + f)*t)*pow(a + b - e - f,2)));
      pmat[MI(1,4,5)] = 1 + (c*e)/(exp((e + f)*t)*pow(a + b - e - f,2)) -
	(pow(a,2) + pow(b,2) + 2*c*e + pow(e,2) + c*f + 2*e*f + pow(f,2) +
	 a*(2*b - c - 2*(e + f)) - b*(c + 2*(e + f)))/
	(exp((a + b)*t)*pow(a + b - e - f,2)) +
	(c*(1 + e*t))/(exp((e + f)*t)*(-a - b + e + f));
      pmat[MI(2,3,5)] = (e*t)/exp((e + f)*t);
      pmat[MI(2,4,5)] = (-1 + exp((e + f)*t) - e*t)/exp((e + f)*t);
    }
    else if (all_equal(a+b,e+f) && all_equal(c+d,g) && !all_equal(a+b,c+d)){
      pmat[MI(0,1,5)] = (a*(-exp(-((a + b)*t)) + exp(-((c + d)*t))))/(a + b - c - d);
      pmat[MI(0,2,5)] = -((a*c*(-exp((a + b)*t) + exp((c + d)*t)*
				(1 + a*t + b*t - c*t - d*t)))/
			  (pow(a + b - c - d,2)*exp((a + b + c + d)*t)));
      pmat[MI(0,3,5)] = (a*c*e*(exp((a + b)*t)*(-2 + a*t + b*t - c*t - d*t) +
				exp((c + d)*t)*(2 + a*t + b*t - c*t - d*t)))/
	(pow(a + b - c - d,3)*exp((a + b + c + d)*t));
      pmat[MI(0,4,5)] = 1 - (pow(a,2)*(b - 2*c - d) + pow(b - c - d,3) +
			     a*(2*pow(b,2) - 5*b*c + 3*pow(c,2) - 4*b*d + 5*c*d + 2*pow(d,2) +
				2*c*e))/(pow(a + b - c - d,3)*exp((a + b)*t)) +
	(a*c*e)/(pow(a + b - c - d,3)*exp((c + d)*t)) +
	(a*c*(a + b - c - d - e)*t)/(pow(a + b - c - d,2)*exp((a + b)*t)) -
	(a*(pow(a,2) + pow(b,2) + c*d + pow(d,2) - c*e - pow(c,2)*e*t -
	    c*d*e*t + b*(-c - 2*d + c*e*t) + a*(2*b - c - 2*d + c*e*t)))/
	(pow(a + b - c - d,3)*exp((c + d)*t));
      pmat[MI(1,2,5)] = (c*(exp(-((a + b)*t)) - exp(-((c + d)*t))))/(-a - b + c + d);
      pmat[MI(1,3,5)] = -((c*e*(-exp((c + d)*t) + exp((a + b)*t)*
				(1 - a*t - b*t + c*t + d*t)))/
			  (pow(a + b - c - d,2)*exp((a + b + c + d)*t)));
      pmat[MI(1,4,5)] = 1 - (c*(-a - b + c + d + e))/(pow(a + b - c - d,2)*exp((a + b)*t)) -
	(pow(a,2) + pow(b,2) + a*(2*b - c - 2*d) + c*d + pow(d,2) -
	 b*(c + 2*d) - c*e)/(pow(a + b - c - d,2)*exp((c + d)*t)) +
	(c*e*t)/((-a - b + c + d)*exp((c + d)*t));
      pmat[MI(2,3,5)] = (e*(-exp(-((a + b)*t)) + exp(-((c + d)*t))))/(a + b - c - d);
      pmat[MI(2,4,5)] = (a - a/exp((a + b)*t) + ((b - c - d)*exp((c + d)*t)*
						 (-1 + exp((a + b)*t)) +
						 e*(-exp((a + b)*t) + exp((c + d)*t)))/
			 exp((a + b + c + d)*t))/(a + b - c - d);
    }
    else if (all_equal(a+b,g) && all_equal(c+d,e+f) && !all_equal(a+b,c+d)){
      pmat[MI(0,1,5)] = (a*(-exp(-((a + b)*t)) + exp(-((c + d)*t))))/(a + b - c - d);
      pmat[MI(0,2,5)] = (a*c*(exp((c + d)*t) + exp((a + b)*t)*(-1 + a*t + b*t - c*t - d*t)))/
	(pow(a + b - c - d,2)*exp((a + b + c + d)*t));
      pmat[MI(0,3,5)] = (a*c*e*(exp((a + b)*t)*(-2 + a*t + b*t - c*t - d*t) +
				exp((c + d)*t)*(2 + a*t + b*t - c*t - d*t)))/
	(pow(a + b - c - d,3)*exp((a + b + c + d)*t));
      pmat[MI(0,4,5)] = 1 - (pow(a,2)*(b - d) + pow(b - c - d,3) +
			     a*(2*pow(b,2) - 3*b*c + pow(c,2) - 4*b*d + 3*c*d + 2*pow(d,2) +
				2*c*e))/(pow(a + b - c - d,3)*exp((a + b)*t)) +
	(a*c*(a + b - c - d + e))/(pow(a + b - c - d,3)*exp((c + d)*t)) -
	(a*c*e*t)/(pow(a + b - c - d,2)*exp((a + b)*t)) +
	(a*(-pow(a,2) - 2*a*b - pow(b,2) + 2*a*c + 2*b*c - pow(c,2) +
	    2*a*d + 2*b*d - 2*c*d - pow(d,2) + c*e +
	    c*(a + b - c - d)*(-a - b + c + d - e)*t))/
	(pow(a + b - c - d,3)*exp((c + d)*t));
      pmat[MI(1,2,5)] = (c*t)/exp((c + d)*t);
      pmat[MI(1,3,5)] = -((c*e*(-exp((c + d)*t) + exp((a + b)*t)*
				(1 - a*t - b*t + c*t + d*t)))/
			  (pow(a + b - c - d,2)*exp((a + b + c + d)*t)));
      pmat[MI(1,4,5)] = 1 - (c*e)/(pow(a + b - c - d,2)*exp((a + b)*t)) -
	(pow(a,2) + pow(b,2) + pow(c,2) + 2*a*(b - c - d) + 2*c*d +
	 pow(d,2) - 2*b*(c + d) - c*e)/
	(pow(a + b - c - d,2)*exp((c + d)*t)) -
	(c*(-a - b + c + d - e)*t)/((-a - b + c + d)*exp((c + d)*t));
      pmat[MI(2,3,5)] = (e*(-exp(-((a + b)*t)) + exp(-((c + d)*t))))/(a + b - c - d);
      pmat[MI(2,4,5)] = (a + b - c - d + e/exp((a + b)*t) - a/exp((c + d)*t) -
			 b/exp((c + d)*t) + c/exp((c + d)*t) + d/exp((c + d)*t) -
			 e/exp((c + d)*t))/(a + b - c - d);
    }
    else if (all_equal(a+b,c+d) && all_equal(a+b,e+f) && !all_equal(a+b,g)){
	pmat[MI(0,1,5)] = (a*t)*e1;
	pmat[MI(0,2,5)] = (a*c*pow(t,2))/(2/e1);
	pmat[MI(0,3,5)] = (a*c*e*(-2*e1 + 2*e4 -
				  (2*(a + b - g)*t)*e1 -
				  (pow(a + b - g,2)*pow(t,2))*e1))/
	    (2*pow(a + b - g,3));
	pmat[MI(0,4,5)] = 1 + (-1 + (a*c*e)/pow(a + b - g,3))*e1 -
	    (a*c*e)/(1/e4*pow(a + b - g,3)) -
	    (a*(pow(a,2) + pow(b,2) - c*e + 2*a*(b - g) - 2*b*g + pow(g,2))*t)/
	    (1/e1*pow(a + b - g,2)) -
	    (a*c*(a + b - e - g)*pow(t,2))/(2/e1*(a + b - g));
	pmat[MI(1,2,5)] = (c*t)*e1;
	pmat[MI(1,3,5)] = (c*e*(-1 + e4/e1 - a*t - b*t + g*t))/
	    (1/e1*pow(a + b - g,2));
	pmat[MI(1,4,5)] = 1 - (c*e)/(1/e4*pow(a + b - g,2)) +
	    (-pow(a,2) - 2*a*b - pow(b,2) + c*e + 2*a*g + 2*b*g - pow(g,2) -
	     c*(a + b - g)*(a + b - e - g)*t)/(1/e1*pow(a + b - g,2));
	pmat[MI(2,3,5)] = (e*(-e1 + e4))/(a + b - g);
	pmat[MI(2,4,5)] = 1 - e/(1/e4*(a + b - g)) +
	    (-a - b + e + g)/(1/e1*(a + b - g));
    }
    else if (all_equal(a+b,c+d) && all_equal(a+b,g) && !all_equal(a+b,e+f)){
	pmat[MI(0,1,5)] = (a*t)*e1;
	pmat[MI(0,2,5)] = -((a*c*(-1/e1 + 1/e3*
				  (1 + a*t + b*t - e*t - f*t)))/
			    (1/(e1*e3)*pow(a + b - e - f,2)));
	pmat[MI(0,3,5)] = (a*c*e*(-2*e1 + 2*e3 -
				  (2*(a + b - e - f)*t)*e1 -
				  (pow(a + b - e - f,2)*pow(t,2))*e1))/
	    (2*pow(a + b - e - f,3));
	pmat[MI(0,4,5)] = 1 - (a*c*(a + b - f))/(1/e3*pow(a + b - e - f,3)) -
	    (pow(a,3) + pow(b - e - f,3) + pow(a,2)*(3*b - c - 3*(e + f)) +
	     a*(3*pow(b,2) + 3*pow(e,2) + 6*e*f + f*(c + 3*f) -
		b*(c + 6*(e + f))))/(1/e1*pow(a + b - e - f,3)) -
	    (a*(pow(a,2) + pow(b,2) + pow(e,2) + c*f + 2*e*f + pow(f,2) +
		a*(2*b - c - 2*(e + f)) - b*(c + 2*(e + f)))*t)/
	    (1/e1*pow(a + b - e - f,2)) +
	    (a*c*e*pow(t,2))/(2/e1*(a + b - e - f));
	pmat[MI(1,2,5)] = (c*(-e1 + e3))/(a + b - e - f);
	pmat[MI(1,3,5)] = (c*e*(1/e1 + 1/e3*(-1 - a*t - b*t + e*t + f*t)))/
	    (1/(e1*e3)*pow(a + b - e - f,2));
	pmat[MI(1,4,5)] = 1 - (c*(a + b - f))/(1/e3*pow(a + b - e - f,2)) -
	    (pow(a,2) + pow(b,2) + pow(e,2) + c*f + 2*e*f + pow(f,2) +
	     c*pow(e,2)*t + c*e*f*t + a*(2*b - c - 2*e - 2*f - c*e*t) -
	     b*(c + 2*(e + f) + c*e*t))/(1/e1*pow(a + b - e - f,2));
	pmat[MI(2,3,5)] = (e*(e1 - e3))/(-a - b + e + f);
	pmat[MI(2,4,5)] = (a + b - e + e*e1 - a*e3 -
			   b*e3 - f + f*e3)/(a + b - e - f);
    }
    else if (all_equal(a+b,e+f) && all_equal(a+b,g) && !all_equal(a+b,c+d)){
	pmat[MI(0,1,5)] = (a*(e2 - e1))/(a + b - c - d);
	pmat[MI(0,2,5)] = -((a*c*(-1/e1 + 1/e2*
				  (1 + a*t + b*t - c*t - d*t)))/
			    (pow(a + b - c - d,2)/(e1*e2)));
	pmat[MI(0,3,5)] = (a*c*e*(-2*e1 + 2*e2 -
				  (2*(a + b - c - d)*t)*e1 -
				  (pow(a + b - c - d,2)*pow(t,2))*e1))/
	    (2*pow(a + b - c - d,3));
	pmat[MI(0,4,5)] = 1 - (pow(a,2)*(b - 2*c - d) + pow(b - c - d,3) +
			       a*(2*pow(b,2) - 5*b*c + 3*pow(c,2) - 4*b*d + 5*c*d + 2*pow(d,2) -
				  c*e))/(pow(a + b - c - d,3)/e1) -
	    (a*(pow(a,2) + pow(b,2) + a*(2*b - c - 2*d) + c*d + pow(d,2) -
		b*(c + 2*d) + c*e))/(pow(a + b - c - d,3)/e2) +
	    (a*c*(a + b - c - d + e)*t)/(pow(a + b - c - d,2)/e1) +
	    (a*c*e*pow(t,2))/(2*(a + b - c - d)/e1);
	pmat[MI(1,2,5)] = (c*(e1 - e2))/(-a - b + c + d);
	pmat[MI(1,3,5)] = (c*e*(1/e1 + 1/e2*(-1 - a*t - b*t + c*t + d*t)))/
	    (pow(a + b - c - d,2)/(e1*e2));
	pmat[MI(1,4,5)] = 1 - (pow(a,2) + pow(b,2) + a*(2*b - c - 2*d) + c*d + pow(d,2) -
			       b*(c + 2*d) + c*e)/(pow(a + b - c - d,2)/e2) -
	    (c*(c + d - e + c*e*t + d*e*t - a*(1 + e*t) - b*(1 + e*t)))/
	    (pow(a + b - c - d,2)/e1);
	pmat[MI(2,3,5)] = (e*t)*e1;
	pmat[MI(2,4,5)] = (-1 + 1/e1 - e*t)*e1;
    }
    else if (all_equal(c+d,e+f) && all_equal(c+d,g) && !all_equal(a+b,c+d)){
	pmat[MI(0,1,5)] = (a*(e2 - e1))/(a + b - c - d);
	pmat[MI(0,2,5)] = (a*c*(1/e2 + 1/e1*(-1 + a*t + b*t - c*t - d*t)))/
	    (pow(a + b - c - d,2)/(e1*e2));
	pmat[MI(0,3,5)] = (a*c*e*(-2*e1 + 2*e2 -
				  (2*(a + b - c - d)*t)*e2 +
				  (pow(a + b - c - d,2)*pow(t,2))*e2))/
	    (2*pow(a + b - c - d,3));
	pmat[MI(0,4,5)] = 1 + (pow(a,2)*(b - d) + pow(b - c - d,3) +
			       a*(2*pow(b,2) - 3*b*c + pow(c,2) - 4*b*d + 3*c*d + 2*pow(d,2) -
				  c*e))/(pow(-a - b + c + d,3)/e1) -
	    (a*c*e)/(pow(a + b - c - d,3)/e2) +
	    (a*c*(1 + e*t))/(pow(a + b - c - d,2)/e2) -
	    (a*(2 + c*t*(2 + e*t)))/(2*(a + b - c - d)/e2);
	pmat[MI(1,2,5)] = (c*t)*e2;
	pmat[MI(1,3,5)] = (c*e*pow(t,2))/(2/e2);
	pmat[MI(1,4,5)] = (-2 + 2/e2 - c*t*(2 + e*t))/(2/e2);
	pmat[MI(2,3,5)] = (e*t)*e2;
	pmat[MI(2,4,5)] = (-1 + 1/e2 - e*t)*e2;
    }
    else if (all_equal(a+b,c+d) && all_equal(a+b,e+f) && all_equal(a+b,g)){
	pmat[MI(0,1,5)] = (a*t)*e1;
	pmat[MI(0,2,5)] = (a*c*pow(t,2))/(2/e1);
	pmat[MI(0,3,5)] = (a*c*e*pow(t,3))/(6./e1);
	pmat[MI(0,4,5)] = (-6 + 6/e1 - a*t*(6 + c*t*(3 + e*t)))/(6./e1);
	pmat[MI(1,2,5)] = (c*t)*e1;
	pmat[MI(1,3,5)] = (c*e*pow(t,2))/(2/e1);
	pmat[MI(1,4,5)] = (-2 + 2/e1 - c*t*(2 + e*t))/(2/e1);
	pmat[MI(2,3,5)] = (e*t)*e1;
	pmat[MI(2,4,5)] = (-1 + 1/e1 - e*t)*e1;
    }
    else {
	pmat[MI(0,1,5)] = (a*(e2 - e1))/(a + b - c - d);
	pmat[MI(0,2,5)] = a*c*(1/((a + b - c - d)/e1*(a + b - e - f)) -
			       1/((a + b - c - d)/e2*(c + d - e - f)) +
			       1/(1/e3*(-a - b + e + f)*(-c - d + e + f)));
	pmat[MI(0,3,5)] = a*c*e*(-(1/((a + b - c - d)/e1*(a + b - e - f)*
				      (a + b - g))) + 1/
				 ((a + b - c - d)/e2*(c + d - e - f)*(c + d - g)) -
				 1/(1/e3*(-a - b + e + f)*(-c - d + e + f)*(e + f - g)) +
				 1/(1/e4*(a + b - g)*(c + d - g)*(e + f - g)));
	pmat[MI(0,4,5)] = 1 - (a*(c*(d - f) + (d - e - f)*(d - g)))/
	    ((a + b - c - d)/e2*(c + d - e - f)*(c + d - g)) -
	    (a*c*e)/(1/e4*(a + b - g)*(c + d - g)*(e + f - g)) -
	    (a*c*(f - g))/
	    (1/e3*(a + b - e - f)*(c + d - e - f)*(e + f - g)) -
	    (pow(a,2)*(b - d) + (b - c - d)*(b - e - f)*(b - g) +
	     a*(2*pow(b,2) + c*f + d*(e + f + g) - b*(c + 2*d + e + f + g)))/
	    ((a + b - c - d)/e1*(a + b - e - f)*(a + b - g));
	pmat[MI(1,2,5)] = (c*(-e2 + e3))/(c + d - e - f);
	pmat[MI(1,3,5)] = c*e*(1/(1/e2*(c + d - e - f)*(c + d - g)) +
			       1/(1/e3*(-c - d + e + f)*(e + f - g)) -
			       1/(1/e4*(c + d - g)*(-e - f + g)));
	pmat[MI(1,4,5)] = 1 - (c*(d - f) + (d - e - f)*(d - g))/
	    (1/e2*(c + d - e - f)*(c + d - g)) -
	    (c*e)/(1/e4*(c + d - g)*(e + f - g)) +
	    (c*(-f + g))/(1/e3*(c + d - e - f)*(e + f - g));
	pmat[MI(2,3,5)] = (e*(-1 + e4/e3))/(1/e3*(e + f - g));
	pmat[MI(2,4,5)] = 1 - e/(1/e4*(e + f - g)) + (-f + g)/(1/e3*(e + f - g));
    }
}

/* Five states 0,1,2,3,4 with adjacent progression, and 1-3, 2-4 jumps */

void p5q1_6_7_11_12(Matrix pmat, double t, Matrix qmat, int *degen){
  double a = qmat[MI(0,1,5)], b=qmat[MI(1,2,5)], c=qmat[MI(1,3,5)], d=qmat[MI(2,3,5)], e=qmat[MI(2,4,5)];
  double e1 = exp(-(a*t));
  double e2 = exp(-((b + c)*t));
  double e3 = exp(-((d + e)*t));
  pmat[MI(0,0,5)] = e1;
  pmat[MI(1,0,5)] = 0;
  pmat[MI(1,1,5)] = e2;
  pmat[MI(2,0,5)] = 0;
  pmat[MI(2,1,5)] = 0;
  pmat[MI(2,2,5)] = e3;
  pmat[MI(3,0,5)] = 0;
  pmat[MI(3,1,5)] = 0;
  pmat[MI(3,2,5)] = 0;
  pmat[MI(3,3,5)] = 1;
  pmat[MI(3,4,5)] = 0;
  pmat[MI(4,0,5)] = 0;
  pmat[MI(4,1,5)] = 0;
  pmat[MI(4,2,5)] = 0;
  pmat[MI(4,3,5)] = 0;
  pmat[MI(4,4,5)] = 1;
  if (all_equal(a,b+c) && !all_equal(a,d+e)){
      pmat[MI(0,1,5)] = (a*t)*e1;
      pmat[MI(0,2,5)] = (a*b*(-e1 + e3 + e1*(-a*t + d*t + e*t)))/
	  (pow(-a + d + e,2));
      pmat[MI(0,3,5)] = (-(b*e) + a*(d + e))/(a*(d + e)) +
	  (-pow(a,3) + b*e*(d + e) - a*(pow(d,2) + 2*d*e + e*(2*b + e)) +
	   pow(a,2)*(b + 2*(d + e)))*e1/(a*pow(-a + d + e,2)) -
	  (a*b*d)*e3/((d + e)*pow(-a + d + e,2)) -
	  ((pow(a,2) + b*e - a*(b + d + e))*t)*e1/((a - d - e));
      pmat[MI(0,4,5)] = (b*e*(pow(d + e,2)*(-1 + 1/e1) -
			      a*(d + e)*(-2 + 2/e1 + d*t + e*t) +
			      pow(a,2)*(1/e1 - e3/e1 + (d + e)*t)))/
	  (a*(d + e)*pow(-a + d + e,2)/e1);
      pmat[MI(1,2,5)] = (b*(-e1 + e3))/((a - d - e));
      pmat[MI(1,3,5)] = (-(b*e) + a*(d + e))/(a*(d + e)) +
	  (-pow(a,2) - b*e + a*(b + d + e))*e1/(a*(a - d - e)) +
	  (b*d)/((d + e)*(-a + d + e)/e3);
      pmat[MI(1,4,5)] = (b*e*(d + e - d/e1 - e/e1 +
			      a*(1/e1 - e3/e1)))/
	  (a*(a - d - e)*(d + e)/e1);
      pmat[MI(2,3,5)] = (d - d*e3)/(d + e);
      pmat[MI(2,4,5)] = (e - e*e3)/(d + e); 
  }
  else if (!all_equal(a,b+c) && all_equal(a,d+e)){
      pmat[MI(0,1,5)] = (a*(-1 + e2/e1))/((a - b - c)/e1);
      pmat[MI(0,2,5)] = (a*b*(-1 + e2/e1 - a*t + b*t + c*t))/
	  (pow(-a + b + c,2)/e1);
      pmat[MI(0,3,5)] = (b*pow(b + c,2)*d*(-1 + 1/e1) +
			 pow(a,3)*c*(1/e1 - e2/e1) -
			 a*(b + c)*(-(pow(c,2)*(-1 + 1/e1)) + pow(b,2)*d*t +
				    b*(c - 2*d - c/e1 + 2*d/e1 + c*d*t)) +
			 pow(a,2)*(pow(c,2)*(1 - 2/e1 + e2/e1) +
				   pow(b,2)*d*t + b*(d*(1/e1 - e2/e1) +
						     c*(1 - 2/e1 + e2/e1 + d*t))))/
	  (a*(b + c)*pow(-a + b + c,2)/e1);
      pmat[MI(0,4,5)] = -((b*(a - d)*(-(pow(b + c,2)*(-1 + 1/e1)) +
				      a*(b + c)*(-2 + 2/e1 + b*t + c*t) -
				      pow(a,2)*(1/e1 - e2/e1 + (b + c)*t)))/
			  (a*(b + c)*pow(-a + b + c,2)/e1));
      pmat[MI(1,2,5)] = -((b*(-1 + e2/e1))/((-a + b + c)/e1));
      pmat[MI(1,3,5)] = (a*c + b*d)/(a*b + a*c) - (b*d)/(a*(-a + b + c)/e1) +
	  (-(a*c) + pow(c,2) + b*(c - d))/((a - b - c)*(b + c)/e2);
      pmat[MI(1,4,5)] = (b*(a - d)*(b + c - b/e1 - c/e1 +
				    a*(1/e1 - e2/e1)))/
	  (a*(a - b - c)*(b + c)/e1);
      pmat[MI(2,3,5)] = (d - d*e1)/a;
      pmat[MI(2,4,5)] = ((a - d)*(-1 + 1/e1))/(a/e1);
  }
  else if (all_equal(b+c,d+e) && !all_equal(a,d+e)){
      pmat[MI(0,1,5)] = (a*(-1 + e2/e1))/((a - b - c)/e1);
      pmat[MI(0,2,5)] = (a*b*(1/e2 + 1/e1*(-1 + a*t - b*t - c*t)))/
	  (pow(-a + b + c,2)/(e1*e2));
      pmat[MI(0,3,5)] = (pow(c,2) + b*(c + d))/pow(b + c,2) -
	  (-(a*c) + pow(c,2) + b*(c + d))/(pow(-a + b + c,2)/e1) +
	  (a*b*d)/((b + c)*pow(-a + b + c,2)/e2) -
	  (a*(pow(c,2) + pow(b,2)*d*t + b*(c + d + c*d*t)))/
	  ((a - b - c)*pow(b + c,2)/e2);
      pmat[MI(0,4,5)] = b*(b + c - d)*(pow(b + c,-2) - 1/(pow(-a + b + c,2)/e1) -
				       (a*(a*(1 + b*t + c*t) - (b + c)*(2 + b*t + c*t)))/
				       (pow(b + c,2)*pow(-a + b + c,2)/e2));
      pmat[MI(1,2,5)] = (b*t)*e2;
      pmat[MI(1,3,5)] = (pow(c,2)*(-1 + 1/e2) - pow(b,2)*d*t +
			 b*(d*(-1 + 1/e2) + c*(-1 + 1/e2 - d*t)))/
	  (pow(b + c,2)/e2);
      pmat[MI(1,4,5)] = (b*(b + c - d)*(-1 + 1/e2 - b*t - c*t))/
	  (pow(b + c,2)/e2);
      pmat[MI(2,3,5)] = (d - d*e2)/(b + c);
      pmat[MI(2,4,5)] = -(((b + c - d)*(-1 + e2))/(b + c));
  }
  else if (all_equal(a,b+c) && all_equal(a,d+e)){
      pmat[MI(0,1,5)] = (a*t)*e1;
      pmat[MI(0,2,5)] = (a*b*pow(t,2))/(2/e1);
      pmat[MI(0,3,5)] = (2*b*d*(-1 + 1/e1) - 2*pow(a,3)*t -
			 2*a*b*(-1 + 1/e1 + d*t) +
			 pow(a,2)*(-2 + 2/e1 + b*t*(2 - d*t)))/
	  (2*pow(a,2)/e1);
      pmat[MI(0,4,5)] = -(b*(a - d)*(2 - 2/e1 + 2*a*t + pow(a,2)*pow(t,2)))/
	  (2*pow(a,2)/e1);
      pmat[MI(1,2,5)] = (b*t)*e1;
      pmat[MI(1,3,5)] = (pow(a,2)*(-1 + 1/e1) + b*d*(-1 + 1/e1) -
			 a*b*(-1 + 1/e1 + d*t))/(pow(a,2)/e1);
      pmat[MI(1,4,5)] = (b*(a - d)*(-1 + 1/e1 - a*t))/(pow(a,2)/e1);
      pmat[MI(2,3,5)] = (d - d*e1)/a;
      pmat[MI(2,4,5)] = ((a - d)*(-1 + 1/e1))/(a/e1);
  }
  else {
      pmat[MI(0,1,5)] = (a*(-1 + e2/e1))/((a - b - c)/e1);
      pmat[MI(0,2,5)] = -((a*b*(d + e + a*e2/e1 - d*e2/e1 -
				e*e2/e1 - a*e3/e1 +
				b*(-1 + e3/e1) + c*(-1 + e3/e1)))/
			  ((-a + b + c)*(b + c - d - e)*(-a + d + e)/e1));
      pmat[MI(0,3,5)] = (b*d + c*(d + e))/((b + c)*(d + e)) +
	  (-(a*c) + b*d + c*(d + e))/((-a + b + c)*(a - d - e)/e1) -
	  (a*(b*(c - d) + c*(c - d - e)))/
	  ((a - b - c)*(b + c)*(b + c - d - e)/e2) -
	  (a*b*d)/((a - d - e)*(b + c - d - e)*(d + e)/e3);
      pmat[MI(0,4,5)] = b*e*(1/((b + c)*(d + e)) - 1/((-a + b + c)*(-a + d + e)/e1) +
			     a/((a - b - c)*(b + c)*(b + c - d - e)/e2) +
			     a/((a - d - e)*(d + e)*(-b - c + d + e)/e3));
      pmat[MI(1,2,5)] = (b*(-e2 + e3))/(b + c - d - e);
      pmat[MI(1,3,5)] = (b*d + c*(d + e))/((b + c)*(d + e)) +
	  (b*(-c + d) + c*(-c + d + e))/
	  ((b + c)*(b + c - d - e)/e2) -
	  (b*d)/((b + c - d - e)*(d + e)/e3);
      pmat[MI(1,4,5)] = -((b*e*(((d + e)*(-1 + 1/e2))*e2 +
				b*(-1 + e3) + c*(-1 + e3)))/
			  ((b + c)*(b + c - d - e)*(d + e)));
      pmat[MI(2,3,5)] = (d - d*e3)/(d + e);
      pmat[MI(2,4,5)] = (e - e*e3)/(d + e);
  }
}
