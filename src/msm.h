#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Applic.h>

/* index to treat a vector as a matrix. ith row, jth column. Fills columns first, as in R */
#define MI(i, j, nrows) ( (int) ((j)*(nrows) + (i)) )
/* index to treat a vector as a 3-dimensional array. Left-most index varies fastest, as in R */
#define MI3(i, j, k, n1, n2) ( (int) ((k)*((n1)*(n2)) + (j)*(n1) + (i)) )
#define MI4(i, j, k, m, n1, n2, n3) ( (int) ((m)*((n1)*(n2)*(n3)) + (k)*((n1)*(n2)) + (j)*(n1) + (i)) )

/* Macros to switch quickly between C and S memory handling. Currently not used */

#define USE_CALLOC
/* #define USE_SALLOC */

#ifdef USE_CALLOC
#define MSM_ALLOC(length, type) Calloc((length), type)
#define MSM_FREE(var) Free((var))
#else
#define MSM_ALLOC(length, type) (type *) S_alloc((length), sizeof(type))
#define MSM_FREE(var)
#endif

typedef double * Array3;
typedef double * Array4;
typedef double * Matrix;
typedef int * iMatrix;
typedef double * vector;
typedef int * ivector;

struct msmdata {
    /* for non-hidden model */
    int *fromstate;
    int *tostate;
    double *timelag;
    int *nocc;
    int *noccsum;
    int *whicha;
    int *obstypea;

    /* for hidden model */
    int *subject;
    double *time;
    double *obs; /* observed state or any other HMM observation */
    int *obstype;
    int *obstrue;
    int *pcomb;
    int *firstobs;

    int nagg;
    int n;
    int npts;
    int ntrans;
    int npcombs;
    int nout;
};

struct qmodel {
    int nst;
    int npars;
    int nopt;
    double *intens;
    double *dintens;
    int iso;
    int *perm;
    int *qperm;
    int expm;
    int nliks;
};

struct cmodel {
    int ncens;
    int *censor;
    int *states;
    int *index;
};

struct hmodel {
    int hidden;
    int mv;
    int ematrix;
    int *models;
    int totpars;
    int *npars;
    int *firstpar;
    double *pars;
    double *dpars;
    int nopt;
    double *initp;
};

typedef struct msmdata msmdata;
typedef struct qmodel qmodel;
typedef struct cmodel cmodel;
typedef struct hmodel hmodel;

int repeated_entries(vector vec, int n);
double logit(double x);
double expit(double x);
double identity(double x);
int all_equal(double x, double y);
void MatrixExpPadeR(double *ExpAt, double *A, int *n, double *t);
void AnalyticP(Matrix pmat, double t, int nstates, int iso, int *perm, int *qperm, Matrix qmat, int *degen);
double pijdeath(int r, int s, Matrix pmat, Matrix qmat, int n);
void Pmat(Matrix pmat, double t, Matrix qmat, int nstates, int exacttimes, int iso, int *perm, int *qperm, int expm);
void DPmat(Array3 dpmat, double t, Array3 dqmat, Matrix qmat, int n, int np, int exacttimes);
void dpijdeath(int r, int s, Array3 dpmat, Matrix pmat, Array4 dqmat, Matrix qmat, int n, int npars, Matrix dcontrib);
