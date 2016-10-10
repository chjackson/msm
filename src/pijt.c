/* *****************************************************************

   pijt.c: Linear algebra routines for msm. To obtain Markov
   transition probabilities (the P-matrix) from transition intensities
   (the Q-matrix).

   Copyright (c) Chris Jackson 2001-2005

   Pade approximants method of calculating matrix exponentials is
   based on matexp.cc from JAGS, Copyright (c) 2004, Martyn Plummer.
   http://www-fis.iarc.fr/~martyn/software/jags

   ******************************************************************  */

#include "msm.h"
/* Crude benchmarks have shown that using the eigensystem routines
   from LINPACK, and the matrix inversion routines from LAPACK, is the
   faster combination. However...
   Occasional crashes on extremely flat optimisations.
   Valgrind reveals memory errors within LINPACK Fortran eigen code.
   Crashes don't occur if use LAPACK code.
*/
#define _USE_LAPACK_EIGEN_
#define _USE_LAPACK_INVERSE_
#define MEXP_PADE 1
#define MEXP_SERIES 2
#include "R_ext/Lapack.h"
#include "R_ext/Rdynload.h"

#define NODERIVDEBUG

/* Interface to expm package. */
typedef enum {Ward_2, Ward_1, Ward_buggy_octave} precond_type;
void (*expm)(double *x, int n, double *z, precond_type precond_kind);
void R_init_msm(DllInfo *dll)
{
    expm = (void (*) (double*, int, double*, precond_type)) R_GetCCallable("expm", "expm");
}

/* Set A to be an n x n identity matrix */

void FormIdentity(Matrix A, int n)
{
    int i;
    for (i = 0; i < (n*n); ++i)
	A[i] = 0;
    for (i = 0; i < n; ++i)
	A[MI(i, i, n)] = 1;
}

/* Invert a matrix by calling LINPACK QR decomposition routines */
/* DGETRF/DGETRI method (current)  is from LAPACK. This is the fastest.   */
/* DQR method (used in <=0.4.1) is from LINPACK */

void MatInvDGE(Matrix A, Matrix Ainv, int n)
{
    int i, j;
    Matrix temp = (Matrix) Calloc(n*n, double);
    Matrix work = (Matrix) Calloc(n*n, double);
    int nsq=n*n, info, *pivot = Calloc(n, int);
    for (i=0; i<nsq; ++i)
      temp[i] = A[i];
    F77_CALL(dgetrf) (&n, &n, temp, &n, pivot, &info);
    if (info < 0)
	REprintf("error code %d from Lapack routine dgetrf\n", info);
    F77_CALL(dgetri) (&n, temp, &n, pivot, work, &nsq, &info);
    if (info < 0)
	REprintf("error code %d from Lapack routine dgetri\n", info);
    for (i=0; i<n; ++i)
	for (j=0; j<n; ++j)
	    Ainv[MI(i,j,n)] = temp[MI(i,j,n)];
    Free(work); Free(pivot); Free(temp);
}

void MatInvDQR(Matrix A, Matrix Ainv, int n)
{
    int i, rank;
    Matrix temp = (Matrix) Calloc(n*n, double);
    Matrix work = (Matrix) Calloc(n*n, double);
    Matrix qraux = (Matrix) Calloc(n*n, double);
    Matrix ident = (Matrix) Calloc(n*n, double);
    int nsq=n*n, info, *pivot = Calloc(n, int);
    double tol=1e-07;
    for (i=0; i<nsq; ++i)
      temp[i] = A[i];
    F77_CALL(dqrdc2) (temp, &n, &n, &n, &tol, &rank, qraux, pivot, work);
    FormIdentity(ident, n);
    F77_CALL(dqrcf) (temp, &n, &rank, qraux, ident, &n, Ainv, &info);
    if (info < 0)
	REprintf("error code %d from Linpack routine dqrcf\n", info);
    Free(temp); Free(work); Free(qraux); Free(ident);Free(pivot);
}

void MatInv(Matrix A, Matrix Ainv, int n)
{
#ifdef _USE_LAPACK_INVERSE_
	MatInvDGE(A, Ainv, n);
#endif
#ifdef _USE_LINPACK_INVERSE_
	MatInvDQR(A, Ainv, n);
#endif
}

/* Multiply two matrices together */

void MultMat(Matrix A, Matrix B, int arows, int acols, int bcols, Matrix AB)
{
    int i, j, k;
    for (i = 0; i < arows; i++) {
	for (j = 0; j < bcols; j++) {
	    AB[MI(i, j, bcols)] = 0;
	    for (k = 0; k < acols; k++)
		AB[MI(i, j, bcols)] += A[MI(i, k, acols)] * B[MI(k, j, bcols)];
	}
    }
}

/* Copy one matrix A to another B */

void CopyMat(Matrix A, Matrix B, int arows, int acols)
{
    int i;
    for (i = 0; i < arows*acols; i++)
	B[i] = A[i];
}

/* Pre-multiply a general matrix by a diagonal matrix (given by a vector) */

void MultMatDiag(vector diag, Matrix B, int n, Matrix AB)
{
    int i, j;
    for (i = 0; i < (n*n); ++i)
	AB[i] = 0;
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    AB[MI(i, j, n)] += diag[i] * B[MI(i, j, n)];
	}
    }
}

/* Calculate a matrix exponential using a power series approximation
   Adapted from mexp in Jim Lindsey's "rmutil" package
   exp(A)  =  I + A + AA/2  + AAA/2*3 + AAAA/2*3*4 + ...
   Overflow correction:
   exp(A)  =  exp(A/2^3) exp(A/2^3) exp(A/2^3)
*/

void MatrixExpSeries(Matrix A, int n, Matrix expmat, double t)
{
    int i, j;
    int order = 20;   /* number of terms in series */
    int overflow_correct = 3;
    Matrix Apower = Calloc(n*n, double);
    Matrix Temp = Calloc(n*n, double);
    Matrix AA = Calloc(n*n, double);
    for (i=0; i<(n*n); ++i)
	AA[i] = A[i] * (t / pow(2, overflow_correct));
    FormIdentity(expmat, n);
    FormIdentity(Apower, n);
    for (i=1; i<=order; i++) {
	MultMat(AA, Apower, n, n, n, Temp);
	for (j=0; j<(n*n); ++j){
	    Apower[j] = Temp[j] / i;
	    expmat[j] += Apower[j];
	}
    }
    for (i=0; i<overflow_correct; ++i){
	MultMat(expmat, expmat, n, n, n, Temp);
	CopyMat(Temp, expmat, n, n);
    }
    Free(Apower); Free(Temp); Free(AA);
}

/* Pade approximants method of calculating matrix exponentials
   From matexp.cc from JAGS by Martyn Plummer  Copyright (c) 2004
   http://www-fis.iarc.fr/~martyn/software/jags
*/

static void solve(double *X, double const *A, double const *B, int n)
{
    /* Solve AX = B, where all matrices are square */
    int N = n*n;
    double *Acopy = Calloc(N, double);
    double *work = Calloc(N, double);
    int *ipiv = Calloc(N, int);
    int info = 0;
    static int c_1 = 1;
    F77_CALL(dcopy)(&N, A, &c_1, Acopy, &c_1);
    F77_CALL(dcopy)(&N, B, &c_1, X, &c_1);
    F77_CALL(dgesv)(&n, &n, Acopy, &n, ipiv, X, &n, &info);
    if (info < 0)
	REprintf("argument %d of Lapack routine dgesv had illegal value\n", -info);
    if (info > 0)
	REprintf("Lapack routine dgesv: system is exactly singular\n");
    Free(Acopy);
    Free(ipiv);
    Free(work);
}

static void
padeseries (double *Sum, double *A, int m, int order,
            double scale, double *Temp)
{
    int i, j, r;
    int N = m*m;
    FormIdentity(Sum, m);
    for (j = order; j >= 1; --j) {
	double s = (order-j+1) / (j*(2*order-j+1) * scale);
	MultMat(Sum, A, m, m, m, Temp);
	for (i = 0; i < N; ++i) {
	    Sum[i] = Temp[i] * s;
	}
	for (r = 0; r < m; ++r) {
	    Sum[r*m+r] += 1;
	}
    }
}

void
MatrixExpPade(double *ExpAt, double *A, int n, double t)
{
  /* Calculate exp(A*t) by diagonal Pade approximation with scaling and
     squaring */
    int i, j;
  int order = 8;
  int N = n*n;
  double *workspace =  Calloc( 4*N, double);
  double * Temp = workspace;
  double * At = workspace + N;
  double * Num = workspace + 2*N;
  double * Denom = workspace + 3*N;
  double l1 = F77_CALL(dlange)("1", &n, &n, At, &n, 0); /* L-1 norm */
  double linf = F77_CALL(dlange)("i", &n, &n, At, &n, Temp); /* L-Infinity norm */
  double K = (log(l1) + log(linf))/log(4);
  int npower = (R_FINITE(K) ? (int)(K)+4 : NA_INTEGER);
  double scale = 1;

  /* Multiply by t */

  for (i = 0; i < N; ++i) {
    At[i] = A[i] * t;
  }

  /* Scale the matrix by a power of 2 */
  /*
     The expression below is not clear because it is optimized.  The
     idea is that sqrt(l1 * linf) is an upper bound on the L2 norm of
     the matrix At (i.e the largest eigenvalue). We want to take the
     log, to base 2 of this to get the smallest K, st ||At/2^K|| <= 1.
  */
  if (npower < 0) {
      npower = 0;
  }
  for (i = 0; i < npower; ++i) {
    scale *= 2;
  }

  /* Calculate exp(A/scale) by Pade series  */

  padeseries (Num, At, n, order, scale, Temp);
  for (i = 0; i < N; ++i) {
    At[i] = -At[i];
  }
  padeseries (Denom, At, n, order, scale, Temp);
  solve(ExpAt, Denom, Num, n);

  /* Now repeatedly square the result */
  for (i = 0; i < npower; ++i) {
    for (j = 0; j < N; ++j) {
      Temp[j] = ExpAt[j];
    }
    MultMat(Temp, Temp, n, n, n, ExpAt);
  }
  Free(workspace);
}

/* Tests if a vector has any non-unique entries */

int repeated_entries(vector vec, int n)
{
    int i, j;
    for (i=1; i<n; ++i)
	for (j=0; j<i; ++j)
	    if (vec[j] == vec[i])
		return 1;
    return 0;
}

void Eigen(Matrix mat, int n, vector revals, vector ievals, Matrix evecs, int *err)
{
    int i, nsq=n*n;
#ifdef _USE_LINPACK_EIGEN_
    int matz = 1;
#endif
#ifdef _USE_LAPACK_EIGEN_
    int lwork = -1;
    char jobVL[1], jobVR[1];
    double *left=0, tmp;
#endif
    Matrix work = (Matrix) Calloc(nsq, double);
    iMatrix worki = (iMatrix) Calloc(nsq, int);
    Matrix temp = (Matrix) Calloc(nsq, double);
    for (i=0; i<nsq; ++i) {
	/* Check whether any of the elements of Q have overflowed.  If
	   so, Fortran eigen function will hang in a infinite loop, so
	   bail out before this happens.  */
	if (!R_FINITE(mat[i])) { 
	    error("numerical overflow in calculating likelihood\n");
	}
	temp[i] = mat[i];
    }
#ifdef _USE_LINPACK_EIGEN_
    F77_CALL(rg) (&n, &n, temp, revals, ievals, &matz, evecs, worki, work, err);
#endif
#ifdef _USE_LAPACK_EIGEN_
    jobVL[0] = 'N'; jobVR[0] = 'V';
    /* calculate optimal size of workspace */
    F77_CALL(dgeev)(jobVL, jobVR, &n, temp, &n, revals, ievals, left, &n, evecs, &n, &tmp, &lwork, err);
    lwork = (int) tmp;
    work = (Matrix) Realloc(work, lwork, double);
    /* calculate eigensystem */
    F77_CALL(dgeev)(jobVL, jobVR, &n, temp, &n, revals, ievals, left, &n, evecs, &n, work, &lwork, err);
#endif
    Free(work); Free(worki); Free(temp);
}

/* Compute exponential of a general matrix */
/* First try to use eigensystem decomposition */
/* If matrix has repeated eigenvalues, then use Pade approximation
   with scaling and squaring, or (less robust) power series.  */

/* TODO: Would it be more economical in the long run to test for
   invertibility by calculating determinant, thus reducing number of
   cases where it is necessary to use power series or Pade? */

void MatrixExpMSM(Matrix mat, int n, Matrix expmat, double t,
	       int degen, /* currently not used */
	       int method)
{
    int i, err=0, complex_evals=0, nsq=n*n;
    Matrix work = (Matrix) Calloc(nsq, double);
    vector revals = (vector) Calloc(n, double);
    vector ievals = (vector) Calloc(n, double);
    Matrix evecs = (Matrix) Calloc(nsq, double);
    Matrix evecsinv = (Matrix) Calloc(nsq, double);
    /* calculate eigensystem */
    if (!degen)
	Eigen(mat, n, revals, ievals, evecs, &err);
    /* Check for eigenvalues with nonzero imaginary part */
    for (i=0; i<n; ++i)
      if (!all_equal(ievals[i],0)){
	complex_evals = 1;
	break;
      }
    if (repeated_entries (revals, n) || (err != 0) || degen || complex_evals){
	if (method == MEXP_SERIES)
	    MatrixExpSeries(mat, n, expmat, t);
	else
	    MatrixExpPade(expmat, mat, n, t);
    }
    else {
	for (i=0; i<n; ++i)
	    revals[i] = exp(revals[i] * t);
	MatInv(evecs, evecsinv, n);
	MultMatDiag(revals, evecsinv, n, work);
	MultMat(evecs, work, n, n, n, expmat);
    }
    Free(work);  Free(revals);  Free(ievals); Free(evecs);  Free(evecsinv);
}

/* Exponential of a matrix.  If matrix represents one of certain
   Markov model structures, then use appropriate analytic formulae.
   Interface to R MatrixExp function, and also used in Pmat */

void MatrixExpR(double *mat, int *n, double *expmat, double *t,
		int *method, int *iso, int *perm, int *qperm,
		int *degen){
    if (*iso > 0)
	AnalyticP(expmat, *t, *n, *iso, perm, qperm, mat, degen);
    else
	MatrixExpMSM(mat, *n, expmat, *t, *degen, *method);
}

void MatrixExpEXPM(double *mat, int *n, double *expmat, double *t,
		int *method, int *iso, int *perm, int *qperm,
		   int *degen, int *err){
    int i;
    int nsq = (*n)*(*n);
    double *matt = Calloc(nsq, double);
    if (*iso > 0)
	AnalyticP(expmat, *t, *n, *iso, perm, qperm, mat, degen);
    else {
	for (i=0; i<((*n)*(*n)); ++i) {
	    matt[i] = (*t) * mat[i];
	    /* Check whether any of the elements of Q have overflowed.  If
	       so, Fortran eigen function will hang in a infinite loop, so
	       bail out before this happens.  */
	/* Could we return loglik = -Inf instead of halting by setting
	   error code?  return all zeros for pmat in that case?
	   Doesn't help convergence with test case in test/rory.r: lik
	   stays at zero.
	 */
	    if (!R_FINITE(matt[i])){
//		*err = -1; return;
		error("numerical overflow in calculating likelihood\n");
	    }
	}
	expm(matt, *n, expmat, Ward_2);
    }
    Free(matt);
}

/* Returns i-j transition intensity time t given vectors of intensities and transition indicators */

/* Calculates the whole transition matrix in time t given an intensity matrix */

void Pmat(Matrix pmat, double t, Matrix qmat, int nstates, int exacttimes, int iso, ivector perm, ivector qperm, int use_expm)
{
    int i,j,method=MEXP_PADE,degen=0,err=0;
    double pii;
    if (exacttimes) {
	for (i=0; i<nstates; ++i) {
	  pii = exp(t * qmat[MI(i, i, nstates)] );
	  for (j=0; j<nstates; ++j) {
	    pmat[MI(i, j, nstates)] = ( i==j  ?  pii  : pii * qmat[MI(i, j, nstates)] );
	  }
	}
    }
    else {
	if (use_expm)
	    MatrixExpEXPM(qmat, &nstates, pmat, &t, &method, &iso, perm, qperm, &degen, &err);
	else
	    MatrixExpR(qmat, &nstates, pmat, &t, &method, &iso, perm, qperm, &degen);
	/* Floating point fuzz sometimes causes trouble */
	for (i=0; i<nstates; ++i)
	    for (j=0; j<nstates; ++j) {
		/* if (err==-1) pmat[MI(i, j, nstates)] = 0;  */ /* leads to zero likelihood */
		if (pmat[MI(i, j, nstates)] < DBL_EPSILON) pmat[MI(i, j, nstates)] = 0;
		if (pmat[MI(i, j, nstates)] > 1 - DBL_EPSILON) pmat[MI(i, j, nstates)] = 1;
	    }
    }
}


double pijdeath(int r, int s, Matrix pmat, Matrix qmat, int n)
{
    int j;
    double contrib;
    if (r == s) return 1;  /* absorbing-same absorbing transition has probability 1 */
    else {    /* sum over unobserved state at the previous instant */
	contrib = 0;
	for (j = 0; j < n; ++j)
	    if (j != s) {
		contrib += pmat[MI(r, j, n)] * qmat[MI(j,s,n)];
	    }
    }
    return contrib;
}




/***************************************************

 CODE FOR DERIVATIVES OF P MATRIX.

***************************************************/

/*

qij exp (qii t)
dqij exp(qii t) + dqii qij t exp(qii t)
exp(qii t) ( dqij + dqii qij t )

or exp(qii t) if diag

*/

void DPmatEXACT(Array3 dqmat, Matrix qmat, int n, int npars, Array3 dpmat, double t)
{
    int i,j,p;
    for (i=0; i<n; ++i) { /* rows */
	for (j=0; j<n; ++j) { /* columns */
	    for (p=0; p<npars; ++p) {
		if (i==j)
		    dpmat[MI3(i,j,p,n,n)] = dqmat[MI3(i,i,p,n,n)] * t * exp(qmat[MI(i,i,n)]*t);
		else
		    dpmat[MI3(i,j,p,n,n)] = exp(qmat[MI(i,i,n)]*t) *
			( dqmat[MI3(i,j,p,n,n)] + dqmat[MI3(i,i,p,n,n)] * qmat[MI(i,j,n)] * t);
	    }
	}
    }
}

void DMatrixExpSeries(Matrix DA, Matrix A, int n, int npars, Array3 DexpA, double t)
{
    int i, j, k, p;
    int order = 20;
    int nsq = n*n;
    double *DAp;
    vector tpower = Calloc(order+1, double); /* values of t^i / i! */
    Matrix DApower = Calloc(nsq, double); /* cumulative outer sum */
    Array3 Apower = Calloc(nsq*(order+1), double);  /* storage for successive powers of A  */
    Matrix Temp = Calloc(nsq, double); /* workspace  */
    Matrix Inner = Calloc(nsq, double); /* one term of inner sum */
    Matrix CInner = Calloc(nsq, double); /* cumulative inner sum */

    FormIdentity(&Apower[0], n);  /* A to the power 0 is identity */
    tpower[0] = 1;
    for (i=1; i<=order; ++i) {
	/* calculate ith power of A and store it in next nsq entries of Apower.  */
	MultMat(A, &Apower[(i-1)*nsq], n, n, n, &Apower[i*nsq]);
	tpower[i] = tpower[i-1] * t/i;
    }
    for (p=0; p<npars; ++p) {
	DAp = &(DA[MI3(0,0,p,n,n)]);
	for (k=0; k<nsq; ++k)
	    DexpA[p*nsq + k] =  DAp[k]*tpower[1];
	for (i=2; i<=order; ++i) {
	    for (k=0; k<nsq; ++k)
		CInner[k] = 0;
	    for (j=0; j <=(i-1); ++j) {
		MultMat(&Apower[j*nsq], DAp, n, n, n, Temp);
		MultMat(Temp, &Apower[(i-1-j)*nsq], n, n, n, Inner);
		for (k=0; k<nsq; ++k)
		    CInner[k] += Inner[k];
	    }
	    for (k=0; k<nsq; ++k)
		DexpA[p*nsq + k] += CInner[k]*tpower[i]; /* Fill in next nsq entries of DexpA with derivative WRT pth parameter.  */
	}
    }
    Free(tpower); Free(DApower); Free(Apower); Free(Temp); Free(Inner); Free(CInner);
}

void dpijdeath(int r, int s, Array3 dpmat, Matrix pmat, Array3 dqmat, Matrix qmat, int n, int npars, vector dcontrib)
{
    int k, p=0;
    for (p=0; p<npars; ++p) {  /* p indexes distinct pars */
	dcontrib[p] = 0;
	for (k=0; k<n; ++k)
	    if (k != s)
		dcontrib[p] += dpmat[MI3(r,k,p,n,n)] * qmat[MI(k,s,n)] + pmat[MI(r,k,n)] * dqmat[MI3(k,s,p,n,n)];
    }
}

/* Derivatives of P matrix (Kalbfleisch & Lawless) */

void DPmat(Array3 dpmat, double t, Array3 dqmat, Matrix qmat, int n, int npars, int exacttimes)
{
    int i, j, p, err=0;
    double eit, ejt;
    Matrix DQ;
    vector revals = (vector) Calloc(n, double);
    vector ievals = (vector) Calloc(n, double);
    Matrix evecs = (Matrix) Calloc(n*n, double);
    Matrix evecsinv = (Matrix) Calloc(n*n, double);
    Matrix work = (Matrix) Calloc(n*n, double);
    Matrix G = (Matrix) Calloc(n*n, double);
    Matrix V = (Matrix) Calloc(n*n, double);

    if (exacttimes) {
	DPmatEXACT(dqmat, qmat, n, npars, dpmat, t);
    }
    else {
	Eigen(qmat, n, revals, ievals, evecs, &err);
	if (err > 0)
	    REprintf("error code %d from EISPACK eigensystem routine rg\n", err);
	if (repeated_entries (revals, n)) {
	    DMatrixExpSeries(dqmat, qmat, n, npars, dpmat, t);
	}
	else {
	    MatInv(evecs, evecsinv, n);
	    for (p=0; p<npars; ++p) {
		DQ = &(dqmat[MI3(0,0,p,n,n)]);
		MultMat(DQ, evecs, n, n, n, work);
		MultMat(evecsinv, work, n, n, n, G);
		for (i=0; i<n; ++i) {
		    eit = exp(revals[i] * t);
		    for (j=0; j<n; ++j) {
			if (i==j)
			    V[MI(i,j,n)] = G[MI(i,i,n)] * t * eit;
			else {
			    ejt = exp(revals[j] * t);
			    V[MI(i,j,n)] = G[MI(i,j,n)] * (eit - ejt) / (revals[i] - revals[j]);
			}
		    }
		}
		MultMat(V, evecsinv, n, n, n, work);
		MultMat(evecs, work, n, n, n, &(dpmat[MI3(0, 0, p, n, n)]));
#ifdef DERIVDEBUG
		for (i=0; i<n; ++i) {
		    for (j=0; j<n; ++j) {
			printf("DQ[%d,%d]=%lf,",i,j,DQ[MI(i,j,n)]);
		    }
		    printf("\n");
		}
		printf("\n");
		for (i=0; i<n; ++i) {
		    for (j=0; j<n; ++j) {
			printf("dpmat[%d,%d]=%lf,",i,j,dpmat[MI3(i, j, p, n, n)]);
		    }
		    printf("\n");
		}
#endif
	    }
	}
    }
    Free(revals); Free(ievals); Free(evecs); Free(evecsinv); Free(work); Free(G); Free(V);
}
