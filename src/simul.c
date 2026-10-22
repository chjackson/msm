/* Low-level engines for simulating from a continuous-time Markov
   multi-state process, with covariates that are piecewise-constant in
   time.  Used by sim.msm (exact transition times for one individual)
   and simmulti.msm (states at fixed observation times for many
   individuals).

   All data preparation (covariate linear predictors, censoring,
   misclassification, assembly of output data frames) is done in R.
   */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>

/* Transition intensity model, shared by all intervals and subjects.
   q      nst x nst matrix of baseline intensities (diagonal ignored)
   qindex nst x nst integer matrix: column of eta holding the linear
          predictor for transition r->s, or -1 if r->s is not allowed
   eta    nlp x nq matrix of linear predictors, X %*% beta.  nq may be 0
          if there are no covariates
   absorb integer indicator of absorbing states */

typedef struct {
    int nst;
    const double *q;
    const int *qindex;
    const double *eta;
    R_xlen_t nlp;
    int nq;
    const int *absorb;
    double *rate; /* workspace, length nst */
} simmodel;

/* Growable record of a simulated path: times of transitions and the
   states entered. */

typedef struct {
    int n;
    int cap;
    double *time;
    int *state;
} simpath;

static void path_init(simpath *p, int cap)
{
    p->n = 0;
    p->cap = (cap > 4) ? cap : 4;
    p->time = (double *) R_alloc(p->cap, sizeof(double)); // freed when .Call returns
    p->state = (int *) R_alloc(p->cap, sizeof(int));
}

static void path_add(simpath *p, double t, int s)
{
    if (p->n == p->cap) {
	int newcap = 2 * p->cap;
	double *nt = (double *) R_alloc(newcap, sizeof(double));
	int *ns = (int *) R_alloc(newcap, sizeof(int));
	memcpy(nt, p->time, p->n * sizeof(double));
	memcpy(ns, p->state, p->n * sizeof(int));
	p->time = nt; p->state = ns; p->cap = newcap;
    }
    p->time[p->n] = t;
    p->state[p->n] = s;
    p->n++;
}

/* Simulate one path from state s (0-based) at time tstart, through nint
   consecutive intervals.  Interval k ends at cut[k+1], or at tend for
   the last one, and its intensities are q * exp(eta[lprow0 + k, ]).
   Within an interval the intensities are fixed, so the sojourn is
   exponential, and memorylessness lets the clock restart at the end of
   the interval.  Transitions are appended to path, and the simulation
   stops on entering an absorbing state.  Returns the final state. */

static int sim_path(const simmodel *m, int s, double tstart,
		    const double *cut, int nint, double tend,
		    R_xlen_t lprow0, simpath *path)
{
    int k, r, col;
    double u = tstart, tk, w, tot, cum, target;
    for (k = 0; k < nint && !m->absorb[s]; k++) {
	tk = (k + 1 < nint) ? cut[k+1] : tend;
	for (;;) {
	    tot = 0;
	    for (r = 0; r < m->nst; r++) {
		col = m->qindex[s + m->nst*r];
		if (col < 0) m->rate[r] = 0;
		else {
		    m->rate[r] = m->q[s + m->nst*r];
		    if (m->nq > 0)
			m->rate[r] *= exp(m->eta[lprow0 + k + m->nlp*col]);
		}
		tot += m->rate[r];
	    }
	    if (tot <= 0) return s;
	    w = exp_rand() / tot;
	    if (u + w >= tk) break;
	    u += w;
	    target = unif_rand() * tot;
	    cum = 0;
	    for (r = 0; r < m->nst; r++) {
		if (m->rate[r] <= 0) continue;
		cum += m->rate[r];
		if (target < cum) break;
	    }
	    if (r == m->nst) /* guard against rounding */
		for (r = m->nst - 1; m->rate[r] <= 0; r--) ;
	    s = r;
	    path_add(path, u, s);
	    if (m->absorb[s]) return s;
	}
	u = tk;
    }
    return s;
}

static void model_init(simmodel *m, SEXP qbase, SEXP qindex, SEXP lp,
		       SEXP absorb)
{
    m->nst = nrows(qbase);
    m->q = REAL(qbase);
    m->qindex = INTEGER(qindex);
    m->eta = REAL(lp);
    m->nlp = nrows(lp);
    m->nq = (m->nlp > 0) ? ncols(lp) : 0;
    m->absorb = INTEGER(absorb);
    m->rate = (double *) R_alloc(m->nst, sizeof(double));
}

static SEXP named_list3(const char *n1, SEXP x1, const char *n2, SEXP x2,
			const char *n3, SEXP x3)
{
    int n = (x3 == R_NilValue) ? 2 : 3;
    SEXP ans = PROTECT(allocVector(VECSXP, n));
    SEXP nms = PROTECT(allocVector(STRSXP, n));
    SET_VECTOR_ELT(ans, 0, x1); SET_STRING_ELT(nms, 0, mkChar(n1));
    SET_VECTOR_ELT(ans, 1, x2); SET_STRING_ELT(nms, 1, mkChar(n2));
    if (n == 3) { SET_VECTOR_ELT(ans, 2, x3); SET_STRING_ELT(nms, 2, mkChar(n3)); }
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/* Engine for sim.msm: one path with exact transition times.
   cuts    times at which the covariates change.  The process starts at
           mintime, in the interval beginning at cuts[0], and the last
           interval ends at maxtime.  Row k of lp gives the linear
           predictors in the interval starting at cuts[k]
   start   starting state (1-based)
   Returns a list with the times of each transition and the (1-based)
   states entered. */

SEXP msm_sim_path(SEXP qbase, SEXP qindex, SEXP lp, SEXP cuts, SEXP start,
		  SEXP mintime, SEXP maxtime, SEXP absorb)
{
    simmodel m;
    simpath path;
    int j;
    model_init(&m, qbase, qindex, lp, absorb);
    path_init(&path, 16);
    GetRNGstate();
    sim_path(&m, asInteger(start) - 1, asReal(mintime), REAL(cuts),
	     LENGTH(cuts), asReal(maxtime), 0, &path);
    PutRNGstate();
    SEXP rtime = PROTECT(allocVector(REALSXP, path.n));
    SEXP rstate = PROTECT(allocVector(INTSXP, path.n));
    for (j = 0; j < path.n; j++) {
	REAL(rtime)[j] = path.time[j];
	INTEGER(rstate)[j] = path.state[j] + 1;
    }
    SEXP ans = named_list3("time", rtime, "state", rstate, NULL, R_NilValue);
    UNPROTECT(2);
    return ans;
}

/* Engine for simmulti.msm: states at fixed observation times for many
   subjects.  Indices into data rows are 0-based here, 1-based in the
   output.
   lp       nobs x nq matrix of linear predictors at each observation.
            Row j applies between observation j and the next
   times    observation times, sorted within subject
   firstobs integer vector of length npts+1: rows of subject i are
            firstobs[i], ..., firstobs[i+1]-1
   start    starting state of each subject (1-based)
   death    integer indicator of "death" states (exact entry time known)
   dropabs  if 1, retain only the first observation in an absorbing state

   Returns a list with components
   row    data row each output observation corresponds to (1-based)
   time   observation time (the exact entry time for a death state)
   state  simulated true state (1-based)
*/

SEXP msm_simmulti_engine(SEXP qbase, SEXP qindex, SEXP lp, SEXP times,
			 SEXP firstobs, SEXP start, SEXP absorb, SEXP death,
			 SEXP dropabs)
{
    simmodel m;
    simpath path;
    int nobs = LENGTH(times);
    int npts = LENGTH(firstobs) - 1;
    int drop = asInteger(dropabs);
    double *t = REAL(times);
    int *fo = INTEGER(firstobs), *st0 = INTEGER(start);
    int *isdeath = INTEGER(death);
    int i, j, jdeath, s, cur, idx, absorbed, nout = 0;
    double abstime;

    model_init(&m, qbase, qindex, lp, absorb);
    path_init(&path, 16);

    /* Each subject contributes at most one extra row (a death time) */
    int maxout = nobs + npts;
    int *orow = (int *) R_alloc(maxout, sizeof(int));
    int *ostate = (int *) R_alloc(maxout, sizeof(int));
    double *otime = (double *) R_alloc(maxout, sizeof(double));
    int *pstate = (int *) R_alloc(nobs > 0 ? nobs : 1, sizeof(int));

    GetRNGstate();
    for (i = 0; i < npts; i++) {
	int f = fo[i], l = fo[i+1];
	if (l <= f) continue;
	if ((i % 1000) == 0) R_CheckUserInterrupt();
	path.n = 0;
	s = sim_path(&m, st0[i] - 1, t[f], &t[f], l - 1 - f, t[l-1], f, &path);
	absorbed = m.absorb[s];
	abstime = (path.n > 0) ? path.time[path.n - 1] : t[f];

	/* State at each observation time */
	cur = st0[i] - 1; idx = 0;
	for (j = f; j < l; j++) {
	    while (idx < path.n && path.time[idx] <= t[j])
		cur = path.state[idx++];
	    pstate[j] = cur;
	}

	if (absorbed && isdeath[s]) {
	    /* Keep observations before death, then one row at the exact
	       death time, which takes the place of the next scheduled
	       observation. */
	    for (jdeath = f; jdeath < l && t[jdeath] < abstime; jdeath++) {
		orow[nout] = jdeath + 1; otime[nout] = t[jdeath];
		ostate[nout] = pstate[jdeath] + 1; nout++;
	    }
	    if (jdeath == l) jdeath = l - 1;
	    orow[nout] = jdeath + 1; otime[nout] = abstime;
	    ostate[nout] = s + 1; nout++;
	    if (!drop)
		for (j = jdeath + 1; j < l; j++) {
		    orow[nout] = j + 1; otime[nout] = t[j];
		    ostate[nout] = s + 1; nout++;
		}
	} else {
	    for (j = f; j < l; j++) {
		/* keep the first observation at or after absorption */
		if (drop && absorbed && j > f && t[j-1] > abstime) break;
		orow[nout] = j + 1; otime[nout] = t[j];
		ostate[nout] = pstate[j] + 1; nout++;
	    }
	}
    }
    PutRNGstate();

    SEXP rrow = PROTECT(allocVector(INTSXP, nout));
    SEXP rtime = PROTECT(allocVector(REALSXP, nout));
    SEXP rstate = PROTECT(allocVector(INTSXP, nout));
    for (j = 0; j < nout; j++) {
	INTEGER(rrow)[j] = orow[j];
	REAL(rtime)[j] = otime[j];
	INTEGER(rstate)[j] = ostate[j];
    }
    SEXP ans = named_list3("row", rrow, "time", rtime, "state", rstate);
    UNPROTECT(3);
    return ans;
}
