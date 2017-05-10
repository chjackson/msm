#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "expm.h"

/* .C calls */
extern void MatrixExpR(void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP msmCEntry(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"MatrixExpR", (DL_FUNC) &MatrixExpR, 9},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"msmCEntry", (DL_FUNC) &msmCEntry, 8},
    {NULL, NULL, 0}
};

void R_init_msm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    expm = (void (*) (double*, int, double*, precond_type)) R_GetCCallable("expm", "expm");
}
