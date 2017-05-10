/* Interface to expm package. */
typedef enum {Ward_2, Ward_1, Ward_buggy_octave} precond_type;
void (*expm)(double *x, int n, double *z, precond_type precond_kind);
