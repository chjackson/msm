
typedef double (*hmmfn)(double x, double *pars);
typedef void (*dhmmfn)(double x, double *pars, double *d);

double hmmCat(double x, double *pars);
double hmmIdent(double x, double *pars);
double hmmUnif(double x, double *pars);
double hmmNorm(double x, double *pars);
double hmmLNorm(double x, double *pars);
double hmmExp(double x, double *pars);
double hmmGamma(double x, double *pars);
double hmmWeibull(double x, double *pars);
double hmmPois(double x, double *pars);
double hmmBinom(double x, double *pars);
double hmmTNorm(double x, double *pars);
double hmmMETNorm(double x, double *pars);
double hmmMEUnif(double x, double *pars);
double hmmNBinom(double x, double *pars);
double hmmBeta(double x, double *pars);
double hmmT(double x, double *pars);

void DhmmCat(double x, double *pars, double *d);
void DhmmIdent(double x, double *pars, double *d);
void DhmmUnif(double x, double *pars, double *d);
void DhmmNorm(double x, double *pars, double *d);
void DhmmLNorm(double x, double *pars, double *d);
void DhmmExp(double x, double *pars, double *d);
void DhmmGamma(double x, double *pars, double *d);
void DhmmWeibull(double x, double *pars, double *d);
void DhmmPois(double x, double *pars, double *d);
void DhmmBinom(double x, double *pars, double *d);
void DhmmTNorm(double x, double *pars, double *d);
void DhmmMETNorm(double x, double *pars, double *d);
void DhmmMEUnif(double x, double *pars, double *d);
void DhmmNBinom(double x, double *pars, double *d);
void DhmmBeta(double x, double *pars, double *d);
void DhmmT(double x, double *pars, double *d);
