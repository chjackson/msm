### PACKAGE GLOBAL CONSTANTS ###

### List of allowed hidden Markov model distributions
### and names of parameters for each distribution
### MUST BE KEPT IN THE SAME ORDER as the C variable HMODELS in src/lik.c
.msm.HMODELPARS <- list(
                        categorical=c("ncats","basecat","p"),
                        identity = NULL,
                        uniform = c("lower", "upper"),
                        normal = c("mean", "sd"),
                        lognormal = c("meanlog", "sdlog"),
                        exponential = c("rate"),
                        gamma = c("shape","rate"),
                        weibull = c("shape","scale"),
                        poisson = c("rate"),
                        binomial = c("size","prob"),
                        truncnorm = c("mean", "sd", "lower", "upper"),
                        metruncnorm = c("mean", "sd", "lower", "upper", "sderr", "meanerr"),
                        meuniform = c("lower", "upper", "sderr", "meanerr"),
                        nbinom = c("disp","prob"),
                        beta = c("shape1","shape2"),
                        t = c("mean","scale","df")
                        )

## TODO - e.g. non-central beta, cauchy, chisq, noncentral chisq, F,
## non-central F, geometric, hypergeometric, logistic, noncentral t.

.msm.HMODELS <- names(.msm.HMODELPARS)

### Models with analytic derivatives available
.msm.HMODELS.DERIV <- c("categorical","identity","uniform","normal","lognormal","exponential",
                        "gamma","weibull","poisson","binomial","nbinom","beta","t")

### Models with expected information matrix available
.msm.HMODELS.INFO <- c("categorical","identity")

### Parameter in each distribution that can have covariates on it
.msm.LOCPARS <- c(categorical="p", identity=NA, uniform=NA, normal="mean", lognormal="meanlog",
                  exponential="rate", gamma="rate", weibull="scale",
                  poisson="rate", binomial="prob", truncnorm="mean",
                  metruncnorm="meanerr", meuniform="meanerr", nbinom="prob", beta=NA, t="mean")

### Link functions for generalised regressions.
### MUST BE KEPT IN SAME ORDER as LINKFNS in lik.c
.msm.LINKFNS <- c("identity", "log", "qlogis")
.msm.INVLINK <- c(identity="identity", log="exp", qlogis="plogis")

### Parameters which are always fixed, never estimated
.msm.AUXPARS <- c("lower", "upper", "which", "size", "meanerr", "ncats", "basecat", "p0", "pbase", "initpbase", "initp0")

### Parameters which should be defined as integer
.msm.INTEGERPARS <- c("size")

### Defined ranges for parameters
.msm.PARRANGES <- list(qbase=c(0, Inf), lower=c(-Inf,Inf), upper=c(-Inf, Inf),
                       mean=c(-Inf, Inf), sd=c(0, Inf),
                       meanlog=c(-Inf,Inf), sdlog=c(0, Inf), rate=c(0, Inf), shape=c(0, Inf), scale=c(0, Inf), shape1=c(0,Inf), shape2=c(0,Inf),
                       prob=c(0, 1), meanerr=c(-Inf, Inf), sderr=c(0, Inf), disp=c(0, Inf),
                       p=c(-Inf,Inf), # handled separately using multinomial logit
                       initp=c(-Inf,Inf), # handled separately using multinomial logit
                       df=c(0, Inf),
                       qcov=c(-Inf,Inf),hcov=c(-Inf,Inf),initpcov=c(-Inf,Inf)
                       )
for (i in .msm.AUXPARS) .msm.PARRANGES[[i]] <- c(-Inf, Inf)
.msm.PARRANGES <- do.call("rbind",.msm.PARRANGES)
colnames(.msm.PARRANGES) <- c("lower","upper")

### Transforms to optimise some parameters on a different scale
### Univariate transforms only: doesn't include multinomial logit transform used for misclassification p and initial state probs

.msm.TRANSFORMS <-
    do.call("rbind",
          apply(.msm.PARRANGES, 1,
                 function(x) {
                     if (identical(x, c(lower=0, upper=Inf))) c(fn="log",inv="exp")
                     else if (identical(x, c(lower=0, upper=1))) c(fn="qlogis",inv="plogis")
                     else NULL
                 }
                 ) )

### Distinct labelled (1 and) 2 and 3 state directed graphs.
### graphs with common "iso" are isomorphic (i.e. identical when states are unlabelled)
### "perm" is permutation of states needed to transform graph into the first in the list of isomorphisms
### This database is used to determine the appropriate method for calculating the analytic P matrix.

### The numbered label gives the indices into the matrix of rates (vectorised by reading across rows)
### e.g. the model with qmatrix of the form
### *,1,1
### 0,*,1
### 0,0,*    is "1-2-4"
### well-disease, well-death, disease-death transitions allowed.

.msm.graphs <-
  list(
       "1" = list(),
       "2" = list(
         "1" = list(iso=1, perm=c(1,2)),
         "2" = list(iso=1, perm=c(2,1)),
         "1-2" = list(iso=2, perm=c(1,2))
         ),
       "3" = list(
         "1-2" = list(iso=1, perm=c(1,2,3)),
         "3-4" = list(iso=1, perm=c(3,1,2)),
         "5-6" = list(iso=1, perm=c(2,3,1)),

         "1-4" = list(iso=2, perm=c(1,2,3)),
         "1-5" = list(iso=2, perm=c(2,3,1)),
         "2-3" = list(iso=2, perm=c(2,1,3)),
         "2-6" = list(iso=2, perm=c(1,3,2)),
         "3-6" = list(iso=2, perm=c(3,2,1)),
         "4-5" = list(iso=2, perm=c(3,1,2)),

         "1-6" = list(iso=3, perm=c(1,2,3)),
         "2-4" = list(iso=3, perm=c(3,1,2)),
         "3-5" = list(iso=3, perm=c(2,3,1)),

         "1-2-4" = list(iso=4,perm=c(1,2,3)),
         "1-2-6" = list(iso=4,perm=c(1,3,2)),
         "1-5-6" = list(iso=4,perm=c(2,3,1)),
         "2-3-4" = list(iso=4,perm=c(2,1,3)),
         "3-4-5" = list(iso=4,perm=c(3,1,2)),
         "3-5-6" = list(iso=4,perm=c(3,2,1)),

         "1-3-5" = list(iso=5,perm=c(1,2,3)),
         "1-3-6" = list(iso=5,perm=c(2,1,3)),
         "1-4-6" = list(iso=5,perm=c(3,1,2)),
         "2-3-5" = list(iso=5,perm=c(1,3,2)),
         "2-4-5" = list(iso=5,perm=c(2,3,1)),
         "2-4-6" = list(iso=5,perm=c(3,2,1)),

         "1-2-4-6" = list(iso=6,perm=c(1,2,3)),
         "1-3-5-6" = list(iso=6,perm=c(2,3,1)),
         "2-3-4-5" = list(iso=6,perm=c(3,1,2))
         ),
       "4" = list(
         "1-5-9" = list(iso=1,perm=c(1,2,3,4)),
         "1-3-5-6-9" = list(iso=2,perm=c(1,2,3,4))
         ),

       "5" = list(
         "1-6-11-16" = list(iso=1,perm=c(1,2,3,4,5)),
         "1-4-6-8-11-12-16" = list(iso=2,perm=c(1,2,3,4,5)),
         "1-6-7-11-12" = list(iso=3,perm=c(1,2,3,4,5))
         )

       )

## Tasks to be performed in C

.msm.CTASKS <- c("lik","deriv","info","viterbi","lik.subj","deriv.subj","dpmat")
