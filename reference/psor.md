# Psoriatic arthritis data

A series of observations of grades of psoriatic arthritis, as indicated
by numbers of damaged joints.

## Format

A data frame containing 806 observations, representing visits to a
psoriatic arthritis (PsA) clinic from 305 patients. The rows are grouped
by patient number and ordered by examination time. Each row represents
an examination and contains additional covariates.

|            |             |                                                         |
|------------|-------------|---------------------------------------------------------|
| `ptnum`    | (numeric)   | Patient identification number                           |
| `months`   | (numeric)   | Examination time in months                              |
| `state`    | (numeric)   | Clinical state of PsA. Patients in states 1, 2, 3 and 4 |
|            |             | have 0, 1 to 4, 5 to 9 and 10 or more damaged joints,   |
|            |             | respectively.                                           |
| `hieffusn` | (numeric)   | Presence of five or more effusions                      |
| `ollwsdrt` | (character) | Erythrocyte sedimentation rate of less than 15 mm/h     |

## References

Gladman, D. D. and Farewell, V.T. (1999) Progression in psoriatic
arthritis: role of time-varying clinical indicators. J. Rheumatol.
26(11):2409-13

## Examples

``` r
## Four-state progression-only model with high effusion and low
## sedimentation rate as covariates on the progression rates.  High
## effusion is assumed to have the same effect on the 1-2, 2-3, and 3-4
## progression rates, while low sedimentation rate has the same effect
## on the 1-2 and 2-3 intensities, but a different effect on the 3-4. 

data(psor)
#> Warning: data set ‘psor’ not found
psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
psor.msm <- msm(state ~ months, subject=ptnum, data=psor, 
                qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)),
                fixedpars=FALSE, control = list(REPORT=1,trace=2), method="BFGS")
#> initial  value 1184.216999 
#> iter   2 value 1127.501356
#> iter   3 value 1122.654955
#> iter   4 value 1121.606113
#> iter   5 value 1120.763406
#> iter   6 value 1119.769934
#> iter   7 value 1116.747874
#> iter   8 value 1116.596341
#> iter   9 value 1114.972649
#> iter  10 value 1114.899884
#> iter  11 value 1114.899464
#> iter  11 value 1114.899461
#> iter  11 value 1114.899461
#> final  value 1114.899461 
#> converged
#> Used 38 function and 11 gradient evaluations
qmatrix.msm(psor.msm)
#>         State 1                    State 2                   
#> State 1 -0.09594 (-0.1216,-0.0757)  0.09594 ( 0.0757, 0.1216)
#> State 2 0                          -0.16431 (-0.2076,-0.1300)
#> State 3 0                          0                         
#> State 4 0                          0                         
#>         State 3                    State 4                   
#> State 1 0                          0                         
#> State 2  0.16431 ( 0.1300, 0.2076) 0                         
#> State 3 -0.25438 (-0.3396,-0.1905)  0.25438 ( 0.1905, 0.3396)
#> State 4 0                          0                         
sojourn.msm(psor.msm)
#>         estimates        SE        L         U
#> State 1 10.423724 1.2597644 8.225277 13.209771
#> State 2  6.086186 0.7266461 4.816349  7.690816
#> State 3  3.931084 0.5796054 2.944488  5.248254
hazard.msm(psor.msm)
#> $ollwsdrt
#>                          HR         L        U
#> State 1 - State 2 0.5651903 0.3853452 0.828971
#> State 2 - State 3 0.5651903 0.3853452 0.828971
#> State 3 - State 4 1.6407662 0.8154000 3.301587
#> 
#> $hieffusn
#>                         HR        L        U
#> State 1 - State 2 1.645956 1.148294 2.359299
#> State 2 - State 3 1.645956 1.148294 2.359299
#> State 3 - State 4 1.645956 1.148294 2.359299
#> 
```
