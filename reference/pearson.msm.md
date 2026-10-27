# Pearson-type goodness-of-fit test

Pearson-type goodness-of-fit test for multi-state models fitted to
panel-observed data.

## Usage

``` r
pearson.msm(
  x,
  transitions = NULL,
  timegroups = 3,
  intervalgroups = 3,
  covgroups = 3,
  groups = NULL,
  boot = FALSE,
  B = 500,
  next.obstime = NULL,
  N = 100,
  indep.cens = TRUE,
  maxtimes = NULL,
  pval = TRUE
)
```

## Arguments

- x:

  A fitted multi-state model, as returned by
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- transitions:

  This should be an integer vector indicating which interval transitions
  should be grouped together in the contingency table. Its length should
  be the number of allowed interval transitions, excluding transitions
  from absorbing states to absorbing states.

  The allowed interval transitions are the set of pairs of states
  \\(a,b)\\ for which it is possible to observe \\a\\ at one time and
  \\b\\ at any later time. For example, in a "well-disease-death" model
  with allowed *instantaneous* 1-2, 2-3 transitions, there are 5 allowed
  *interval* transitions. In numerical order, these are 1-1, 1-2, 1-3,
  2-2 and 2-3, excluding absorbing-absorbing transitions.

  Then, to group transitions 1-1,1-2 together, and transitions 2-2,2-3
  together, specify

  `transitions = c(1,1,2,3,3)`.

  Only transitions from the same state may be grouped. By default, each
  interval transition forms a separate group.

- timegroups:

  Number of groups based on quantiles of the time since the start of the
  process.

- intervalgroups:

  Number of groups based on quantiles of the time interval between
  observations, within time groups

- covgroups:

  Number of groups based on quantiles of \\\sum_r q\_{irr}\\, where
  \\q\_{irr}\\ are the diagonal entries of the transition intensity
  matrix for the *i*th transition. These are a function of the covariate
  effects and the covariate values at the *i*th transition: \\q\_{irr}\\
  is minus the sum of the off-diagonal entries \\q\_{rs}^{(0)} exp
  (\beta\_{rs}^T z_i)\\ on the *r*th row.

  Thus `covgroups` summarises the impact of covariates at each
  observation, by calculating the overall rate of progression through
  states at that observation.

  For time-inhomogeneous models specified using the `pci` argument to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), if the only
  covariate is the time period, `covgroups` is set to 1, since
  `timegroups` ensures that transitions are grouped by time.

- groups:

  A vector of arbitrary groups in which to categorise each transition.
  This can be an integer vector or a factor. This can be used to
  diagnose specific areas of poor fit. For example, the contingency
  table might be grouped by arbitrary combinations of covariates to
  detect types of individual for whom the model fits poorly.

  The length of `groups` should be `x$data$n`, the number of
  observations used in the model fit, which is the number of
  observations in the original dataset with any missing values excluded.
  The value of `groups` at observation \\i\\ is used to categorise the
  transition which *ends* at observation i. Values of `groups` at the
  first observation for each subject are ignored.

- boot:

  Estimate an "exact" p-value using a parametric bootstrap.

  All objects used in the original call to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md) which
  produced `x`, such as the `qmatrix`, should be in the working
  environment, or else an “object not found” error will be given. This
  enables the original model to be refitted to the replicate datasets.

  Note that `groups` cannot be used with bootstrapping, as the simulated
  observations will not be in the same categories as the original
  observations.

- B:

  Number of bootstrap replicates.

- next.obstime:

  This is a vector of length `nrow(model.frame(x))` (the number of
  observations used in the model fit) giving the time to the next
  *scheduled* observation following each time point. This is only used
  when times to death are known exactly.

  For individuals who died (entered an absorbing state) before the next
  scheduled observation, and the time of death is known exactly,
  `next.obstime` would be *greater* than the observed death time.

  If the individual did not die, and a scheduled observation did follow
  that time point, `next.obstime` should just be the same as the time to
  that observation.

  `next.obstime` is used to determine a grouping of the time interval
  between observations, which should be based on scheduled observations.
  If exact times to death were used in the grouping, then shorter
  intervals would contain excess deaths, and the goodness-of-fit
  statistic would be biased.

  If `next.obstime` is unknown, it is multiply-imputed using a
  product-limit estimate based on the intervals to observations other
  than deaths. The resulting tables of transitions are averaged over
  these imputations. This may be slow.

- N:

  Number of imputations for the estimation of the distribution of the
  next scheduled observation time, when there are exact death times.

- indep.cens:

  If `TRUE`, then times to censoring are included in the estimation of
  the distribution to the next scheduled observation time. If `FALSE`,
  times to censoring are assumed to be systematically different from
  other observation times.

- maxtimes:

  A vector of length `x$data$n`, or a common scalar, giving an upper
  bound for the next scheduled observation time. Used in the multiple
  imputation when times to death are known exactly. If a value greater
  than `maxtimes` is simulated, then the next scheduled observation is
  taken as censored. This should be supplied, if known. If not supplied,
  this is taken to be the maximum interval occurring in the data, plus
  one time unit. For observations which are not exact death times, this
  should be the time since the previous observation.

- pval:

  Calculate a p-value using the improved approximation of Titman (2009).
  This is optional since it is not needed during bootstrapping, and it
  is computationally non-trivial. Only available currently for
  non-hidden Markov models for panel data without exact death times.
  Also not available for models with censoring, including
  time-homogeneous models fitted with the `pci` option to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

## Value

A list whose first two elements are contingency tables of observed
transitions \\O\\ and expected transitions \\E\\, respectively, for each
combination of groups. The third element is a table of the deviances
\\(O - E)^2 / E\\ multiplied by the sign of \\O - E\\. If the expected
number of transitions is zero then the deviance is zero. Entries in the
third matrix will be bigger in magnitude for groups for which the model
fits poorly.  

- list("\\test\\"):

  the fourth element of the list, is a data frame with one row
  containing the Pearson-type goodness-of-fit test statistic `stat`. The
  test statistic is the sum of the deviances. For panel-observed data
  without exact death times, misclassification or censored observations,
  `p` is the p-value for the test statistic calculated using the
  improved approximation of Titman (2009).

  For these models, for comparison with older versions of the package,
  `test` also presents `p.lower` and `p.upper`, which are theoretical
  lower and upper limits for the p-value of the test statistic, based on
  chi-squared distributions with `df.lower` and `df.upper` degrees of
  freedom, respectively. `df.upper` is the number of independent cells
  in the contingency table, and `df.lower` is `df.upper` minus the
  number of estimated parameters in the model.

- list("\\intervalq\\"):

  (not printed by default) contains the definition of the grouping of
  the intervals between observations. These groups are defined by
  quantiles within the groups corresponding to the time since the start
  of the process.

- list("\\sim\\"):

  If there are exact death times, this contains simulations of the
  contingency tables and test statistics for each imputation of the next
  scheduled sampling time. These are averaged over to produce the
  presented tables and test statistic. This element is not printed by
  default.

  With exact death times, the null variance of the test statistic
  (formed by taking mean of simulated test statistics) is less than
  twice the mean (Titman, 2008), and the null distribution is not
  chi-squared. In this case, `p.upper` is an upper limit for the true
  asymptotic p-value, but `p.lower` is not a lower limit, and is not
  presented.

- list("\\boot\\"):

  If the bootstrap has been used, the element will contain the bootstrap
  replicates of the test statistics (not printed by default).

- list("\\lambda\\"):

  If the Titman (2009) p-value has been calculated, this contains the
  weights defining the null distribution of the test statistic as a
  weighted sum of chi-squared(1) random variables (not printed by
  default).

## Details

This method (Aguirre-Hernandez and Farewell, 2002) is intended for data
which represent observations of the process at arbitrary times
("snapshots", or "panel-observed" data). For data which represent the
exact transition times of the process,
[`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md)
can be used to assess fit, though without a formal test.

When times of death are known exactly, states are misclassified, or an
individual's final observation is a censored state, the modification by
Titman and Sharples (2008) is used. The only form of censoring supported
is a state at the end of an individual's series which represents an
unknown transient state (i.e. the individual is only known to be alive
at this time). Other types of censoring are omitted from the data before
performing the test.

See the references for further details of the methods. The method used
for censored states is a modification of the method in the appendix to
Titman and Sharples (2008), described at
<https://chjackson.github.io/msm/misc/robustcensoring.pdf> (Titman,
2007).

Groupings of the time since initiation, the time interval and the impact
of covariates are based on equally-spaced quantiles. The number of
groups should be chosen that there are not many cells with small
expected numbers of transitions, since the deviance statistic will be
unstable for sparse contingency tables. Ideally, the expected numbers of
transitions in each cell of the table should be no less than about 5.
Conversely, the power of the test is reduced if there are too few
groups. Therefore, some sensitivity analysis of the test results to the
grouping is advisable.

Saved model objects fitted with previous versions of R (versions less
than 1.2) will need to be refitted under the current R for use with
`pearson.msm`.

## References

Aguirre-Hernandez, R. and Farewell, V. (2002) A Pearson-type
goodness-of-fit test for stationary and time-continuous Markov
regression models. *Statistics in Medicine* 21:1899-1911.

Titman, A. and Sharples, L. (2008) A general goodness-of-fit test for
Markov and hidden Markov models. *Statistics in Medicine*
27(12):2177-2195

Titman, A. (2009) Computation of the asymptotic null distribution of
goodness-of-fit tests for multi-state models. *Lifetime Data Analysis*
15(4):519-533.

Titman, A. (2008) Model diagnostics in multi-state models of biological
systems. PhD thesis, University of Cambridge.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),
[`prevalence.msm`](https://chjackson.github.io/msm/reference/prevalence.msm.md),
[`scoreresid.msm`](https://chjackson.github.io/msm/reference/scoreresid.msm.md),

## Author

Andrew Titman <a.titman@lancaster.ac.uk>, Chris Jackson
<chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r

psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
psor.msm <- msm(state ~ months, subject=ptnum, data=psor,
                qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
pearson.msm(psor.msm, timegroups=2, intervalgroups=2, covgroups=2)
#> $Observed
#>           Time Interval           Cov 1-1 1-2 1-3 1-4 2-2 2-3 2-4 3-3 3-4
#> 1   [0.4,9.13)        1 [-0.81,-0.49)  20  11   1   0  11   8   2   4   4
#> 2 [9.13,57.06]        1 [-0.81,-0.49)   4   6   0   1   8   7   1   7  10
#> 3   [0.4,9.13)        2 [-0.81,-0.49)  26   2   1   2   9   1   0   4   3
#> 4 [9.13,57.06]        2 [-0.81,-0.49)  22   0   0   1  13   1   4   6   1
#> 5   [0.4,9.13)        1         -0.49  10  12   5   2  10   5   1   4   2
#> 6 [9.13,57.06]        1         -0.49  10   6   4   0  15   8   6   9  10
#> 7   [0.4,9.13)        2         -0.49  38   8   1   2  11   2   0   2   0
#> 8 [9.13,57.06]        2         -0.49  29   3   4   0  12   2   1  10   4
#> 
#> $Expected
#>           Time Interval           Cov       1-1       1-2       1-3         1-4
#> 1   [0.4,9.13)        1 [-0.81,-0.49) 30.153882  1.711383 0.1197351 0.015000170
#> 2 [9.13,57.06]        1 [-0.81,-0.49)  9.659011  1.093258 0.1922818 0.055449096
#> 3   [0.4,9.13)        2 [-0.81,-0.49) 24.446257  4.787454 1.1078563 0.658432763
#> 4 [9.13,57.06]        2 [-0.81,-0.49) 13.434261  5.565404 1.6534384 2.346895857
#> 5   [0.4,9.13)        1         -0.49 26.931707  1.917908 0.1409192 0.009465852
#> 6 [9.13,57.06]        1         -0.49 17.132490  2.390072 0.4047735 0.072664384
#> 7   [0.4,9.13)        2         -0.49 34.360164 10.126221 3.2263497 1.287265020
#> 8 [9.13,57.06]        2         -0.49 14.959359  9.248204 5.9052553 5.887182158
#>         2-2      2-3       2-4       3-3       3-4
#> 1 18.230441 2.383741 0.3858178  6.177131 1.8228694
#> 2 12.698334 2.268765 1.0329011 12.581727 4.4182729
#> 3  6.235941 2.057100 1.7069592  2.522098 4.4779016
#> 4  6.363320 2.856426 8.7802543  1.249858 5.7501417
#> 5 14.356291 1.524817 0.1188923  5.326063 0.6739374
#> 6 22.208605 5.504201 1.2871938 15.066175 3.9338249
#> 7  7.462942 3.827703 1.7093551  1.338348 0.6616517
#> 8  3.967780 4.405898 6.6263211  4.473459 9.5265407
#> 
#> $`Deviance*sign(O-E)`
#>           Time Interval           Cov          1-1        1-2          1-3
#> 1   [0.4,9.13)        1 [-0.81,-0.49)  -3.41917244 50.4144530   6.47150280
#> 2 [9.13,57.06]        1 [-0.81,-0.49)  -3.31549503 22.0223437  -0.19228180
#> 3   [0.4,9.13)        2 [-0.81,-0.49)   0.09875205 -1.6229713  -0.01050044
#> 4 [9.13,57.06]        2 [-0.81,-0.49)   5.46154919 -5.5654044  -1.65343837
#> 5   [0.4,9.13)        1         -0.49 -10.64480283 52.9997341 167.54749538
#> 6 [9.13,57.06]        1         -0.49  -2.96935320  5.4523792  31.93305610
#> 7   [0.4,9.13)        2         -0.49   0.38557452 -0.4464465  -1.53629750
#> 8 [9.13,57.06]        2         -0.49  13.17834649 -4.2213659  -0.61470632
#>            1-4       2-2        2-3          2-4        3-3        3-4
#> 1  -0.01500017 -2.867692 13.2322971  6.753405695 -0.7673300  2.6002400
#> 2  16.09000826 -1.738365  9.8664183 -0.001047999 -2.4762640  7.0515511
#> 3   2.73346460  1.225160 -0.5432215 -1.706959185  0.8660222 -0.4877716
#> 4  -0.77299061  6.921783 -1.2065135 -2.602524980 18.0531236 -3.9240505
#> 5 418.58102797 -1.321878  7.9202290  6.529864884 -0.3301580  2.6092071
#> 6  -0.07266438 -2.339813  1.1316836 17.255010981 -2.4424567  9.3543769
#> 7   0.39462826  1.676387 -0.8727163 -1.709355069  0.3271069 -0.6616517
#> 8  -5.88718216 16.260111 -1.3137723 -4.777234415  6.8275242 -3.2060590
#> 
#> $test
#>      stat p df.lower p.lower df.upper p.upper
#>  1010.483 0       42       0       48       0
#> 
# More 1-2, 1-3 and 1-4 observations than expected in shorter time
# intervals - the model fits poorly.
# A random effects model might accommodate such fast progressors.
```
