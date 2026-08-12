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
#>           Time Interval           Cov      1-1       1-2       1-3         1-4
#> 1   [0.4,9.13)        1 [-0.81,-0.49) 30.15388  1.711383 0.1197352 0.015000173
#> 2 [9.13,57.06]        1 [-0.81,-0.49)  9.65901  1.093259 0.1922819 0.055449102
#> 3   [0.4,9.13)        2 [-0.81,-0.49) 24.44626  4.787454 1.1078571 0.658432845
#> 4 [9.13,57.06]        2 [-0.81,-0.49) 13.43426  5.565404 1.6534395 2.346896197
#> 5   [0.4,9.13)        1         -0.49 26.93171  1.917908 0.1409193 0.009465849
#> 6 [9.13,57.06]        1         -0.49 17.13249  2.390072 0.4047737 0.072664358
#> 7   [0.4,9.13)        2         -0.49 34.36016 10.126221 3.2263515 1.287264637
#> 8 [9.13,57.06]        2         -0.49 14.95936  9.248203 5.9052593 5.887180924
#>         2-2      2-3       2-4       3-3       3-4
#> 1 18.230440 2.383742 0.3858178  6.177131 1.8228685
#> 2 12.698333 2.268766 1.0329010 12.581729 4.4182708
#> 3  6.235940 2.057101 1.7069591  2.522100 4.4779001
#> 4  6.363318 2.856427 8.7802550  1.249859 5.7501405
#> 5 14.356290 1.524817 0.1188923  5.326063 0.6739368
#> 6 22.208604 5.504203 1.2871932 15.066178 3.9338221
#> 7  7.462941 3.827705 1.7093543  1.338349 0.6616512
#> 8  3.967779 4.405902 6.6263193  4.473463 9.5265366
#> 
#> $`Deviance*sign(O-E)`
#>           Time Interval           Cov          1-1        1-2          1-3
#> 1   [0.4,9.13)        1 [-0.81,-0.49)  -3.41917220 50.4144389   6.47149697
#> 2 [9.13,57.06]        1 [-0.81,-0.49)  -3.31549478 22.0223394  -0.19228194
#> 3   [0.4,9.13)        2 [-0.81,-0.49)   0.09875221 -1.6229716  -0.01050059
#> 4 [9.13,57.06]        2 [-0.81,-0.49)   5.46155139 -5.5654043  -1.65343949
#> 5   [0.4,9.13)        1         -0.49 -10.64480259 52.9997257 167.54741678
#> 6 [9.13,57.06]        1         -0.49  -2.96935296  5.4523782  31.93303690
#> 7   [0.4,9.13)        2         -0.49   0.38557488 -0.4464466  -1.53629913
#> 8 [9.13,57.06]        2         -0.49  13.17835135 -4.2213650  -0.61470850
#>            1-4       2-2        2-3          2-4        3-3        3-4
#> 1  -0.01500017 -2.867691 13.2322860  6.753406750 -0.7673305  2.6002432
#> 2  16.09000614 -1.738365  9.8664089 -0.001047993 -2.4762654  7.0515597
#> 3   2.73346393  1.225161 -0.5432223 -1.706959075  0.8660200 -0.4877708
#> 4  -0.77299088  6.921790 -1.2065149 -2.602525535 18.0530980 -3.9240493
#> 5 418.58119783 -1.321878  7.9202243  6.529869316 -0.3301582  2.6092113
#> 6  -0.07266436 -2.339813  1.1316812 17.255023705 -2.4424585  9.3543918
#> 7   0.39462880  1.676389 -0.8727176 -1.709354344  0.3271063 -0.6616512
#> 8  -5.88718092 16.260122 -1.3137747 -4.777232650  6.8275078 -3.2060557
#> 
#> $test
#>      stat p df.lower p.lower df.upper p.upper
#>  1010.483 0       42       0       48       0
#> 
# More 1-2, 1-3 and 1-4 observations than expected in shorter time
# intervals - the model fits poorly.
# A random effects model might accommodate such fast progressors.
```
