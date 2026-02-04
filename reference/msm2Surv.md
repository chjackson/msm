# Convert data for \`msm' to data for \`survival', \`mstate' or \`flexsurv' analysis

Converts longitudinal data for a
[`msm`](https://chjackson.github.io/msm/reference/msm.md) model fit,
where observations represent the exact transition times of the process,
to counting process data. This enables, for example, flexible parametric
multi-state models to be fitted with
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
from the flexsurv package, or semiparametric models to be implemented
with [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) and the
mstate package.

## Usage

``` r
msm2Surv(data, subject, time, state, covs = NULL, Q)
```

## Arguments

- data:

  Data frame in the format expected by a
  [`msm`](https://chjackson.github.io/msm/reference/msm.md) model fit
  with `exacttimes=TRUE` or all `obstype=2`. Each row represents an
  observation of a state, and the time variable contains the exact and
  complete transition times of the underlying process. This is explained
  in more detail in the help page for
  [`msm`](https://chjackson.github.io/msm/reference/msm.md), section
  `obstype=2`.

- subject:

  Name of the subject ID in the data (character format, i.e. quoted).

- time:

  Name of the time variable in the data (character).

- state:

  Name of the state variable in the data (character).

- covs:

  Vector of covariate names to carry through (character). If not
  supplied, this is taken to be all remaining variables in the data.

- Q:

  Transition intensity matrix. This should have number of rows and
  number of columns both equal to the number of states. If an
  instantaneous transition is not allowed from state \\r\\ to state
  \\s\\, then `Q` should have \\(r,s)\\ entry 0, otherwise it should be
  non-zero. The diagonal entries are ignored.

## Value

A data frame of class `"msdata"`, with rows representing observed or
censored transitions. There will be one row for each observed transition
in the original data, and additional rows for every potential transition
that could have occurred out of each observed state.

The data frame will have columns called:

- id:

  Subject ID

- from:

  Starting state of the transition

- to:

  Finishing state of the transition

- Tstart:

  The starting time of the transition

- Tstop:

  The finishing time of the transition

- time:

  The time difference = `Tstop` - `Tstart`

- status:

  Event or censoring indicator, with 1 indicating an observed
  transition, and 0 indicating censoring

- trans:

  Transition number

and any remaining columns will represent covariates. Any covariates
whose names clash with the standard variables in the returned data
(`"id"`, `"from"`, `"to"`, `"Tstart"`, `"Tstop"`, `"time"`, `"status"`
or `"trans"`) have `".2"` appended to their names.

The transition matrix in mstate format is stored in the `trans`
attribute of the returned object. See the example code below.

## Details

For example, if the data supplied to
[`msm`](https://chjackson.github.io/msm/reference/msm.md) look like
this:

|        |        |          |       |         |
|--------|--------|----------|-------|---------|
| `subj` | `days` | `status` | `age` | `treat` |
| 1      | 0      | 1        | 66    | 1       |
| 1      | 27     | 2        | 66    | 1       |
| 1      | 75     | 3        | 66    | 1       |
| 1      | 97     | 4        | 66    | 1       |
| 1      | 1106   | 4        | 69    | 1       |
| 2      | 0      | 1        | 49    | 0       |
| 2      | 90     | 2        | 49    | 0       |
| 2      | 1037   | 2        | 51    | 0       |

then the output of `msm2Surv` will be a data frame looking like this:

|      |        |      |          |         |        |          |       |         |         |     |
|------|--------|------|----------|---------|--------|----------|-------|---------|---------|-----|
| `id` | `from` | `to` | `Tstart` | `Tstop` | `time` | `status` | `age` | `treat` | `trans` | 1   |
| 1    | 2      | 0    | 27       | 27      | 1      | 66       | 1     | 1       | 1       | 1   |
| 4    | 0      | 27   | 27       | 0       | 66     | 1        | 2     | 1       | 2       | 3   |
| 27   | 75     | 48   | 1        | 66      | 1      | 3        | 1     | 2       | 4       | 27  |
| 75   | 48     | 0    | 66       | 1       | 4      | 1        | 3     | 4       | 75      | 97  |
| 22   | 1      | 69   | 1        | 5       | 2      | 1        | 2     | 0       | 90      | 90  |
| 1    | 49     | 0    | 1        | 2       | 1      | 4        | 0     | 90      | 90      | 0   |
| 49   | 0      | 2    | 2        | 2       | 3      | 90       | 1037  | 947     | 0       | 49  |
| 0    | 3      | 2    | 2        | 4       | 90     | 1037     | 947   | 0       | 49      | 0   |

At 27 days, subject 1 is observed to move from state 1 to state 2 (first
row, status 1), which means that their potential transition from state 1
to state 4 is censored (second row, status 0).

See the mstate package and the references below for more details of this
data format and using it for semi-parametric multi-state modelling.

## References

Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
Competing risks and multi-state models. *Statistics in Medicine* 26:
2389-2430.

Liesbeth C. de Wreede, Marta Fiocco, Hein Putter (2011). mstate: An R
Package for the Analysis of Competing Risks and Multi-State Models.
*Journal of Statistical Software*, 38(7), 1-30.

Jackson, C. H. (2014). flexsurv: Flexible parametric survival and
multi-state models. R package version 0.5.

## See also

[`msprep`](https://rdrr.io/pkg/mstate/man/msprep.html), in mstate, which
produces data in a similar format, given data in "wide" format with one
row per subject.

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r
msmdat <- data.frame(
 subj = c(1, 1, 1, 1, 1, 2, 2, 2),
 days = c(0, 27, 75, 97, 1106, 0, 90, 1037),
 status = c(1, 2, 3, 4, 4, 1, 2, 2),
 age = c(66, 66, 66, 66, 69, 49, 49, 51),
 treat = c(1, 1, 1, 1, 1, 0, 0, 0)
)
# transitions only allowed to next state up or state 4
Q <- rbind(c(1, 1, 0, 1), 
           c(0, 1, 1, 1),
           c(0, 0, 1, 1),
           c(0, 0, 0, 0))
dat <- msm2Surv(data=msmdat, subject="subj", time="days", state="status", 
         Q=Q)
dat
#> An object of class 'msdata'
#> 
#> Data:
#>   id from to Tstart Tstop time status age treat trans
#> 1  1    1  2      0    27   27      1  66     1     1
#> 2  1    1  4      0    27   27      0  66     1     2
#> 3  1    2  3     27    75   48      1  66     1     3
#> 4  1    2  4     27    75   48      0  66     1     4
#> 5  1    3  4     75    97   22      1  66     1     5
#> 6  2    1  2      0    90   90      1  49     0     1
#> 7  2    1  4      0    90   90      0  49     0     2
#> 8  2    2  3     90  1037  947      0  49     0     3
#> 9  2    2  4     90  1037  947      0  49     0     4
attr(dat, "trans")
#>     to
#> from  1  2  3  4
#>    1 NA  1 NA  2
#>    2 NA NA  3  4
#>    3 NA NA NA  5
#>    4 NA NA NA NA
```
