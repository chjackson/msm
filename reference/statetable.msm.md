# Table of transitions

Calculates a frequency table counting the number of times each pair of
states were observed in successive observation times. This can be a
useful way of summarising multi-state data.

## Usage

``` r
statetable.msm(state, subject, data = NULL)
```

## Arguments

- state:

  Observed states, assumed to be ordered by time within each subject.

- subject:

  Subject identification numbers corresponding to `state`. If not given,
  all observations are assumed to be on the same subject.

- data:

  An optional data frame in which the variables represented by `subject`
  and `state` can be found.

## Value

A frequency table with starting states as rows and finishing states as
columns.

## Details

If the data are intermittently observed (panel data) this table should
not be used to decide what transitions should be allowed in the \\Q\\
matrix, which works in continuous time. This function counts the
transitions between states over a time interval, not in real time. There
can be observed transitions between state \\r\\ and \\s\\ over an
interval even if \\q\_{rs}=0\\, because the process may have passed
through one or more intermediate states in the middle of the interval.

## See also

[`crudeinits.msm`](https://chjackson.github.io/msm/reference/crudeinits.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r
## Heart transplant data
data(cav)
#> Warning: data set ‘cav’ not found

## 148 deaths from state 1, 48 from state 2 and 55 from state 3.
statetable.msm(state, PTNUM, data=cav)
#>     to
#> from    1    2    3    4
#>    1 1367  204   44  148
#>    2   46  134   54   48
#>    3    4   13  107   55

```
