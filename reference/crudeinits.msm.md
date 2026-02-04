# Calculate crude initial values for transition intensities

Calculates crude initial values for transition intensities by assuming
that the data represent the exact transition times of the Markov
process.

## Usage

``` r
crudeinits.msm(
  formula,
  subject,
  qmatrix,
  data = NULL,
  censor = NULL,
  censor.states = NULL
)
```

## Arguments

- formula:

  A formula giving the vectors containing the observed states and the
  corresponding observation times. For example,

  `state ~ time`

  Observed states should be in the set `1, ...{}, n`, where `n` is the
  number of states. Note hidden Markov models are not supported by this
  function.

- subject:

  Vector of subject identification numbers for the data specified by
  `formula`. If missing, then all observations are assumed to be on the
  same subject. These must be sorted so that all observations on the
  same subject are adjacent.

- qmatrix:

  Matrix of indicators for the allowed transitions. An initial value
  will be estimated for each value of qmatrix that is greater than zero.
  Transitions are taken as disallowed for each entry of `qmatrix` that
  is 0.

- data:

  An optional data frame in which the variables represented by `subject`
  and `state` can be found.

- censor:

  A state, or vector of states, which indicates censoring. See
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- censor.states:

  Specifies the underlying states which censored observations can
  represent. See
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

## Value

The estimated transition intensity matrix. This can be used as the
`qmatrix` argument to
[`msm`](https://chjackson.github.io/msm/reference/msm.md).

## Details

Suppose we want a crude estimate of the transition intensity \\q\_{rs}\\
from state \\r\\ to state \\s\\. If we observe \\n\_{rs}\\ transitions
from state \\r\\ to state \\s\\, and a total of \\n_r\\ transitions from
state \\r\\, then \\q\_{rs} / \\\\ q\_{rr}\\ can be estimated by
\\n\_{rs} / n_r\\. Then, given a total of \\T_r\\ years spent in state
\\r\\, the mean sojourn time \\1 / q\_{rr}\\ can be estimated as \\T_r /
n_r\\. Thus, \\n\_{rs} / T_r\\ is a crude estimate of \\q\_{rs}\\.

If the data do represent the exact transition times of the Markov
process, then these are the exact maximum likelihood estimates.

Observed transitions which are incompatible with the given `qmatrix` are
ignored. Censored states are ignored.

## See also

[`statetable.msm`](https://chjackson.github.io/msm/reference/statetable.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r
data(cav)
#> Warning: data set ‘cav’ not found
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
statetable.msm(state, PTNUM, data=cav)
#>     to
#> from    1    2    3    4
#>    1 1367  204   44  148
#>    2   46  134   54   48
#>    3    4   13  107   55
crudeinits.msm(state ~ years, PTNUM, data=cav, qmatrix=twoway4.q)
#>            [,1]        [,2]       [,3]       [,4]
#> [1,] -0.1173149  0.06798932  0.0000000 0.04932559
#> [2,]  0.1168179 -0.37584883  0.1371340 0.12189692
#> [3,]  0.0000000  0.04908401 -0.2567471 0.20766310
#> [4,]  0.0000000  0.00000000  0.0000000 0.00000000
```
