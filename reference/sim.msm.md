# Simulate one individual trajectory from a continuous-time Markov model

Simulate the exact times of transition of a continuous-time Markov
process up to a given time.

## Usage

``` r
sim.msm(
  qmatrix,
  maxtime,
  covs = NULL,
  beta = NULL,
  obstimes = 0,
  start = 1,
  mintime = 0
)
```

## Arguments

- qmatrix:

  The transition intensity matrix of the Markov process. The diagonal of
  `qmatrix` is ignored, and computed as appropriate so that the rows sum
  to zero. For example, a possible `qmatrix` for a three state
  illness-death model with recovery is:

  `rbind( c( 0, 0.1, 0.02 ), c( 0.1, 0, 0.01 ), c( 0, 0, 0 ) )`

- maxtime:

  Maximum time for the simulated process.

- covs:

  Matrix of time-dependent covariates, with one row for each observation
  time and one column for each covariate.

- beta:

  Matrix of linear covariate effects on log transition intensities. The
  rows correspond to different covariates, and the columns to the
  transition intensities. The intensities are ordered by reading across
  rows of the intensity matrix, starting with the first, counting the
  positive off-diagonal elements of the matrix.

- obstimes:

  Vector of times at which the covariates are observed.

- start:

  Starting state of the process. Defaults to 1.

- mintime:

  Starting time of the process. Defaults to 0.

## Value

A list with components,

- states:

  Simulated states through which the process moves. This ends with
  either an absorption before `obstime`, or a transient state at
  `obstime`.

- times:

  Exact times at which the process changes to the corresponding states

- qmatrix:

  The given transition intensity matrix

## Details

The effect of time-dependent covariates on the transition intensity
matrix for an individual is determined by assuming that the covariate is
a step function which remains constant in between the individual's
observation times.

## See also

[`simmulti.msm`](https://chjackson.github.io/msm/reference/simmulti.msm.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r

qmatrix <- rbind(
                 c(-0.2,   0.1,  0.1 ),
                 c(0.5,   -0.6,  0.1 ),
                 c(0,  0,  0)
                 )
sim.msm(qmatrix, 30)
#> $states
#> [1] 1 2 1 2 1 2 1 1
#> 
#> $times
#> [1]  0.0000000  0.2709619  2.5593627 13.9154581 16.4764457 17.4488506 19.5114268
#> [8] 30.0000000
#> 
#> $qmatrix
#>      [,1] [,2] [,3]
#> [1,] -0.2  0.1  0.1
#> [2,]  0.5 -0.6  0.1
#> [3,]  0.0  0.0  0.0
#> 
```
