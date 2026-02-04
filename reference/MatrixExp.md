# Matrix exponential

Calculates the exponential of a square matrix.

## Usage

``` r
MatrixExp(mat, t = 1, method = NULL, ...)
```

## Arguments

- mat:

  A square matrix

- t:

  An optional scaling factor for `mat`. If a vector is supplied, then an
  array of matrices is returned with different scaling factors.

- method:

  Under the default of `NULL`, this simply wraps the
  [`expm`](https://rdrr.io/pkg/expm/man/expm.html) function from the
  expm package. This is recommended. Options to
  [`expm`](https://rdrr.io/pkg/expm/man/expm.html) can be supplied to
  `MatrixExp`, including `method`.

  Otherwise, for backwards compatibility, the following options, which
  use code in the msm package, are available: `"pade"` for a Pade
  approximation method, `"series"` for the power series approximation,
  or `"analytic"` for the analytic formulae for simpler Markov model
  intensity matrices (see below). These options are only used if `mat`
  has repeated eigenvalues, thus the usual eigen-decomposition method
  cannot be used.

- ...:

  Arguments to pass to [`expm`](https://rdrr.io/pkg/expm/man/expm.html).

## Value

The exponentiated matrix \\\exp(mat)\\. Or, if `t` is a vector of length
2 or more, an array of exponentiated matrices.

## Details

See the [`expm`](https://rdrr.io/pkg/expm/man/expm.html) documentation
for details of the algorithms it uses.

Generally the exponential \\E\\ of a square matrix \\M\\ can often be
calculated as

\$\$E = U \exp(D) U^{-1}\$\$

where \\D\\ is a diagonal matrix with the eigenvalues of \\M\\ on the
diagonal, \\\exp(D)\\ is a diagonal matrix with the exponentiated
eigenvalues of \\M\\ on the diagonal, and \\U\\ is a matrix whose
columns are the eigenvectors of \\M\\.

This method of calculation is used if `"pade"` or `"series"` is supplied
but \\M\\ has distinct eigenvalues. I If \\M\\ has repeated eigenvalues,
then its eigenvector matrix may be non-invertible. In this case, the
matrix exponential is calculated using the Pade approximation defined by
Moler and van Loan (2003), or the less robust power series
approximation,

\$\$\exp(M) = I + M + M^2/2 + M^3 / 3! + M^4 / 4! + ...\$\$

For a continuous-time homogeneous Markov process with transition
intensity matrix \\Q\\, the probability of occupying state \\s\\ at time
\\u + t\\ conditional on occupying state \\r\\ at time \\u\\ is given by
the \\(r,s)\\ entry of the matrix \\\exp(tQ)\\.

If `mat` is a valid transition intensity matrix for a continuous-time
Markov model (i.e. diagonal entries non-positive, off-diagonal entries
non-negative, rows sum to zero), then for certain simpler model
structures, there are analytic formulae for the individual entries of
the exponential of `mat`. These structures are listed in the PDF manual
and the formulae are coded in the msm source file `src/analyticp.c`.
These formulae are only used if `method="analytic"`. This is more
efficient, but it is not the default in `MatrixExp` because the code is
not robust to extreme values. However it is the default when calculating
likelihoods for models fitted by
[`msm`](https://chjackson.github.io/msm/reference/msm.md).

The implementation of the Pade approximation used by `method="pade"` was
taken from JAGS by Martyn Plummer (<https://mcmc-jags.sourceforge.io>).

## References

Cox, D. R. and Miller, H. D. *The theory of stochastic processes*,
Chapman and Hall, London (1965)

Moler, C and van Loan, C (2003). Nineteen dubious ways to compute the
exponential of a matrix, twenty-five years later. *SIAM Review* **45**,
3–49.
