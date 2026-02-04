# Changelog

## Version 1.8.2 (2024-11-07)

CRAN release: 2024-11-07

- Fix a further memory error affecting CRAN checks of some dependent
  packages ([\#106](https://github.com/chjackson/msm/issues/106)).

## Version 1.8.1 (2024-10-04)

CRAN release: 2024-10-04

- Fix memory error that crashed reverse dependency checks
  ([\#106](https://github.com/chjackson/msm/issues/106)).

## Version 1.8 (2024-09-06)

CRAN release: 2024-09-08

- Subject-level weights supported for likelihood calculation with new
  [`msm()`](https://chjackson.github.io/msm/reference/msm.md) argument
  `subject.weights`.

- Tidying methods added for `msm` objects and most of msm’s output
  functions (`prevalence.msm`, `qmatrix.msm` and related functions, and
  `totlos.msm` and related functions). These methods convert the output
  to a tidy data frame (tibble), in the manner of the
  [broom](https://broom.tidymodels.org/) package.

  To use these methods, just call `tidy(x)`, where `x` is the result of
  calling, e.g. `msm`, `prevalence.msm`, or `qmatrix.msm`.

  Hence the msm package now imports the `generics` and `tibble`
  packages.

- Subjects with one observation are no longer dropped in HMMs, since
  they provide information about the distribution of the outcome given
  the hidden state.

- `ppass.msm` now supports `pci` models and other time-inhomogeneous
  models. Thanks to Jon Fintzi for working on this.

- New function `hmodel2list` to extract HMM constructor function calls
  from fitted HMMs. Thanks to Will Hulme for working on this.

- Objects returned by `totlos.msm`, `efpt.msm` and `envisits.msm` now
  have class `"msm.estbystate"`.

- Objects returned by `prevalence.msm` now have class
  `"msm.prevalence"`.

- Fix of bugs for models containing a covariate named `"baseline"` or
  `"Baseline"`.

## Version 1.7.1 (2023-11-23)

CRAN release: 2023-11-23

- Fix of a bug in the Viterbi algorithm for the calculation of the
  fitted state at the initial time for each subject.

- Bug fix for `pmatrix.piecewise.msm` given just intensity matrices
  instead of a fitted model.

- Auto-generated initial values set to a small positive number rather
  than zero when there are no observed data for a particular permitted
  transition. This fixes consistency checks (e.g. for `qconstraint`) in
  this situation.

- Fix when covariates come into the data as one-column matrices instead
  of vectors.

- Modernised to use roxygen, and pkgdown website created at
  <https://chjackson.github.io/msm>

- Some internal functions (e.g. `deriv.msm`) that clash with base S3
  generics renamed

## Version 1.7 (2022-11-27)

CRAN release: 2022-11-28

- `rpexp` is now more efficient. Thanks to Mark Clements. Note that the
  values simulated by `rpexp`, `sim.msm` and `simmulti.msm` for a
  specific seed will now be different in models where the intensities
  are piecewise constant.

- Bug fix: in HMMs with partially known initial states specified through
  `obstrue` and `censor`, `initprobs` is now accounted for.

- Bug fix for bootstrapping with character subject IDs.

- `crudeinits.msm` handled when state is a factor.

- Fix for `hranges` on multiple parameters.

- `obstrue` handled in forward-backward algorithm (`viterbi.msm`).

## Version 1.6.9

CRAN release: 2021-09-27

- … and earlier versions: see
  [inst/NEWS](https://github.com/chjackson/msm/blob/master/inst/NEWS) in
  the source for changes.
