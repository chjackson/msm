User-visible changes only.  For internal changes, see Github commits.

# Version 1.7.1  (2023-11-23)

* Fix of a bug in the Viterbi algorithm for the calculation of the fitted state at the initial time for each subject.

* Bug fix for `pmatrix.piecewise.msm` given just intensity matrices instead of a fitted model.

* Auto-generated initial values set to a small positive number rather than zero when there are no observed data for a particular permitted transition.  This fixes consistency checks (e.g. for `qconstraint`) in this situation.

* Fix when covariates come into the data as one-column matrices instead of vectors.

* Modernised to use roxygen, and pkgdown website created at https://chjackson.github.io/msm

* Some internal functions (e.g. `deriv.msm`) that clash with base S3 generics renamed


# Version 1.7  (2022-11-27)

* `rpexp` is now more efficient.  Thanks to Mark Clements.  Note that the values simulated by `rpexp`, `sim.msm` and `simmulti.msm` for a specific seed will now be different in models where the intensities are piecewise constant.

* Bug fix: in HMMs with partially known initial states specified through `obstrue` and `censor`, `initprobs` is now accounted for.

* Bug fix for bootstrapping with character subject IDs.

* `crudeinits.msm` handled when state is a factor.

* Fix for `hranges` on multiple parameters.

* `obstrue` handled in forward-backward algorithm (`viterbi.msm`).


# Version 1.6.9 

* ... and earlier versions: see [inst/NEWS](https://github.com/chjackson/msm/blob/master/inst/NEWS) in the source for changes.
