# tmle3shift 0.2.2

* Change default value of argument `max_shifted_ratio` from 2 to 5 throughout.
* Bump up minimum versions of dependencies `sl3`, `hal9001`, and `haldensify`.

# tmle3shift 0.2.1

* Change all instances of default shift being `shift_additive_bounded` to use
  `shift_additive` instead, since the former is incompatible with the simpler
  identifying functional under the modified treatment policy definition.

# tmle3shift 0.1.9

* Change `Spec_shift` object to use point treatment helpers from `tmle3`.
* Support missing data functionality based on updates to `tmle3`.
* Change conditional density estimation in tests and vignettes to rely on
  `Lrnr_density_semiparametric` over `Lrnr_haldensify` to reduce time.
* Remove reliance on `Lrnr_rfcde` since package appears unstable.

# tmle3shift 0.1.8

* Change default for MSMs to use identity-based weighting over variance-based,
  improving stability of estimates.
* Bug fixes for cross-validation: `cv_fold` -> `fold_number`.

# tmle3shift 0.1.7

* Addition of `cvtmle = TRUE` in `Updater` objects as the default.
* Bug fixes to bounding of conditional densities at the natural exposure value.
* Tweaks to `cvtmle` option and `fold_number` argument to match `tmle3` updates.

# tmle3shift 0.1.6

* Addition of HAL-based conditional density estimation as an option.
* Addition of random forest-based conditional density estimation as an option.
* Addition of bounding to `Spec` objects for variable importance measures.
* Addition of bounding for conditional density and likelihood objects.
* Minor bug fixes, including to vignettes.

# tmle3shift 0.1.5

* An initial stable release of the package.
