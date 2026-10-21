# Changelog

## RobustSandwichRV 0.5.0

- Added version-pinned GitHub installation instructions and a runnable
  `demo("quick-start")` with a practitioner-oriented explanation of the
  table.
- Added fixed Bartlett and quadratic-spectral HAC, including article
  aliases `HAC_Bartlett_L4` and `HAC_QS_BW5`, frozen row ordering and
  explicit bandwidths.
- Added raw two-way CR0/CR1 via `TWCR0_raw`, `TWCR1_raw` and a rule
  constructor. CR1 preserves separate component corrections and the
  actual residual df.
- Extended all three current RV estimands, benchmark deletion refits,
  independent covariance evaluation and
  [`validate_sandwich_reference()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/validate_sandwich_reference.md)
  to these rules.
- Added a continuous signed-form lower envelope for indefinite operators
  and signed-quadratic direction optimization/root finding. Nonpositive
  target variances are rejected, with no PSD repair.
- Full current-article validation now covers all 112 cases and all
  estimators.

## RobustSandwichRV 0.2.0

- Added
  [`rv_sand()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/rv_sand.md)
  for unrestricted scalar-index loss robustness, with an analytical
  continuous-envelope lower bound and independently refitted uppers.
- Added
  [`rv_cal()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/rv_cal.md)
  with restoring-increment and strength-dependent profile calibration
  from the current manuscript, allowing every positive ratio below the
  cap. Direct fixed caps are also supported.
- Added
  [`rv_benchmarks()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/rv_benchmarks.md)
  retaining all signed deletion increments, common-chi corrections,
  partial R-squared values and binding controls. Empty, incomplete and
  complete-zero pools are distinguished. Blocks are not rank-normalized
  in the new increment calibration; the profile requires scalar
  benchmarks.
- Added
  [`sandwich_rv()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/sandwich_rv.md)
  and a four-row practitioner table on a common test cutoff.
- Retained the legacy point and symmetric-band APIs for compatibility.

## RobustSandwichRV 0.1.0

- Initial pure-R implementation.
- Added classical, HC0, HC1, CR0, and CR1 scalar-confounder evaluation.
- Added protected-first, diagonal-first, and full point-conditioned
  workflows.
- Added benchmark-bounded calibration with certified lower and upper
  bounds.
- Added scalar and rank-normalized grouped observed-control benchmarks.
- Added explicit full-regression witness validation and nested-band
  reuse.
- Added two-direction reduction for arbitrary omitted-confounder blocks.
- Added frozen cross-language numerical fixtures, `sensemakr`
  comparison, and `sandwich` covariance validation.
