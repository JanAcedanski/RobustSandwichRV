## Test environment

* Windows 11 x64, R 4.6.1 (ucrt), 16 October 2026.

## Local check result

R CMD check --no-manual RobustSandwichRV_0.5.0.tar.gz

0 errors | 0 warnings | 0 notes

The full check includes installation, examples, tests and rebuilt vignettes.

## Numerical validation

* All 112 current-article cases replayed in R: 37 HC1, 20 one-way CR1,
  24 fixed Bartlett HAC L4, 23 fixed quadratic-spectral HAC BW5, and 8 raw
  two-way CR1 cases. All prescribed scalar-control deletion benchmarks were
  recomputed on the maintained samples. Profile lower bounds and actual
  point witnesses agree with the article within 5e-8; every retained upper
  index was independently refitted in R.
* Fresh R searches without supplied Julia indices: 18 calibrated point
  results and 6 unrestricted RV_sand intervals, covering all estimator
  families and all seven datasets.
* Independent covariance checks against sandwich cover HC0/HC1, one-way
  CR0/CR1, raw two-way CR0/CR1, Bartlett and quadratic-spectral HAC with
  both unadjusted and n/df corrections.
* Tests verify signed quadratic gradients, a negative-to-positive variance
  path root, valid indefinite lower envelopes, fixed row ordering, separate
  component corrections, zero and unavailable benchmark pools, and RNG state.

No PSD repair, automatic HAC tuning or prewhitening is applied. Nonpositive
observed/candidate target variances and invalid deletion variances retain
separate meanings. Ordinary floating-point certificates are not formal
interval arithmetic or claims about finite-sample test size.

Julia is only used by development export tools, not by installation or R
analysis. The unrestricted analytical envelope may be looser than specialized
article bounds; unresolved minima remain computational intervals.
