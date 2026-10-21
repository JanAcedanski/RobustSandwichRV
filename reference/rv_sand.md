# Unrestricted sandwich robustness value for loss of significance

Computes lower and upper bounds on the weakest added scalar index that
removes significance under the specified sandwich rule. Strength is the
maximum of treatment-side and outcome-side partial R-squared.

## Usage

``` r
rv_sand(model, treatment = NULL, vcov = "HC1", cluster = NULL,
  alpha = 0.05, critical = NULL, null = 0, n_starts = 8L,
  max_iterations = 80L, seed = 1234L, tolerance = 5e-8,
  max_radius = 1 - 1e-6, search_steps = 12L, run_search = TRUE,
  supplied_witnesses = list())
```

## Arguments

- model:

  An unweighted full-rank `lm` model with an intercept, or an
  `rv_cache`.

- treatment:

  Name of one target model-matrix column; required for `lm`.

- vcov:

  HC0/HC1, one-way CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw,
  TWCR1_raw, or a fixed covariance rule.

- cluster:

  One- or two-way cluster labels aligned with rows retained by the fit,
  or a formula evaluated in its model frame.

- alpha:

  Two-sided significance level.

- critical:

  Optional fixed positive cutoff. Defaults to a t cutoff with augmented
  residual degrees of freedom for HC/HAC, number of clusters minus one
  for one-way CR, and min(G1,G2)-1 for two-way CR.

- null:

  Null value of the target coefficient.

- n_starts:

  Number of screened residual directions.

- max_iterations:

  Maximum local optimizer iterations.

- seed:

  Local random seed, leaving the caller's RNG state unchanged.

- tolerance:

  Positive numerical refit and interval-closing tolerance, at most
  0.001.

- max_radius:

  Largest strength queried by the radial search; strictly between zero
  and one. A previously verified witness above this value may still be
  retained.

- search_steps:

  Number of radial queries.

- run_search:

  Run fallback search after direct checks.

- supplied_witnesses:

  Vectors or result objects; each candidate is independently rechecked.

## Details

An analytical continuous score envelope supplies the lower bound. Failed
local searches never increase this bound. Every finite upper bound is
the measured strength of an actual index verified by a full augmented QR
regression. A point value is returned only when the verified bounds
agree within tolerance. The bounds are computational intervals, not
confidence intervals or formal interval arithmetic.

The observed test and cutoff are held fixed. If the observed test
already does not reject (absolute t no greater than the cutoff), the
loss RV is zero by convention. Boundary attainment and strict loss are
recorded separately. The current implementation supports unweighted OLS
with HC0/HC1, one-way CR0/CR1, raw two-way CR0/CR1, and fixed Bartlett
or quadratic-spectral HAC. Weights and order stay fixed; adaptive tuning
and PSD repair are not performed. Nonpositive candidate variances are
inadmissible. See
[`fixed_hac_rule`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/fixed_score_rules.md)
for rule constructors and article aliases.

## Value

A `robust_rv_result`: `RV_lower`, `RV_upper`, `rv` (missing unless
resolved), status, tolerances, runtime, and the independently refitted
`best_witness_z`. A missing upper bound means that no admissible loss
index was found; it does not establish infeasibility.

## See also

[`rv_cal`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/rv_cal.md),
[`sandwich_rv`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/sandwich_rv.md)

## Examples

``` r
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
rv_sand(fit, "wt", run_search = FALSE)
#> RV_sand 
#>   Covariance: HC1 
#>   Baseline decision: reject 
#>   Status: LOWER_BOUND_ONLY / sand 
#>   Certificate: analytical_lower_bound_only 
#>   Certified lower bound: 0.02032113 
#>   Certified upper bound: not found 
#>   Runtime: 0.002 seconds
```
