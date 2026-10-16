# Increment and strength-dependent profile calibrated robustness values

Computes the weakest added scalar index satisfying 0 \< omega \<= the
calibrated cap that removes significance in the maintained full model.

## Usage

``` r
rv_cal(model, treatment = NULL,
  calibration = c("increment", "profile", "cap"), benchmarks = NULL,
  vcov = "HC1", cluster = NULL, k = 1, M = NULL,
  alpha = 0.05, critical = NULL, null = 0, ...)
```

## Arguments

- model:

  An unweighted `lm` model, or an `rv_cache`. Observed-control
  calibration requires stored `lm` metadata.

- treatment:

  Target model-matrix column.

- calibration:

  Restoring-increment, strength-dependent profile, or direct fixed cap.

- benchmarks:

  Separate control units as in
  [`rv_benchmarks`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_benchmarks.md).

- vcov:

  HC0/HC1, one-way CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw,
  TWCR1_raw, or a fixed covariance rule.

- cluster:

  One- or two-way cluster labels aligned with fitted rows, or a formula.

- k:

  Positive profile multiplier. The reference is one; two is a
  sensitivity scenario.

- M:

  Positive fixed cap when calibration is `"cap"`.

- alpha:

  Two-sided test level.

- critical:

  Optional fixed test cutoff; see
  [`rv_sand`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_sand.md)
  for defaults.

- null:

  Null coefficient value.

- ...:

  Named numerical arguments from
  [`rv_sand`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_sand.md),
  including previously found witnesses.

## Details

Increment calibration uses `M_delta` from `rv_benchmarks`. Profile
calibration imposes a0\*exp(k\*min(L,beta\*(sqrt(x)+sqrt(y)))); a0 is
the zero-strength augmented-model ratio, with the actual finite-sample
correction. Both permit ratios below one. A complete zero profile
benchmark gives the constant cap a0 for every k. An empty or invalid
pool is never replaced by a zero benchmark.

The profile lower bound uses the crossing of the protected CH bound and
the envelope a0\*exp(k\*min(L,2\*beta\*sqrt(r))). The unrestricted
analytical envelope can strengthen both calibrated lower bounds.
Endpoint and fallback candidates are tested against their actual
pair-specific cap and independently refitted. Failure to attain the
relaxation returns an interval or lower bound only. These are
sensitivity assumptions about omitted variables, not consequences of the
observed controls.

## Value

A `robust_rv_result`, including the complete benchmark records and
profile parameters when applicable. `CALIBRATION_UNAVAILABLE` has
missing bounds. For an already nonrejecting observed test and an
available calibration, loss RV is zero by convention. See
[`rv_sand`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_sand.md)
for numerical and inferential interpretation.

## See also

[`sandwich_rv`](https://JanAcedanski.github.io/RobustSandwichRV/reference/sandwich_rv.md),
[`rv_benchmarks`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_benchmarks.md)

## Examples

``` r
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
rv_cal(fit, "wt", "increment", c("hp", "disp"), run_search = FALSE)
#> RV_cal_delta 
#>   Covariance: HC1 
#>   Baseline decision: reject 
#>   Status: POINT_CERTIFICATE / increment 
#>   Benchmark: zero 
#>   Certificate: verified_interval_within_tolerance 
#>   Certified lower bound: 0.2691982 
#>   Certified upper bound: 0.2691982 
#>   Reported RV: 0.2691982 
#>   Witness omega / robust t: 0.908662 / -2.051831 
#>   Runtime: 0.014 seconds
rv_cal(fit, "wt", "profile", c("hp", "disp"), k = 1, run_search = FALSE)
#> RV_cal_profile 
#>   Covariance: HC1 
#>   Baseline decision: reject 
#>   Status: POINT_CERTIFICATE / profile 
#>   Benchmark: zero 
#>   Certificate: verified_interval_within_tolerance 
#>   Certified lower bound: 0.2691982 
#>   Certified upper bound: 0.2691982 
#>   Reported RV: 0.2691982 
#>   Witness omega / robust t: 0.908662 / -2.051831 
#>   Runtime: 0.009 seconds
```
