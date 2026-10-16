# Compare CH, unrestricted sandwich and both calibrated robustness values

Returns a practitioner table with RV_CH, RV_sand, RV_cal_delta and
RV_cal_profile on the same partial-R-squared scale and test cutoff.

## Usage

``` r
sandwich_rv(model, treatment, benchmarks = NULL, vcov = "HC1",
  cluster = NULL, k = 1, alpha = 0.05, critical = NULL, null = 0, ...)
# S3 method for class 'sandwich_rv_analysis'
print(x, ...)
# S3 method for class 'sandwich_rv_analysis'
summary(object, ...)
# S3 method for class 'sandwich_rv_analysis'
as.data.frame(x, row.names = NULL,
  optional = FALSE, ...)
```

## Arguments

- model:

  Unweighted `lm` with an intercept and full-rank design.

- treatment:

  One target model-matrix column.

- benchmarks:

  Observed control units; NULL selects scalar terms.

- vcov:

  HC0/HC1, one-way CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw,
  TWCR1_raw, or a fixed covariance rule.

- cluster:

  One- or two-way cluster labels aligned with fitted rows, or a formula.

- k:

  Positive profile multiplier; default one.

- alpha:

  Two-sided significance level.

- critical:

  Optional fixed cutoff common to the comparison.

- null:

  Null coefficient value.

- ...:

  Numerical solver arguments from
  [`rv_sand`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_sand.md);
  unused in methods.

- x,object:

  A `sandwich_rv_analysis`.

- row.names,optional:

  Standard data-frame method arguments.

## Details

The calibrated searches are run first so their independently verified
indices can also improve the unrestricted upper bound, after rechecking.
The two calibrated classes are not generally nested. A point estimate is
shown only when bounds close within tolerance. The table's CH row uses
the sandwich cutoff; the conventional default-cutoff reference is
retained separately as `RV_CH_standard`.

## Value

A `sandwich_rv_analysis` with `table`, three full `results`, benchmark
diagnostics, `RV_CH`, `RV_CH_standard`, and k.
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
the four-row table. Its benchmark column labels positive, zero and
unavailable calibrations separately from computational status.

## Examples

``` r
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
sandwich_rv(fit, "wt", c("hp", "disp"), n_starts = 4,
            max_iterations = 20, search_steps = 3)
#> Robustness values for loss of significance (partial R-squared scale)
#>         measure estimate     lower     upper            status benchmark
#>           RV_CH 0.242656 0.2426556 0.2426556        ANALYTICAL      <NA>
#>         RV_sand       NA 0.0203211 0.0825404  CERTIFIED_BOUNDS      <NA>
#>    RV_cal_delta 0.269198 0.2691982 0.2691982 POINT_CERTIFICATE      zero
#>  RV_cal_profile 0.269198 0.2691982 0.2691982 POINT_CERTIFICATE      zero
#>  critical
#>   2.05183
#>   2.05183
#>   2.05183
#>   2.05183
#> Profile multiplier k = 1 
#> Bounds are computational intervals; see status for unavailable or unresolved values.
```
