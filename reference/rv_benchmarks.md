# Restoring-increment and profile benchmarks from observed controls

Deletes each prescribed control unit and refits the same target on the
same rows. Records signed changes in the target's
sandwich-to-conventional SE ratio and explanatory strengths in the
reduced model.

## Usage

``` r
rv_benchmarks(model, treatment, benchmarks = NULL, vcov = "HC1", cluster = NULL)
```

## Arguments

- model:

  An unweighted full-rank `lm` model with an intercept.

- treatment:

  Target model-matrix column name.

- benchmarks:

  Character vector of separate controls, or a named list of whole
  blocks. NULL selects eligible scalar terms.

- vcov:

  HC0/HC1, one-way CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw,
  TWCR1_raw, or a fixed covariance rule.

- cluster:

  One- or two-way cluster labels aligned with the fitted sample, or a
  formula.

## Details

For increment calibration, M equals omega in the full model plus the
largest positive restoring increment, or zero added increment if all are
nonpositive. Whole blocks are deleted together and their changes are not
divided by rank. Scalar and block pools must be submitted separately.

For profile calibration, the reduced ratio is adjusted to the full
model's chi correction, chi(v) = v\*a(v). L is the positive maximum log
increment; beta is the maximum positive log increment divided by
sqrt(x)+sqrt(y), measured in the reduced model. The profile requires
scalar units. A zero denominator and zero increment contributes zero;
nonzero change at zero strength makes the profile unavailable.

Empty and incomplete pools are unavailable. Negative changes remain in
the records. Variances must be positive and finite; a reduced SE ratio
at or below 100 times machine epsilon is treated as numerical zero.
Target and intercept deletion, overlapping units and mixed scalar/block
pools are rejected. Deletions use stored model rows and matrices,
avoiding changes to missing-value selection.

## Value

An `rv_calibration` with `records`, pool status, `M_delta`, `L`, `beta`,
zero-benchmark indicators through the status and parameter values, and
binding control names.

## See also

[`rv_cal`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/rv_cal.md)

## Examples

``` r
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
rv_benchmarks(fit, "wt", c("hp", "disp"))
#> $records
#>   benchmark rank omega_full omega_minus omega_minus_common_chi delta_omega
#> 1        hp    1   0.908662   0.9466929              0.9466929 -0.03803099
#> 2      disp    1   0.908662   1.0291912              1.0291912 -0.12052926
#>   log_increment       r2_d         r2_y slope    status profile_valid error
#> 1   -0.04100166 0.02400574 0.2095480108     0 AVAILABLE          TRUE      
#> 2   -0.12455541 0.63534282 0.0002926464     0 AVAILABLE          TRUE      
#> 
#> $status
#> [1] "COMPLETE_ZERO_INCREMENT"
#> 
#> $complete
#> [1] TRUE
#> 
#> $scalar
#> [1] TRUE
#> 
#> $profile_available
#> [1] TRUE
#> 
#> $omega_full
#> [1] 0.908662
#> 
#> $M_delta
#> [1] 0.908662
#> 
#> $delta_positive
#> [1] 0
#> 
#> $L
#> [1] 0
#> 
#> $beta
#> [1] 0
#> 
#> $binding_delta
#> character(0)
#> 
#> $binding_level
#> character(0)
#> 
#> $binding_slope
#> [1] "hp"   "disp"
#> 
#> $treatment
#> [1] "wt"
#> 
#> $vcov
#> [1] "HC1"
#> 
#> $covariance_rule
#> Fixed sandwich covariance rule: HC1 
#>   Observations: 32 
#>   PSD repair: none
#> 
#> attr(,"class")
#> [1] "rv_calibration"
```
