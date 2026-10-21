# Sandwich robustness values for applied regression

## Start from the reported finding

The three sandwich measures ask how weak an added scalar index can be
while removing significance under the maintained sandwich test. Its
strength is `max(x,y)`, where x and y are the target-side and
outcome-side partial R-squared. The same index generates the coefficient
and SE change.

``` r

set.seed(1)
n <- 80
d <- rnorm(n)
x <- rnorm(n)
w <- rnorm(n)
y <- .8 * d + .4 * x + (1 + abs(d)) * rnorm(n)
fit <- lm(y ~ d + x + w)
rv <- sandwich_rv(fit, "d", benchmarks = c("x", "w"), vcov = "HC1",
                   n_starts = 4, max_iterations = 30, search_steps = 6)
rv
#> Robustness values for loss of significance (partial R-squared scale)
#>         measure estimate     lower    upper            status benchmark
#>           RV_CH 0.371171 0.3711710 0.371171        ANALYTICAL      <NA>
#>         RV_sand       NA 0.0501088 0.123037  CERTIFIED_BOUNDS      <NA>
#>    RV_cal_delta 0.342267 0.3422672 0.342267 POINT_CERTIFICATE  positive
#>  RV_cal_profile 0.341822 0.3418223 0.341822 POINT_CERTIFICATE  positive
#>  critical
#>    1.9921
#>    1.9921
#>    1.9921
#>    1.9921
#> Profile multiplier k = 1 
#> Bounds are computational intervals; see status for unavailable or unresolved values.
```

The CH row uses the same cutoff as the sandwich test. The conventional
reference with its own default cutoff is also stored in
`rv$RV_CH_standard`. To reproduce a reported test exactly, supply its
fixed `critical` value. HC defaults use augmented residual degrees of
freedom; CR defaults use the number of clusters minus one. Observation
order and rows stay fixed.

## Choose credible observed-control benchmarks

Increment calibration caps the ratio omega = sandwich SE / conventional
SE at `omega_full + max(0,max(delta_omega))`, where each delta is the
signed change when a deleted control is restored. The reduced models
calibrate the cap; the omitted-variable search uses the original full
model.

Profile calibration caps omega at
`a0*exp(k*min(L,beta*(sqrt(x)+sqrt(y))))`. It adjusts restoring log
increments to a common finite-sample correction. L is their positive
maximum and beta is the largest positive increment per `sqrt(x)+sqrt(y)`
in the reduced model. The augmented zero-strength origin a0 uses the
candidate model’s correction. The reference multiplier is k=1; k=2 is a
sensitivity scenario.

``` r

rv$benchmark$records
#>   benchmark rank omega_full omega_minus omega_minus_common_chi delta_omega
#> 1         x    1   1.143028    1.085909               1.085909  0.05711896
#> 2         w    1   1.143028    1.205470               1.205470 -0.06244190
#>   log_increment       r2_d       r2_y    slope    status profile_valid error
#> 1    0.05126343 0.08607993 0.02674409 0.112191 AVAILABLE          TRUE      
#> 2   -0.05318858 0.02776279 0.04691219 0.000000 AVAILABLE          TRUE
c(M_delta = rv$benchmark$M_delta,
  L = rv$benchmark$L, beta = rv$benchmark$beta)
#>    M_delta          L       beta 
#> 1.20014673 0.05126343 0.11219099
profile_k2 <- rv_cal(fit, "d", "profile", c("x", "w"), k = 2,
                      n_starts = 4, max_iterations = 30, search_steps = 6)
profile_k2
#> RV_cal_profile 
#>   Covariance: HC1 
#>   Baseline decision: reject 
#>   Status: POINT_CERTIFICATE / profile 
#>   Benchmark: positive 
#>   Certificate: verified_interval_within_tolerance 
#>   Certified lower bound: 0.3323726 
#>   Certified upper bound: 0.3323726 
#>   Reported RV: 0.3323726 
#>   Witness omega / robust t: 1.266437 / 1.992102 
#>   Runtime: 0.011 seconds
```

Both calibrated classes permit 0 \< omega \<= cap, including omega below
one. Neither calibration generally contains the other. Observed controls
inform an explicit similarity assumption about omitted variables.

Negative increments remain in the records. If the complete scalar pool
has no positive log increment, L=beta=0 and the profile becomes the
constant cap a0 for all k. Empty or incomplete pools are unavailable,
not zero benchmarks.

A character vector lists separate scalar controls. A named list
specifies whole blocks, for example `list(village = "factor(village)")`.
Increment calibration uses the full block change without rank division.
Use separate scalar and block pools. Profile calibration requires scalar
units.

## Interpret and verify the computational bounds

A small verified upper bound establishes that a weak index can remove
significance. A large analytical lower bound establishes robustness
within the stated class. If the bounds do not resolve the comparison,
retain the interval. A missing upper endpoint means the search found no
verified index; it is not a proof that the class is empty. These
intervals describe numerical knowledge of a fixed-sample minimum and are
not confidence intervals.

``` r

value <- rv$results$RV_cal_profile
cache <- cache_from_lm(fit, "d")
check <- evaluate_confounder_reference(cache, value$best_witness_z, vcov = "HC1")
c(x = check$r2_d, y = check$r2_y, omega = check$omega, t = check$t_robust)
#>         x         y     omega         t 
#> 0.3418223 0.3418223 1.2031512 1.9921022
value$control
#> $n_starts
#> [1] 4
#> 
#> $max_iterations
#> [1] 30
#> 
#> $seed
#> [1] 1234
#> 
#> $tolerance
#> [1] 5e-08
#> 
#> $max_radius
#> [1] 0.999999
#> 
#> $search_steps
#> [1] 6
#> 
#> $run_search
#> [1] TRUE
```

Every finite upper bound is the actual strength of an independently
refitted index. `rv` is missing unless lower and upper meet within the
requested tolerance. Boundary attainment and literal strict loss have
separate flags. If the observed test already does not reject, loss RV is
zero by convention.

The lower bound uses an analytical continuous score envelope and the
protected CH formula. Its derivation is installed in `LOWER-BOUND.txt`.
Failed local searches never increase it. Different solvers or stronger
analytical bounds may give different widths for the same estimand.
Ordinary floating-point certificates are not formal interval arithmetic.

## Supported models and reproducibility

The current R interfaces support unweighted, full-rank lm models with an
intercept, a scalar target, HC0/HC1, one-way or raw two-way CR0/CR1, and
fixed-kernel HAC. Cluster vectors must align with rows retained by lm. A
cluster formula is evaluated in the stored model frame. Weighted fits,
adaptive HAC tuning and leverage-dependent corrections require further
implementation.

Set a local `seed` and retain the result object to preserve settings,
benchmark records, witnesses and refits. The caller’s RNG state is
unchanged. Julia is not needed for installation or analysis. Older
[`robust_rv()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/robust_rv.md),
[`rv_cal_bounded()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/rv_cal_bounded.md)
and
[`conditional_omega_rv()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/conditional_omega_rv.md)
preserve their original point/symmetric-band objects; the current
manuscript’s entry points are
[`sandwich_rv()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/sandwich_rv.md),
[`rv_sand()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/rv_sand.md)
and
[`rv_cal()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/rv_cal.md).

## Fixed HAC and raw two-way rules

``` r

hac <- rv_cal(fit, "d", "profile", c("x", "w"),
                vcov = "HAC_Bartlett_L4", n_starts = 4,
                max_iterations = 30, search_steps = 3)
hac
#> RV_cal_profile 
#>   Covariance: HAC_Bartlett_L4 
#>   Baseline decision: reject 
#>   Status: POINT_CERTIFICATE / profile 
#>   Benchmark: positive 
#>   Certificate: verified_interval_within_tolerance 
#>   Certified lower bound: 0.3476631 
#>   Certified upper bound: 0.3476631 
#>   Reported RV: 0.3476631 
#>   Witness omega / robust t: 1.163531 / 1.992102 
#>   Runtime: 0.008 seconds
fixed_hac_rule(nobs(fit), "QS", bandwidth = 5)
#> Fixed sandwich covariance rule: HAC_QS_BW5 
#>   Observations: 80 
#>   PSD repair: none
fixed_twoway_cluster_rule(rep(1:10, each = 8), rep(1:8, 10))
#> Fixed sandwich covariance rule: TWCR1_raw 
#>   Observations: 80 
#>   PSD repair: none
```

`HAC_Bartlett_L4` and `HAC_QS_BW5` match the article: fixed retained-row
order, fixed weights, no prewhitening and no adjustment. The generic HAC
constructor also allows an explicit row permutation and an n/df
adjustment. All choices remain fixed during candidate and deletion
refits.

For `TWCR1_raw`, pass `cluster = data.frame(first, second)` or a
two-variable formula with labels aligned to fitted rows. A rule from
`fixed_twoway_cluster_rule(first, second)` can be passed directly as
`vcov`. The three component factors are G1/(G1-1), G2/(G2-1), and
G12/(G12-1), with a common (n-1)/df factor. The default cutoff df is
min(G1,G2)-1.

No PSD repair is applied. Intermediate negative signed forms are allowed
in direction optimization and positive-target bracketing; every reported
index must have a strictly positive refitted target variance. Invalid
observed variances stop analysis, and invalid deletion variances make
the prescribed calibration unavailable. The indefinite lower envelope is
derived separately in the installed `LOWER-BOUND.txt`.
