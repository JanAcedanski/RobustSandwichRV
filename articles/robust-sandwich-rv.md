# Certified Robustness Values with Sandwich Standard Errors

## The object being computed

For a scalar omitted confounder, let

``` math
x=R^2_{D\sim Z\mid X},\qquad
y=R^2_{Y\sim Z\mid D,X},\qquad
\omega_R(Z)=\frac{se_R(\widehat\tau_f)}{se_C(\widehat\tau_f)}.
```

The bounded calibrated value uses an observed-control magnitude $`b`$
and the zero-strength synthetic-index origin $`\omega_{R,00}`$:

``` math
\mathcal W_R(b)=
[\omega_{R,00}e^{-b},\omega_{R,00}e^b].
```

It finds the smallest $`\max(x,y)`$ among attainable confounders that
change the robust test decision and have
$`\omega_R(Z)\in\mathcal W_R(b)`$.

The software first computes an analytical protected lower bound at the
adverse band endpoint. It then either certifies equality constructively
or returns the best independently refitted witness as a certified upper
bound.

The result records `theorem_exact` separately from `numerically_closed`.
The former requires the adverse-endpoint equality certificate; the
latter only means that a verified lower–upper interval is no wider than
the stated numerical tolerance.

## Reproducible HC1 example

``` r

set.seed(20260809)
n <- 120
x1 <- rnorm(n)
x2 <- rnorm(n)
d <- 0.4 * x1 - 0.2 * x2 + rnorm(n)
e <- (0.7 + 0.3 * abs(d)) * rnorm(n)
y <- 0.55 * d + 0.25 * x1 - 0.15 * x2 + e
fit <- lm(y ~ d + x1 + x2)
```

We use `x1` as an observed scalar benchmark.

``` r

benchmark <- benchmark_from_model(fit, "d", "x1", vcov = "HC1")
benchmark
#> Observed-control benchmark calibration
#>   Name: x1 
#>   Type / rank: scalar / 1 
#>   Signed per-rank delta_L: 0.03733079 
#>   Bounded magnitude b: 0.03733079
```

``` r

result <- robust_rv(
  fit,
  treatment = "d",
  vcov = "HC1",
  benchmark = benchmark,
  n_starts = 8,
  max_iterations = 20
)
result
#> Benchmark-bounded robust RV 
#>   Covariance: HC1 
#>   Baseline decision: reject 
#>   Status: EXACT_BOUNDED_RV / ADVERSE_ENDPOINT_PROTECTED_DIAGONAL 
#>   Certificate: adverse_endpoint_equality 
#>   Certified lower bound: 0.2092966 
#>   Certified upper bound: 0.2092966 
#>   Reported RV: 0.2092966 
#>   Witness omega / robust t: 1.400975 / 1.980806 
#>   Runtime: 0.119 seconds
```

`RV_lower` is analytical. A finite `RV_upper` always has an explicit
`best_witness_z` and `direct_refit_valid = TRUE`. `LOWER_BOUND_ONLY`
means that the numerical search did not find an upper witness; it is not
a global infeasibility statement.

## CR1 uses the same algorithm

``` r

cluster <- rep(seq_len(12), each = 10)
cluster_benchmark <- benchmark_from_model(
  fit, "d", "x1", vcov = "CR1", cluster = cluster
)
cluster_result <- robust_rv(
  fit,
  treatment = "d",
  vcov = "CR1",
  cluster = cluster,
  benchmark = cluster_benchmark,
  n_starts = 8,
  max_iterations = 20
)
cluster_result
#> Benchmark-bounded robust RV 
#>   Covariance: CR1 
#>   Baseline decision: reject 
#>   Status: EXACT_BOUNDED_RV / ADVERSE_ENDPOINT_PROTECTED_DIAGONAL 
#>   Certificate: adverse_endpoint_equality 
#>   Certified lower bound: 0.1390936 
#>   Certified upper bound: 0.1390936 
#>   Reported RV: 0.1390936 
#>   Witness omega / robust t: 1.677225 / 2.200985 
#>   Runtime: 0.047 seconds
```

Only the score aggregation changes: HC sums squared observation scores,
while CR first sums scores within clusters. Protected formulas,
bracketing, spherical root finding, fallback logic, and direct
validation are shared.

## Point conditioning

``` r

point_result <- conditional_omega_rv(
  cache_from_lm(fit, "d"),
  omega_target = 1,
  vcov = "HC1",
  method = "protected_first",
  n_starts = 8,
  max_iterations = 20
)
point_result
#> Point-conditioned robust RV 
#>   Covariance: HC1 
#>   Baseline decision: reject 
#>   Status: EXACT_POINT_RV / PROTECTED_DIAGONAL 
#>   Certified lower bound: 0.265243 
#>   Certified upper bound: 0.265243 
#>   Reported RV: 0.265243 
#>   Witness omega / robust t: 1 / 1.980808 
#>   Runtime: 0.045 seconds
```

Point conditioning is useful for theory and diagnostics, but it imposes
one exact signed sandwich ratio. Bounded calibration is the recommended
observed-control workflow.

## Refit the witness explicitly

``` r

if (!is.null(result$best_witness_z)) {
  witness <- evaluate_confounder_reference(
    cache_from_lm(fit, "d"),
    result$best_witness_z,
    vcov = "HC1"
  )
  c(
    x = witness$r2_d,
    y = witness$r2_y,
    omega = witness$omega,
    robust_t = witness$t_robust
  )
}
#>         x         y     omega  robust_t 
#> 0.2092968 0.2092968 1.4009752 1.9808056
```

## Scope

Version 0.1.0 supports unweighted linear regression and one-way
clustering. HC2, HC3, CR2, HAC, weighted least squares, generalized
linear models, and critical values that change with the omitted
confounder are outside the current certificate theory.
