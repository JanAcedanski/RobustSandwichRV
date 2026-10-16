# Fixed HAC and raw two-way cluster covariance rules

Construct fixed score operators for use as the `vcov` argument in
`sandwich_rv`, `rv_sand`, `rv_cal`, `rv_benchmarks` and the confounder
evaluators.

## Usage

``` r
fixed_hac_rule(n, kernel = "bartlett", bandwidth = 4,
  order = seq_len(n), adjust = FALSE, critical_df = NULL)
fixed_twoway_cluster_rule(cluster1, cluster2,
  variant = c("CR1", "CR0"), critical_df = NULL)
# S3 method for class 'rv_score_rule'
print(x, ...)
```

## Arguments

- n:

  Number of rows retained by the regression.

- kernel:

  `"bartlett"` or `"quadratic_spectral"`; `"QS"` is also accepted.

- bandwidth:

  For Bartlett, a positive integer maximum lag smaller than n; for QS, a
  fixed positive bandwidth.

- order:

  Permutation giving the temporal order of the retained rows.

- adjust:

  If TRUE, use an n/residual-df correction. FALSE matches the article's
  unadjusted HAC rules.

- critical_df:

  Optional positive cutoff degrees of freedom or function of residual
  df. Defaults to augmented residual df for HAC and min(G1,G2)-1 for
  two-way clustering. The resulting test cutoff is held fixed during RV
  search.

- cluster1,cluster2:

  Two cluster dimensions aligned with retained observations, with at
  least two groups in each.

- variant:

  CR0 or CR1.

- x:

  A fixed `rv_score_rule`.

- ...:

  Unused.

## Details

HAC freezes row order, bandwidth and kernel weights. Bartlett maximum
lag L uses weight 1-h/(L+1) at h\<=L. QS uses 3/z^2\*(sin(z)/z-cos(z)),
z=6\*pi\*h/(5\*bandwidth), at all retained lags. No prewhitening or
automatic bandwidth selection is performed. The aliases
`"HAC_Bartlett_L4"` and `"HAC_QS_BW5"` reproduce the article's settings
in current fitted-row order.

Raw two-way clustering uses A1+A2-A12. For CR1, each component has its
own Gj/(Gj-1) multiplier in the fixed operator, and the common
(n-1)/residual-df multiplier changes with model rank. The intersection
term is not replaced by HC0. There is no PSD repair. A nonpositive
observed target variance is an error; a nonpositive candidate variance
is inadmissible; an invalid deletion makes the prescribed calibration
pool unavailable.

Alternatively pass `vcov="TWCR1_raw"` or `"TWCR0_raw"` with a two-column
cluster data frame, matrix, or formula. A rule object embeds its order
or cluster partitions, so subsequent cluster arguments are unnecessary.
Rules are reused unchanged for every augmented and reduced model; only
the specified scalar rank correction changes.

## Value

An `rv_score_rule` containing fixed weights or partitions, component
corrections, critical-df convention and metadata. The complete rule is
also stored as `result$covariance_rule` for independent refits.

## See also

[`sandwich_rv`](https://JanAcedanski.github.io/RobustSandwichRV/reference/sandwich_rv.md),
[`validate_sandwich_reference`](https://JanAcedanski.github.io/RobustSandwichRV/reference/validate_sandwich_reference.md)

## Examples

``` r
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
hac <- fixed_hac_rule(nobs(fit), "QS", bandwidth = 5)
rv_cal(fit, "wt", "profile", c("hp", "disp"), vcov = hac,
       n_starts = 4, max_iterations = 20, run_search = FALSE)
#> RV_cal_profile 
#>   Covariance: HAC_QS_BW5 
#>   Baseline decision: reject 
#>   Status: POINT_CERTIFICATE / profile 
#>   Benchmark: positive 
#>   Certificate: verified_interval_within_tolerance 
#>   Certified lower bound: 0.3629549 
#>   Certified upper bound: 0.3629549 
#>   Reported RV: 0.3629549 
#>   Witness omega / robust t: 0.5545122 / -2.051831 
#>   Runtime: 0.013 seconds
fixed_twoway_cluster_rule(rep(1:8, each = 4), rep(1:4, 8))
#> Fixed sandwich covariance rule: TWCR1_raw 
#>   Observations: 32 
#>   PSD repair: none
```
