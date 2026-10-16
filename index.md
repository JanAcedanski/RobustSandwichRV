# RobustSandwichRV

Pure R sensitivity analysis for linear regressions using HC0/HC1, one-
or two-way CR0/CR1, and fixed-kernel HAC standard errors. Version 0.5.0
computes the current manuscript’s **RV_sand**, **RV_cal^Delta**, and
**RV_cal^prof(k)** on the Cinelli-Hazlett partial R-squared scale. Julia
is not required.

## Install version 0.5.0 from GitHub

Run these commands in R (R 4.2 or newer):

``` r

install.packages("remotes")  # only needed once
remotes::install_github("JanAcedanski/RobustSandwichRV@v0.5.0")
```

The [0.5.0
release](https://github.com/JanAcedanski/RobustSandwichRV/releases/tag/v0.5.0)
also includes the source archive `RobustSandwichRV_0.5.0.tar.gz`.

## Quick start: one table from your regression

This complete example uses the built-in `mtcars` data. Here `wt` is the
coefficient being assessed, and `hp` and `disp` are observed controls
used to calibrate the two `RV_cal` measures.

``` r

library(RobustSandwichRV)
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)

rv <- sandwich_rv(
  fit, treatment = "wt", vcov = "HC1",
  benchmarks = c("hp", "disp"), k = 1, alpha = 0.05
)
rv
as.data.frame(rv)
```

For your own application, replace the regression, the target coefficient
name, and the benchmark control names. `vcov = "HC1"` selects
heteroskedasticity-robust standard errors; the HAC and clustering
options are illustrated below. The same example is installed with the
package:

``` r

demo("quick-start", package = "RobustSandwichRV", ask = FALSE)
```

The four rows are CH, unrestricted sandwich, increment calibrated, and
profile calibrated robustness. CH uses the same test cutoff as the
sandwich rows; the conventional default-cutoff reference is retained in
`rv$RV_CH_standard`. The profile reference is `k = 1`; rerun with
`k = 2` to assess that assumption.

`estimate` is populated only when the analytical lower and independently
verified upper bounds meet within the declared numerical tolerance.
Otherwise, report `[lower, upper]`. An `NA` upper bound means no
admissible loss index was found. A failed local search does not
establish that the candidate class is empty. These are computational
intervals, not confidence intervals.

Values use a 0–1 partial-R-squared scale: `0.10` means 10%, not 0.10%.
For each row, the RV is the minimum common upper limit on the omitted
index’s two partial R-squared strengths needed to reach loss of
significance, under that row’s assumptions. Smaller values indicate
greater sensitivity. The two calibrated rows additionally impose the
observed-control restrictions described below. If the original test does
not reject, its loss RV is zero.

In this example, `benchmark = "zero"` indicates that neither control
supplies a positive restoring increment. Here the two calibrated values
coincide; inspect `rv$benchmark$records` for your own data.

## Individual calculations

``` r

rv_sand(fit, "wt", vcov = "HC1")
rv_cal(fit, "wt", "increment", benchmarks = c("hp", "disp"))
rv_cal(fit, "wt", "profile", benchmarks = c("hp", "disp"), k = 1)

# User-specified ratio cap: 0 < omega <= 1.2
rv_cal(fit, "wt", "cap", M = 1.2)
```

For one-way clustering, use `vcov = "CR1"` and cluster labels aligned
with the rows retained by `lm`. A cluster formula works when its
variable is present in the stored model frame. Supply `critical = ...`
to reproduce an application’s reported cutoff exactly; the chosen cutoff
stays fixed during the search.

## HAC and two-way clustering

The article’s HAC settings have direct aliases:

``` r

hac_bartlett <- sandwich_rv(fit, "wt", c("hp", "disp"),
                             vcov = "HAC_Bartlett_L4")
hac_qs <- sandwich_rv(fit, "wt", c("hp", "disp"),
                       vcov = "HAC_QS_BW5")
```

These use the retained row order, fixed Bartlett lag 4 or
quadratic-spectral bandwidth 5, no prewhitening, and no finite-sample
adjustment. Other fixed settings and an explicit row permutation can be
supplied as a rule:

``` r

hac_rule <- fixed_hac_rule(nobs(fit), kernel = "QS", bandwidth = 5,
                           order = seq_len(nobs(fit)), adjust = FALSE)
rv_cal(fit, "wt", "profile", c("hp", "disp"), vcov = hac_rule)
```

For two-way clustering, supply both label vectors in the fitted sample’s
order:

``` r

# firm_retained and year_retained must match the rows retained by fit
rv_two <- sandwich_rv(fit, "wt", c("hp", "disp"),
                       vcov = "TWCR1_raw",
                       cluster = data.frame(firm_retained, year_retained))
# A formula such as cluster = ~ firm + year is also accepted.
```

The raw operator adds the two one-way terms and subtracts their
intersection. CR1 applies a separate G/(G-1) multiplier to each term and
the common (n-1)/residual-df factor. No PSD repair or HC0 intersection
replacement is applied. The default cutoff df is min(G1,G2)-1.
`TWCR0_raw` is also available. Alternatively create a reusable
`fixed_twoway_cluster_rule(firm, year)`.

Nonpositive observed target variances stop the calculation. Nonpositive
candidate variances are inadmissible; an invalid prescribed deletion
makes calibration unavailable. Signed quadratic optimization and path
roots preserve this distinction. The full frozen rule is stored in
`value$covariance_rule`:

``` r

value <- rv_two$results$RV_cal_profile
cache <- cache_from_lm(fit, "wt")
evaluate_confounder_reference(cache, value$best_witness_z,
                               vcov = value$covariance_rule)
```

## What the calibrations assume

Let omega be the target’s sandwich SE divided by its conventional SE.
Increment calibration deletes each observed control, keeps the same
sample, then records the signed restoring change in the target’s ratio.
Its cap is `omega_full + max(0, max(delta_omega))`.

Profile calibration adjusts restoring log increments to a common
finite-sample correction. Its cap at strengths x,y is
`a0 * exp(k * min(L, beta * (sqrt(x) + sqrt(y))))`, where `a0` is the
zero-strength augmented-model origin. `L` and `beta` are calibrated from
scalar controls in their reduced models. Both classes allow every
positive omega below the cap, including values below one. They represent
explicit similarity assumptions about possible omitted variables.

``` r

b <- rv_benchmarks(fit, "wt", c("hp", "disp"))
b$records
b$M_delta
c(L = b$L, beta = b$beta)
```

Negative increments are retained. A complete zero benchmark is a
distinct, restrictive scenario. Empty or incomplete pools return
`CALIBRATION_UNAVAILABLE`; invalid refits are never silently discarded.

A character vector lists **separate** benchmark units. A named list can
define whole blocks, for example `list(village = "factor(village)")`.
The increment uses each whole block’s change without dividing by rank.
Submit scalar and block pools separately. Profile calibration requires
scalar units. By default, `benchmarks = NULL` selects eligible scalar
model terms, excluding the target and intercept; inspect the returned
records before interpreting the calibration.

## Verify and reproduce a result

``` r

value <- rv_cal(fit, "wt", "profile", c("hp", "disp"), seed = 1234)
value$control
value$benchmark$records
value$witness_refit

cache <- cache_from_lm(fit, "wt")
evaluate_confounder_reference(cache, value$best_witness_z, vcov = "HC1")
```

Every finite new upper bound uses the actual refitted strength of a
synthetic index. `witness_literal_decision_change` distinguishes strict
loss from reaching the closed boundary. Numerical point certification
uses ordinary floating point. If the observed test already does not
reject, its loss RV is zero by convention.

The current R release supports unweighted, full-rank `lm` models with an
intercept and one scalar target. HC2/HC3, CR2, adaptive HAC tuning,
weighted fits, and literal multivariate omitted blocks are outside the
new entry points. The unrestricted analytical envelope is documented in
`inst/LOWER-BOUND.txt`; it can be less sharp than specialized research
bounds in the article.

## Earlier interfaces

[`robust_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/robust_rv.md),
[`rv_cal_bounded()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_cal_bounded.md)
and
[`conditional_omega_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/conditional_omega_rv.md)
retain their earlier point/symmetric-band definitions and HC/one-way CR
scope for compatibility. Use
[`sandwich_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/sandwich_rv.md),
[`rv_sand()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_sand.md)
and
[`rv_cal()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_cal.md)
for the current manuscript’s estimands.

## Development checks

``` sh
Rscript tools/check-local.R .
R CMD build .
R CMD check --no-manual RobustSandwichRV_0.5.0.tar.gz
```

The tests cover independent sandwich refits, both calibration formulas,
score gradients, analytical envelopes, missing and zero pools, preserved
random state, and frozen Julia covariance fixtures. Full article replay
scripts in `tools/` add comparisons against stored current-manuscript
matrices and indices.
