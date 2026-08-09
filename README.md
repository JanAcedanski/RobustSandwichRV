
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# RobustSandwichRV

<!-- badges: start -->
[![R-CMD-check](https://github.com/JanAcedanski/RobustSandwichRV/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JanAcedanski/RobustSandwichRV/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`RobustSandwichRV` computes omitted-variable robustness values for
linear regressions using classical, HC0, HC1, CR0, or CR1 standard
errors. It combines Cinelli–Hazlett partial-$R^2$ geometry with the
sample attainability of the sandwich-to-classical standard-error ratio.

The recommended observed-control workflow is **benchmark-bounded
calibration**. It reports

``` text
analytical protected lower bound <= calibrated RV <= verified witness upper bound
```

and labels equality as exact only when supported by a constructive
certificate from the adverse-endpoint theorem. A separately labelled
numerically closed result means that the verified interval is within a
declared tolerance, not that the theorem supplied equality. Every finite
upper endpoint is an explicit confounder checked in an independent
full-regression refit.

The package is implemented entirely in R. Julia is not needed to
install, test, or use it.

## Installation

``` r
install.packages("remotes")
remotes::install_github("JanAcedanski/RobustSandwichRV")
```

## HC1 example

``` r
library(RobustSandwichRV)

fit <- lm(y ~ treatment + age + female + factor(village), data = dat)

result <- robust_rv(
  fit,
  treatment = "treatment",
  vcov = "HC1",
  calibration = "bounded",
  benchmark = "female"
)

result
summary(result)
```

## One-way CR1 example

``` r
clustered <- robust_rv(
  fit,
  treatment = "treatment",
  vcov = "CR1",
  cluster = dat$village,
  calibration = "bounded",
  benchmark = "female"
)
```

If the cluster variable is present in the fitted model frame, a formula
can be used:

``` r
clustered <- robust_rv(
  fit,
  treatment = "treatment",
  vcov = "CR1",
  cluster = ~ village,
  benchmark = "female"
)
```

## Grouped benchmarks

A factor or a vector of observed terms is treated as a block. Its score-
alignment change is divided by the block’s **numerical rank** before it
is used to calibrate one synthetic omitted index.

``` r
group_benchmark <- benchmark_from_model(
  fit,
  treatment = "treatment",
  benchmark = "factor(village)",
  vcov = "HC1"
)

group_result <- robust_rv(
  fit,
  treatment = "treatment",
  vcov = "HC1",
  benchmark = group_benchmark
)
```

## Point-conditioned diagnostics

The earlier signed point-omega object remains explicitly available:

``` r
point <- robust_rv(
  fit,
  treatment = "treatment",
  vcov = "HC1",
  calibration = "point",
  omega_target = 1.05,
  method = "protected_first"
)
```

Methods `"protected_first"`, `"diagonal_first"`, and `"full"` share the
same estimator-agnostic control flow. HC and CR differ only through the
supplied function $u\mapsto\omega_R(x,y,s,u)$.

## External validation

``` r
compare_sensemakr(fit, "treatment")
validate_sandwich_reference(fit, "treatment", vcov = "HC1")
```

The test suite additionally checks immutable numerical fixtures produced
by the independent Julia implementation. These fixtures are plain CSV
files; Julia is never invoked by the R package.

## Reproduce the Darfur example

With the suggested `sensemakr` and `sandwich` packages installed, run:

``` r
demo("darfur", package = "RobustSandwichRV")
```

The demo fits the published Darfur specification and computes bounded
HC1 and village-clustered CR1 results using `female` as the observed
benchmark.

A small single-threaded timing benchmark is installed with the package:

``` r
source(system.file(
  "benchmarks", "benchmark_runtime.R", package = "RobustSandwichRV"
))
```

## Interpretation

“Exact” refers to the realized sample statistic, the selected HC/CR
formula, the null value, and the fixed critical-value convention. It
does **not** imply exact finite-sample test size, confidence-interval
coverage, or absence of sampling uncertainty in a population sensitivity
parameter.
