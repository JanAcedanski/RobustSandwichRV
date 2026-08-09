# Benchmark-bounded calibrated robustness value

Computes the analytical adverse-endpoint lower bound and searches for an
explicit, independently refitted decision-changing confounder whose
sandwich ratio lies anywhere inside the benchmark band.

## Usage

``` r
rv_cal_bounded(
  cache,
  benchmark,
  vcov = "HC1",
  cluster = NULL,
  bound_direction = c("symmetric", "observed"),
  alpha = 0.05,
  null = 0,
  critical = NULL,
  n_starts = 16L,
  seed = 1234L,
  max_iterations = 60L,
  root_tolerance = 1e-09,
  exact_tolerance = 2e-04,
  max_radius = 0.9,
  supplied_witnesses = list(),
  run_witness_search = TRUE
)
```

## Arguments

- cache:

  An `rv_cache`.

- benchmark:

  A nonnegative `b` or an
  [`benchmark_calibration()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/benchmark_calibration.md)
  object.

- vcov:

  HC0, HC1, CR0, or CR1.

- cluster:

  Optional one-way cluster labels.

- bound_direction:

  Symmetric or observed-direction band.

- alpha:

  Two-sided test level.

- null:

  Null coefficient value.

- critical:

  Optional fixed critical value.

- n_starts:

  Number of sphere starts.

- seed:

  Local random seed.

- max_iterations:

  Maximum local sphere iterations.

- root_tolerance:

  Log-omega root tolerance.

- exact_tolerance:

  Numerical interval-closing tolerance.

- max_radius:

  Largest searched partial-R-squared radius.

- supplied_witnesses:

  Optional earlier witness vectors for nested bands.

- run_witness_search:

  If `FALSE`, stop after the endpoint certificate.

## Value

A `robust_rv_result` with certified lower and upper endpoints.
