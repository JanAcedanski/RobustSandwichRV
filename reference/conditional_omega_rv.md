# Point-conditioned attainable robustness value

Point-conditioned attainable robustness value

## Usage

``` r
conditional_omega_rv(
  cache,
  omega_target,
  method = c("protected_first", "diagonal_first", "full"),
  vcov = "HC1",
  cluster = NULL,
  alpha = 0.05,
  null = 0,
  critical = NULL,
  n_starts = 16L,
  seed = 1234L,
  max_iterations = 60L,
  root_tolerance = 1e-09,
  exact_tolerance = 2e-04,
  max_radius = 0.9
)

rv_att_point(
  cache,
  omega_target,
  method = c("protected_first", "diagonal_first", "full"),
  vcov = "HC1",
  cluster = NULL,
  alpha = 0.05,
  null = 0,
  critical = NULL,
  n_starts = 16L,
  seed = 1234L,
  max_iterations = 60L,
  root_tolerance = 1e-09,
  exact_tolerance = 2e-04,
  max_radius = 0.9
)
```

## Arguments

- cache:

  An `rv_cache`.

- omega_target:

  Target sandwich-to-classical standard-error ratio.

- method:

  `"protected_first"`, `"diagonal_first"`, or `"full"`.

- vcov:

  HC0, HC1, CR0, or CR1.

- cluster:

  Optional one-way cluster labels.

- alpha:

  Two-sided test level.

- null:

  Null coefficient value.

- critical:

  Optional fixed critical value.

- n_starts:

  Number of structured/random sphere starts.

- seed:

  Random seed used locally without changing `.Random.seed`.

- max_iterations:

  Maximum local sphere iterations.

- root_tolerance:

  Log-omega root tolerance.

- exact_tolerance:

  Gap used to label a numerically closed interval.

- max_radius:

  Largest partial-R-squared radius searched by the fallback.

## Value

A `robust_rv_result` containing a theorem certificate or certified
numerical bounds.
