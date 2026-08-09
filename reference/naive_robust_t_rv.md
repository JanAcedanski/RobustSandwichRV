# Incorrect robust-t plug-in benchmark

This deliberately inserts the robust short-model t statistic into the
conventional CH formula. It is provided only as a comparison benchmark.

## Usage

``` r
naive_robust_t_rv(
  cache,
  vcov = "HC1",
  cluster = NULL,
  alpha = 0.05,
  null = 0,
  critical = NULL
)
```

## Arguments

- cache:

  An `rv_cache`.

- vcov:

  Covariance estimator used for the baseline decision and critical value
  convention.

- cluster:

  Optional one-way cluster labels.

- alpha:

  Two-sided test level.

- null:

  Null coefficient value.

- critical:

  Optional fixed critical value.
