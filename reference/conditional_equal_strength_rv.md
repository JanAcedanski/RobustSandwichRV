# Equal-strength conditional robustness value

Equal-strength conditional robustness value

## Usage

``` r
conditional_equal_strength_rv(
  cache,
  omega,
  alpha = 0.05,
  null = 0,
  critical = NULL,
  vcov = "HC1",
  cluster = NULL
)
```

## Arguments

- cache:

  An `rv_cache`.

- omega:

  Positive target sandwich-to-classical standard-error ratio.

- alpha:

  Two-sided test level.

- null:

  Null coefficient value.

- critical:

  Optional fixed critical value.

- vcov:

  Covariance estimator used for the baseline decision and critical value
  convention.

- cluster:

  Optional one-way cluster labels.

## Value

A scalar robustness value.
