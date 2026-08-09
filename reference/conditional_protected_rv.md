# Analytical protected conditional robustness value

Minimizes `max(x, y)` over the complete conventional CH decision
surface, before imposing sandwich attainability.

## Usage

``` r
conditional_protected_rv(
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

A list containing the lower bound, protected coordinates, branch, and
baseline decision.
