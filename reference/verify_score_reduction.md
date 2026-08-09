# Verify the two-direction reduction numerically

Verify the two-direction reduction numerically

## Usage

``` r
verify_score_reduction(
  cache,
  U,
  vcov = "HC1",
  cluster = NULL,
  tolerance = 1e-08
)
```

## Arguments

- cache:

  An `rv_cache`.

- U:

  Omitted-confounder matrix or orthonormal basis.

- vcov:

  Covariance estimator.

- cluster:

  Optional cluster labels.

- tolerance:

  Numerical comparison tolerance.
