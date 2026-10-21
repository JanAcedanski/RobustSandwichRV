# Fast FWL evaluation of an omitted confounder

The function accepts either `evaluate_confounder(cache, z, ...)` or
`evaluate_confounder(y, d, X, z, ...)`.

## Usage

``` r
evaluate_confounder(
  y,
  d = NULL,
  X = NULL,
  z = NULL,
  vcov = "HC1",
  cluster = NULL,
  null = 0
)
```

## Arguments

- y:

  An `rv_cache` or outcome vector.

- d:

  Treatment vector, or `z` when `y` is an `rv_cache`.

- X:

  Control matrix.

- z:

  Confounder vector or matrix.

- vcov:

  Classical, a built-in sandwich estimator, or a fixed covariance rule
  (see
  [`fixed_hac_rule`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/fixed_score_rules.md)).

- cluster:

  One-way labels or two-column data frame/matrix for two-way clustering.

- null:

  Null coefficient value.
