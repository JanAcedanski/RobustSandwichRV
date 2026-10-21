# Evaluate a residualized confounder subspace

Evaluate a residualized confounder subspace

## Usage

``` r
evaluate_subspace(
  cache,
  U,
  vcov = "HC1",
  cluster = NULL,
  original_rank = ncol(as.matrix(U)),
  null = 0
)
```

## Arguments

- cache:

  An
  [`prepare_model()`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/prepare_model.md)
  cache.

- U:

  Matrix whose columns span the omitted-confounder subspace.

- vcov:

  Classical, a built-in sandwich estimator, or a fixed covariance rule
  (see
  [`fixed_hac_rule`](https://JanAcedanski.github.io/RVwithSandwichSE/reference/fixed_score_rules.md)).

- cluster:

  One-way labels or two-column data frame/matrix for two-way clustering.

- original_rank:

  Original omitted-block rank used in degrees-of-freedom corrections.
  This may exceed the geometric rank of a reduced two-direction
  representation.

- null:

  Null coefficient value.

## Value

An `rv_confounder_evaluation`.
