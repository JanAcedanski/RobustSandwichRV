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
  [`prepare_model()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/prepare_model.md)
  cache.

- U:

  Matrix whose columns span the omitted-confounder subspace.

- vcov:

  One of `"classical"`, `"HC0"`, `"HC1"`, `"CR0"`, or `"CR1"`.

- cluster:

  Optional one-way cluster labels.

- original_rank:

  Original omitted-block rank used in degrees-of-freedom corrections.
  This may exceed the geometric rank of a reduced two-direction
  representation.

- null:

  Null coefficient value.

## Value

An `rv_confounder_evaluation`.
