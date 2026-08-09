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

  One of `"classical"`, `"HC0"`, `"HC1"`, `"CR0"`, or `"CR1"`.

- cluster:

  Optional one-way cluster labels.

- null:

  Null coefficient value.
