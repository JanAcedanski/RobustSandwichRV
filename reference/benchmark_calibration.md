# Construct a rank-normalized benchmark calibration

Construct a rank-normalized benchmark calibration

## Usage

``` r
benchmark_calibration(
  delta_L = NULL,
  benchmark_type = c("scalar", "group"),
  benchmark_name = "benchmark",
  group_rank = 1L,
  rank_normalized = TRUE,
  L_full = NULL,
  L_reduced = NULL
)
```

## Arguments

- delta_L:

  Signed per-synthetic-index change in log square-root score alignment.

- benchmark_type:

  `"scalar"` or `"group"`.

- benchmark_name:

  Human-readable benchmark name.

- group_rank:

  Numerical rank of the observed benchmark block.

- rank_normalized:

  Whether `delta_L` is already divided by group rank.

- L_full, L_reduced:

  Optional score-alignment values. When supplied, `delta_L` is computed
  as `log(L_full / L_reduced) / (2 * group_rank)`.

## Value

An `rv_benchmark` object.
