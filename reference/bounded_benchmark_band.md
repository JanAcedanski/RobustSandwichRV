# Construct a benchmark sandwich-ratio band

Construct a benchmark sandwich-ratio band

## Usage

``` r
bounded_benchmark_band(
  omega_origin,
  b,
  bound_direction = c("symmetric", "observed"),
  delta_L_signed = NA_real_
)
```

## Arguments

- omega_origin:

  Positive zero-strength synthetic-index ratio.

- b:

  Nonnegative log-scale benchmark magnitude.

- bound_direction:

  `"symmetric"` or `"observed"`.

- delta_L_signed:

  Signed observed change, required for the observed direction mode.
