# Calibrated robustness value

This is the cache-level dispatcher. Bounded calibration is the
recommended default; point conditioning remains available without
reinterpretation.

## Usage

``` r
calibrated_rv(
  cache,
  benchmark,
  calibration = c("bounded", "point"),
  point_omega = NULL,
  point_method = "protected_first",
  ...
)
```

## Arguments

- cache:

  An `rv_cache`.

- benchmark:

  Benchmark object or `b` in bounded mode; omega target in point mode
  when `point_omega` is not supplied.

- calibration:

  `"bounded"` or `"point"`.

- point_omega:

  Explicit point-conditioned omega target.

- point_method:

  Point solver method.

- ...:

  Arguments passed to
  [`rv_cal_bounded()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_cal_bounded.md)
  or
  [`conditional_omega_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/conditional_omega_rv.md).
