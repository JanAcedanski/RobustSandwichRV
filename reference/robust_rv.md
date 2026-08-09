# Main model-facing robustness-value workflow

Main model-facing robustness-value workflow

## Usage

``` r
robust_rv(
  model,
  treatment = NULL,
  vcov = "HC1",
  calibration = c("bounded", "point"),
  benchmark = NULL,
  b = NULL,
  omega_target = NULL,
  cluster = NULL,
  ...
)

observed_control_rv(
  model,
  treatment = NULL,
  vcov = "HC1",
  calibration = c("bounded", "point"),
  benchmark = NULL,
  b = NULL,
  omega_target = NULL,
  cluster = NULL,
  ...
)
```

## Arguments

- model:

  An unweighted `lm` object or an `rv_cache`.

- treatment:

  Name of the scalar treatment model-matrix column. Not used when
  `model` is already an `rv_cache`.

- vcov:

  HC0, HC1, CR0, or CR1.

- calibration:

  `"bounded"` (recommended) or `"point"`.

- benchmark:

  In bounded mode, an observed control term, vector of terms, benchmark
  object, or direct nonnegative `b`. In point mode, an omega target or
  object containing `omega_target`.

- b:

  Optional direct bounded benchmark magnitude.

- omega_target:

  Optional point-conditioned ratio.

- cluster:

  Optional cluster formula or labels.

- ...:

  Numerical solver options.

## Value

A `robust_rv_result`.
