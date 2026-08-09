# Construct the estimator-specific attainable omega evaluator

The returned function has signature `function(x, y, bias_sign, u)` and
is the only estimator-specific input needed by the common point and
bounded calibration algorithms.

## Usage

``` r
conditional_omega_evaluator(cache, vcov = "HC1", cluster = NULL)
```

## Arguments

- cache:

  An `rv_cache`.

- vcov:

  HC0, HC1, CR0, or CR1.

- cluster:

  Optional one-way cluster labels.

## Value

A function computing `se_R / se_classical`.
