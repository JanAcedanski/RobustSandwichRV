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

  HC0/HC1, one-way CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw,
  TWCR1_raw, or a fixed covariance rule.

- cluster:

  One-way labels or a two-column data frame/matrix for two-way
  clustering.

## Value

A function computing `se_R / se_classical`.
