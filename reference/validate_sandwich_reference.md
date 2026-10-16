# Validate a short-model covariance against `sandwich`

Validate a short-model covariance against `sandwich`

## Usage

``` r
validate_sandwich_reference(model, treatment, vcov = "HC1", cluster = NULL)
```

## Arguments

- model:

  An `lm` model.

- treatment:

  Treatment coefficient name.

- vcov:

  HC0/HC1, one-way CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw,
  TWCR1_raw, or a fixed covariance rule.

- cluster:

  One-way labels, a two-column data frame/matrix, or a formula for one
  or two dimensions.

## Value

A data frame comparing standard errors and t statistics.
