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

  HC0, HC1, CR0, or CR1.

- cluster:

  Cluster labels or formula for CR estimators.

## Value

A data frame comparing standard errors and t statistics.
