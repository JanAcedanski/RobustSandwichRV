# Construct a scalar confounder from partial correlations

Construct a scalar confounder from partial correlations

## Usage

``` r
single_parameterization(cache, rho_d, rho_y, u)
```

## Arguments

- cache:

  An `rv_cache`.

- rho_d:

  Signed treatment-side partial correlation.

- rho_y:

  Signed outcome-side partial correlation.

- u:

  Residual direction. It is projected onto the complement of
  `span(X, d, e_short)` and normalized.

## Value

A unit-norm residualized confounder vector.
