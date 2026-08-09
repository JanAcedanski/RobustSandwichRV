# Exact normalized score directions

Exact normalized score directions

## Usage

``` r
alignment_directions(cache, x, y, bias_sign, u)
```

## Arguments

- cache:

  An `rv_cache`.

- x:

  Treatment-side partial R-squared.

- y:

  Outcome-side partial R-squared.

- bias_sign:

  Either `-1` or `1`.

- u:

  Residual direction. It is projected onto the complement of
  `span(X, d, e_short)` and normalized.

## Value

A list containing normalized treatment and outcome residuals.
