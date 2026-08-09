# Prepare a linear-regression sensitivity cache

Computes the Frisch-Waugh-Lovell quantities used throughout the package.
`X` should contain the intercept and all observed controls.
Rank-deficient columns are absorbed by a rank-revealing QR
decomposition.

## Usage

``` r
prepare_model(y, d, X)
```

## Arguments

- y:

  Numeric outcome vector.

- d:

  Numeric treatment vector.

- X:

  Numeric control matrix, normally including an intercept.

## Value

An object of class `rv_cache`.
