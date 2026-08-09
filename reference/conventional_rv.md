# Conventional Cinelli-Hazlett robustness value

Returns the standard loss-of-significance robustness value. A short
model that does not reject has conventional RV zero under this
definition.

## Usage

``` r
conventional_rv(cache, alpha = 0.05, null = 0, critical = NULL)
```

## Arguments

- cache:

  An `rv_cache`.

- alpha:

  Two-sided test level.

- null:

  Null coefficient value.

- critical:

  Optional fixed critical value.
