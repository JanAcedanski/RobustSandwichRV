# Create a sensitivity cache from an `lm` model

Create a sensitivity cache from an `lm` model

## Usage

``` r
cache_from_lm(model, treatment)
```

## Arguments

- model:

  An unweighted [`stats::lm()`](https://rdrr.io/r/stats/lm.html) object.

- treatment:

  Name of the single treatment model-matrix column.

## Value

An `rv_cache` with the original model metadata attached.
