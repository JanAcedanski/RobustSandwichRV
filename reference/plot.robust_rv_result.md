# Plot the analytical decision geometry

Shows the fixed-omega CH decision region at the adverse endpoint
(bounded mode) or target omega (point mode), together with the protected
solution and the best verified witness.

## Usage

``` r
# S3 method for class 'robust_rv_result'
plot(x, cache = NULL, maximum = NULL, grid_size = 100L, ...)
```

## Arguments

- x:

  A `robust_rv_result`.

- cache:

  The `rv_cache` used to compute `x`.

- maximum:

  Largest displayed partial R-squared.

- grid_size:

  Grid resolution.

- ...:

  Additional arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).
