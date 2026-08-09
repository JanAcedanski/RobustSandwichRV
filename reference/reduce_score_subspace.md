# Reduce an omitted-confounder block to at most two active directions

Preserves the treatment residual, full outcome residual, coefficient,
group partial R-squared values, and scalar score vector. The returned
`original_rank` must be retained for classical, HC1, and CR1
finite-sample corrections.

## Usage

``` r
reduce_score_subspace(cache, U)
```

## Arguments

- cache:

  An `rv_cache`.

- U:

  Omitted-confounder matrix or orthonormal basis.

## Value

A list containing the reduced basis and original numerical rank.
