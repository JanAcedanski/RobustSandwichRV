# Reproducible examples

Run the quick start without extra data packages:

```r
demo("quick-start", package = "RobustSandwichRV", ask = FALSE)
```

It fits an `lm` model to the built-in `mtcars` data and reports `RV_sand`,
`RV_cal_delta`, and `RV_cal_profile`, alongside the common-cutoff CH reference.
Replace the regression, target coefficient, and benchmark names with your own.
Point estimates are reported only when the numerical bounds meet; otherwise use
the returned lower and upper bounds. These are computational bounds, not
confidence intervals.

Run the installed Darfur example with:

```r
demo("darfur", package = "RobustSandwichRV")
```

It uses the public `sensemakr::darfur` data and computes benchmark-bounded HC1
and village-clustered CR1 robustness values. The example is entirely R-based;
it does not call Julia or read any files produced by the Julia project.
