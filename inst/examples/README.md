# Reproducible examples

Run the installed Darfur example with:

```r
demo("darfur", package = "RobustSandwichRV")
```

It uses the public `sensemakr::darfur` data and computes benchmark-bounded HC1
and village-clustered CR1 robustness values. The example is entirely R-based;
it does not call Julia or read any files produced by the Julia project.
