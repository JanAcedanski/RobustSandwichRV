# A complete example using data included with R.
library(RobustSandwichRV)

fit <- lm(mpg ~ wt + hp + disp, data = mtcars)

# wt: target coefficient; hp and disp: observed benchmark controls.
rv <- sandwich_rv(
  fit, treatment = "wt", vcov = "HC1",
  benchmarks = c("hp", "disp"), k = 1, alpha = 0.05
)
print(rv)
print(as.data.frame(rv))

# estimate is a point value only if the numerical bounds meet.
# Otherwise report lower and upper; these are computational bounds, not CIs.
# Values are fractions: 0.10 means 10% on the partial-R-squared scale.
