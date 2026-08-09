if (!requireNamespace("sensemakr", quietly = TRUE)) {
  stop("The Darfur demo requires the suggested package 'sensemakr'.")
}

data("darfur", package = "sensemakr")

darfur_model <- stats::lm(
  peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female + village,
  data = darfur
)

# The recommended observed-control workflow: female supplies the
# rank-normalized bounded calibration magnitude.
darfur_hc1 <- robust_rv(
  darfur_model,
  treatment = "directlyharmed",
  vcov = "HC1",
  benchmark = "female",
  seed = 20260809,
  n_starts = 8,
  max_iterations = 30
)
print(darfur_hc1)

# The control flow is unchanged for one-way clustering. Only the score-based
# omega evaluator and the CR1 finite-sample factor differ.
darfur_cr1 <- robust_rv(
  darfur_model,
  treatment = "directlyharmed",
  vcov = "CR1",
  cluster = ~ village,
  benchmark = "female",
  seed = 20260809,
  n_starts = 8,
  max_iterations = 30
)
print(darfur_cr1)

# Independent comparisons against sandwich and the conventional CH result.
validate_sandwich_reference(
  darfur_model, "directlyharmed", vcov = "HC1"
)
compare_sensemakr(darfur_model, "directlyharmed")
