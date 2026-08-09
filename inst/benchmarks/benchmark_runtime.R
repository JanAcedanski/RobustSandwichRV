# Reproducible single-threaded timing smoke benchmark.
# Run from an installed package with:
# source(system.file("benchmarks", "benchmark_runtime.R",
#                    package = "RobustSandwichRV"))

set.seed(20260809)
n <- 1000L
X <- cbind(1, matrix(stats::rnorm(n * 8L), ncol = 8L))
d <- X[, 2L] * 0.3 + stats::rnorm(n)
y <- 0.4 * d + X[, 3L] * 0.2 +
  (0.5 + 0.3 * abs(d)) * stats::rnorm(n)
z <- stats::rnorm(n)
cache <- prepare_model(y, d, X)

fast <- system.time(for (i in seq_len(100L)) {
  evaluate_confounder(cache, z, vcov = "HC1")
})
reference <- system.time(for (i in seq_len(100L)) {
  evaluate_confounder_reference(cache, z, vcov = "HC1")
})

timings <- data.frame(
  implementation = c("fast FWL", "explicit full regression"),
  evaluations = 100L,
  elapsed_seconds = c(fast[["elapsed"]], reference[["elapsed"]]),
  row.names = NULL
)
print(timings)

bounded <- system.time({
  result <- rv_cal_bounded(
    cache, benchmark = 0.05, vcov = "HC1",
    n_starts = 8, max_iterations = 30
  )
})
print(result)
print(data.frame(workflow = "bounded HC1",
                 elapsed_seconds = bounded[["elapsed"]]))
