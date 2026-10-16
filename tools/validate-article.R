# Run from the workspace root after export-article-validation.jl.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
pkgload::load_all("RobustSandwichRV-R", quiet = TRUE)
root <- "outputs/r_package_v05_validation_20261016"
inputs <- file.path(root, "article_inputs")
cases <- read.csv(file.path(inputs, "cases.csv"), check.names = FALSE)
out <- list()
for (ds in unique(cases$dataset)) {
  data <- read.csv(file.path(inputs, paste0(ds, "_input.csv")), check.names = FALSE)
  columns <- read.csv(file.path(inputs, paste0(ds, "_columns.csv")), check.names = FALSE)
  fit <- lm(y ~ . - 1, data)
  subset <- cases[cases$dataset == ds, ]
  for (i in seq_len(nrow(subset))) {
    row <- subset[i, ]; key <- row$key
    cluster <- if (row$vcov == "CR1") scan(file.path(inputs, paste0(key, "_cluster.csv")), quiet = TRUE) else
      if (row$vcov == "TWCR1_raw") read.csv(file.path(inputs, paste0(key, "_cluster.csv")), header = FALSE) else NULL
    pool <- columns$name[columns$scalar %in% c(TRUE, "true", "TRUE") & columns$name != row$target]
    b <- rv_benchmarks(fit, row$target, pool, row$vcov, cluster)
    stopifnot(b$complete, abs(b$omega_full - row$omega) < 5e-8,
              abs(b$L - row$L) < 5e-8, abs(b$beta - row$beta) < 5e-8)
    witnesses <- as.matrix(read.csv(file.path(inputs, paste0(key, "_witness.csv")), header = FALSE))
    candidates <- lapply(seq_len(ncol(witnesses)), function(j) witnesses[, j])
    common <- list(model = fit, treatment = row$target, vcov = row$vcov, cluster = cluster,
      critical = row$critical, n_starts = 4, max_iterations = 40,
      run_search = FALSE, supplied_witnesses = candidates)
    delta <- do.call(rv_cal, c(common, list(calibration = "increment", benchmarks = pool)))
    profile <- do.call(rv_cal, c(common, list(calibration = "profile", benchmarks = pool)))
    sand <- do.call(rv_sand, common)
    stopifnot(abs(profile$RV_lower - row$profile) < 5e-8,
      abs(profile$RV_upper - row$profile) < 5e-8,
      delta$RV_lower <= row$delta_upper + 5e-8,
      delta$RV_upper <= row$delta_upper + 5e-8,
      sand$RV_lower <= row$sand_upper + 5e-8,
      sand$RV_upper <= row$sand_upper + 5e-8)
    out[[length(out) + 1L]] <- data.frame(dataset = ds, key = key, vcov = row$vcov,
      benchmark_max_error = max(abs(c(b$omega_full - row$omega, b$L - row$L, b$beta - row$beta))),
      profile_lower_error = profile$RV_lower - row$profile,
      profile_upper_error = profile$RV_upper - row$profile,
      profile_status = profile$status, delta_lower = delta$RV_lower,
      delta_upper = delta$RV_upper, delta_status = delta$status,
      sand_lower = sand$RV_lower, sand_upper = sand$RV_upper,
      all_refits_valid = all(vapply(list(sand, delta, profile), function(x) x$direct_refit_valid, logical(1))))
    write.csv(do.call(rbind, out), file.path(root, "article_validation.csv"), row.names = FALSE)
    cat("VALIDATED", key, "\n"); flush.console()
  }
}
writeLines(c(paste("PASS:", length(out), "current-article cases across all covariance estimators."),
  "Benchmarks recalculated by R from stored full model matrices and same-sample deletions.",
  "All profile lower bounds and point witnesses match the article within 5e-8.",
  "Every finite upper bound independently refitted by R; supplied Julia indices rechecked.",
  "This replay audits estimands, bounds and refits. It does not assert that local searches in different languages find identical minima."),
  file.path(root, "article_validation.txt"))
