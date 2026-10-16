# Independent R searches: no Julia witness vectors are supplied.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
pkgload::load_all("RobustSandwichRV-R", quiet = TRUE)
root <- "outputs/r_package_v05_validation_20261016"
inputs <- file.path(root, "article_inputs")
cases <- read.csv(file.path(inputs, "cases.csv"), check.names = FALSE)
out <- list()
selected <- unique(c(vapply(split(seq_len(nrow(cases)), cases$dataset), function(ix) ix[1], integer(1)),
  match(c("HAC_Bartlett_L4", "HAC_QS_BW5", "TWCR1_raw"), cases$vcov)))
for (ix in selected) {
  row <- cases[ix, ]; ds <- row$dataset
  data <- read.csv(file.path(inputs, paste0(ds, "_input.csv")), check.names = FALSE)
  columns <- read.csv(file.path(inputs, paste0(ds, "_columns.csv")), check.names = FALSE)
  fit <- lm(y ~ . - 1, data)
  cluster <- if (row$vcov == "CR1") scan(file.path(inputs, paste0(row$key, "_cluster.csv")), quiet = TRUE) else
    if (row$vcov == "TWCR1_raw") read.csv(file.path(inputs, paste0(row$key, "_cluster.csv")), header = FALSE) else NULL
  pool <- columns$name[columns$scalar %in% c(TRUE, "true", "TRUE") & columns$name != row$target]
  for (mode in c("increment", "profile", if (ds %in% c("growthdata", "mroz") ||
    row$vcov %in% c("HAC_Bartlett_L4", "HAC_QS_BW5", "TWCR1_raw")) "sand")) {
    common <- list(model = fit, treatment = row$target, vcov = row$vcov,
      cluster = cluster, critical = row$critical, n_starts = 4,
      max_iterations = 60, search_steps = 6, seed = 1234)
    value <- if (mode == "sand") do.call(rv_sand, common) else
      do.call(rv_cal, c(common, list(calibration = mode, benchmarks = pool)))
    stopifnot(is.finite(value$RV_upper), value$direct_refit_valid,
              value$RV_lower <= value$RV_upper)
    if (mode == "profile") stopifnot(abs(value$RV_lower - row$profile) < 5e-8)
    out[[length(out) + 1L]] <- data.frame(key = row$key, mode = mode,
      lower = value$RV_lower, upper = value$RV_upper, status = value$status,
      gap = value$RV_gap, seconds = value$runtime,
      search_error = value$search_error, refit_valid = value$direct_refit_valid)
    write.csv(do.call(rbind, out), file.path(root, "fresh_search_validation.csv"), row.names = FALSE)
    cat("FRESH", row$key, mode, value$status, value$RV_lower, value$RV_upper, "\n")
    flush.console()
  }
}
