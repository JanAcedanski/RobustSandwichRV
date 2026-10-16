test_that("restoring benchmarks match independent sandwich refits", {
  dat <- make_rv_fixture(n = 96, seed = 180)
  for (vc in c("HC0", "HC1", "CR0", "CR1")) {
    cluster <- if (grepl("CR", vc)) dat$cluster else NULL
    b <- rv_benchmarks(dat$model, "d", c("x1", "x2"), vc, cluster)
    expect_true(b$complete)
    omega <- function(fit) {
      V <- if (grepl("HC", vc)) sandwich::vcovHC(fit, type = vc) else
        sandwich::vcovCL(fit, cluster = cluster,
                         type = if (vc == "CR1") "HC1" else "HC0",
                         cadjust = vc == "CR1")
      sqrt(V["d", "d"] / vcov(fit)["d", "d"])
    }
    full <- omega(dat$model)
    for (i in seq_len(nrow(b$records))) {
      term <- b$records$benchmark[i]
      reduced <- update(dat$model, paste(". ~ . -", term))
      expected <- omega(reduced)
      expect_equal(b$records$omega_minus[i], expected, tolerance = 1e-10)
      expect_equal(b$records$delta_omega[i], full - expected, tolerance = 1e-10)
      factor <- if (vc %in% c("HC0", "CR0"))
        sqrt(df.residual(dat$model) / df.residual(reduced)) else 1
      expect_equal(b$records$log_increment[i], log(full / (expected * factor)),
                   tolerance = 1e-10)
      dfit <- lm(reformulate(setdiff(c("x1", "x2"), term), "d"), model.frame(dat$model))
      restored_dfit <- lm(d ~ x1 + x2, model.frame(dat$model))
      expect_equal(b$records$r2_d[i], 1 - deviance(restored_dfit) / deviance(dfit),
                   tolerance = 1e-10)
      expect_equal(b$records$r2_y[i], 1 - deviance(dat$model) / deviance(reduced),
                   tolerance = 1e-10)
    }
    expect_equal(b$M_delta, full + max(0, b$records$delta_omega))
    expect_equal(b$L, max(0, b$records$log_increment))
  }
})

test_that("the new optimizer gradient and universal envelope are valid", {
  dat <- make_rv_fixture(n = 48, seed = 182)
  for (vc in c("HC1", "CR1")) {
    cluster <- if (vc == "CR1") dat$cluster else NULL
    ctx <- RobustSandwichRV:::.rv_loss_context(dat$model, "d", vc, cluster, .05, 2, 0)
    set.seed(81); v <- rnorm(48); h <- 1e-6
    state <- function(v) RobustSandwichRV:::.rv_loss_direction_objective(ctx, .13, .21, -1, v)
    a <- state(v); numeric <- vapply(seq_along(v), function(i) {
      offset <- numeric(length(v)); offset[i] <- h
      (state(v + offset)$value - state(v - offset)$value) / (2 * h)
    }, numeric(1))
    expect_equal(a$gradient, numeric, tolerance = 1e-7)
    env <- RobustSandwichRV:::.rv_loss_envelope(ctx)
    evaluate <- conditional_omega_evaluator(dat$cache, vc, cluster)
    for (r in c(0, .0001, .01, .1, .7)) for (s in c(-1, 1)) {
      u <- RobustSandwichRV:::.rv_structured_directions(dat$cache, 4, 9)[[1]]
      expect_lte(evaluate(r, r, s, u), env$bound(r) + 1e-10)
      expect_lte(evaluate(r / 2, r, s, u), env$bound(r) + 1e-10)
    }
  }
})

test_that("loss values preserve bounds and independently refitted witnesses", {
  dat <- make_rv_fixture(n = 64, seed = 183, target_t = 4)
  for (vc in c("HC0", "HC1", "CR0", "CR1")) {
    cluster <- if (grepl("CR", vc)) dat$cluster else NULL
    values <- list(
      rv_sand(dat$model, "d", vc, cluster, critical = 2,
              n_starts = 4, max_iterations = 30, search_steps = 4),
      rv_cal(dat$model, "d", "increment", c("x1", "x2"), vc, cluster,
             critical = 2, n_starts = 4, max_iterations = 30, search_steps = 4),
      rv_cal(dat$model, "d", "profile", c("x1", "x2"), vc, cluster,
             critical = 2, n_starts = 4, max_iterations = 30, search_steps = 4))
    for (value in values) {
      expect_gte(value$RV_lower, 0)
      expect_true(is.finite(value$RV_upper))
      expect_lte(value$RV_lower, value$RV_upper)
      expect_true(value$direct_refit_valid)
      a <- evaluate_confounder_reference(dat$cache, value$best_witness_z,
                                          vcov = vc, cluster = cluster)
      expect_equal(value$RV_upper, max(a$r2_d, a$r2_y), tolerance = 1e-10)
      expect_lte(abs(a$t_robust), 2 + 5e-8)
      cap <- if (value$calibration_mode == "profile") {
        p <- value$profile
        p$a0 * exp(p$k * min(p$L, p$beta * (sqrt(a$r2_d) + sqrt(a$r2_y))))
      } else value$omega_upper
      expect_lte(a$omega, cap * (1 + 5e-8))
      expect_identical(is.finite(value$rv), value$exact)
      expect_s3_class(as.data.frame(value), "data.frame")
    }
    protected <- conditional_protected_rv(dat$cache, values[[2]]$omega_upper,
                                           vcov = vc, cluster = cluster, critical = 2)
    expect_gte(values[[2]]$RV_lower + 1e-10, protected$rv)
  }
})

test_that("empty, zero, invalid and block benchmark pools stay distinct", {
  dat <- make_rv_fixture(n = 64, seed = 185, target_t = 4)
  empty <- rv_cal(dat$model, "d", benchmarks = character(), run_search = FALSE)
  expect_identical(empty$status, "CALIBRATION_UNAVAILABLE")
  expect_identical(empty$substatus, "EMPTY_POOL")
  expect_true(is.na(empty$RV_lower))
  expect_error(rv_benchmarks(dat$model, "d", "d"), "treatment")
  expect_error(rv_benchmarks(dat$model, "d", "(Intercept)"), "intercept")
  block <- rv_cal(dat$model, "d", "profile", list(block = c("x1", "x2")),
                   run_search = FALSE)
  expect_identical(block$status, "CALIBRATION_UNAVAILABLE")
  expect_identical(block$benchmark$records$rank, 2L)
  expect_true(is.finite(block$benchmark$M_delta))
  expect_error(rv_benchmarks(dat$model, "d", list(a = "x1", b = c("x1", "x2"))), "overlap")
  # An independent reduced fit can have zero score variance while the full
  # model's variance is positive. The prescribed pool must remain incomplete.
  set.seed(1); d <- rep(c(1, -1, 0, 0), 8); y <- rep(c(0, 0, 1, -1), 8)
  x <- rnorm(32); fit <- lm(y ~ d + x)
  bad <- rv_benchmarks(fit, "d", "x")
  expect_false(bad$complete)
  expect_identical(bad$status, "INCOMPLETE_POOL")
  expect_identical(rv_cal(fit, "d", benchmarks = "x")$status, "CALIBRATION_UNAVAILABLE")
  expect_error(rv_benchmarks(lm((d + x) ~ d + x), "d"), "zero residual")
  # Search deterministic examples until all signed restoring changes are negative.
  found <- FALSE
  for (seed in 1:30) {
    sample <- make_rv_fixture(n = 48, seed = seed)
    b <- rv_benchmarks(sample$model, "d")
    if (b$delta_positive == 0 && b$L == 0) { found <- TRUE; break }
  }
  expect_true(found)
  expect_identical(b$status, "COMPLETE_ZERO_INCREMENT")
  a <- rv_cal(sample$model, "d", "profile", k = 1, run_search = FALSE)
  b2 <- rv_cal(sample$model, "d", "profile", k = 2, run_search = FALSE)
  expect_equal(a$RV_lower, b2$RV_lower)
  expect_equal(a$omega_upper, b2$omega_upper)
  expect_identical(a$benchmark_status, "zero")
})

test_that("search failure, common cutoffs and caller RNG are reported correctly", {
  dat <- make_rv_fixture(n = 64, seed = 183, target_t = 4)
  set.seed(310); before <- .Random.seed
  sand <- rv_sand(dat$model, "d", run_search = FALSE)
  expect_identical(.Random.seed, before)
  expect_identical(sand$status, "LOWER_BOUND_ONLY")
  expect_true(is.na(sand$RV_upper))
  expect_gt(sand$RV_lower, 0)
  result <- sandwich_rv(dat$model, "d", critical = 2, n_starts = 4,
                         max_iterations = 25, search_steps = 3)
  expect_identical(.Random.seed, before)
  expect_equal(nrow(as.data.frame(result)), 4)
  expect_true(all(result$table$critical == 2))
  expect_lte(result$results$RV_sand$RV_upper,
              min(result$results$RV_cal_delta$RV_upper,
                  result$results$RV_cal_profile$RV_upper) + 5e-8)
  nonreject <- make_rv_fixture(n = 64, seed = 183, target_t = .5)
  for (value in list(rv_sand(nonreject$model, "d"),
                     rv_cal(nonreject$model, "d", "profile"))) {
    expect_identical(value$status, "ALREADY_NONREJECTING")
    expect_equal(value$rv, 0)
    expect_equal(value$RV_upper, 0)
  }
  expect_error(rv_sand(dat$model, "d", vcov = "HC3"), "arg")
  expect_error(rv_cal(dat$model, "d", "profile", k = 0), "k must")
  expect_error(rv_cal(dat$model, "d", nonsense = TRUE), "unknown")
})
