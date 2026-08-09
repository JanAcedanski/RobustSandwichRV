test_that("two-direction reduction preserves score-based inference", {
  fixture <- make_rv_fixture(n = 110, seed = 40)
  set.seed(41)
  Z <- matrix(rnorm(fixture$n * 6), ncol = 6)
  for (vcov in c("HC0", "HC1", "CR0", "CR1")) {
    result <- verify_score_reduction(
      fixture$cache, Z, vcov = vcov,
      cluster = if (grepl("CR", vcov)) fixture$cluster else NULL,
      tolerance = 1e-8
    )
    expect_true(result$passed)
    expect_lte(ncol(result$construction$basis), 2L)
    expect_equal(result$reduced$original_rank,
                 result$original$original_rank)
  }
})

test_that("lm-facing API performs bounded and point calibration", {
  fixture <- make_rv_fixture(n = 80, seed = 42)
  bounded <- robust_rv(
    fixture$model, treatment = "d", vcov = "HC1",
    calibration = "bounded", benchmark = "x1",
    n_starts = 8, max_iterations = 15, max_radius = 0.5
  )
  expect_s3_class(bounded, "robust_rv_result")
  expect_equal(bounded$benchmark_name, "x1")
  expect_true(bounded$direct_refit_valid)
  expect_s3_class(summary(bounded), "summary.robust_rv_result")
  expect_s3_class(as.data.frame(bounded), "data.frame")

  point <- robust_rv(
    fixture$model, treatment = "d", vcov = "HC1",
    calibration = "point", omega_target = 1,
    n_starts = 8, max_iterations = 15, max_radius = 0.5
  )
  expect_equal(point$calibration_mode, "point")
})

test_that("input failures are informative", {
  fixture <- make_rv_fixture(n = 60, seed = 44)
  expect_error(cache_from_lm(fixture$model, "does_not_exist"),
               "exactly match")
  offset_model <- lm(y ~ d + x1 + x2 + offset(rep(0, fixture$n)),
                     data = data.frame(
                       y = fixture$y, d = fixture$d,
                       x1 = fixture$x1, x2 = fixture$x2
                     ))
  expect_error(cache_from_lm(offset_model, "d"), "offset")
  expect_error(evaluate_confounder(
    fixture$cache, fixture$cache$d, vcov = "HC1"
  ), "outcome-side|unidentified")
  expect_error(rv_cal_bounded(
    fixture$cache, -0.1, vcov = "HC1"
  ), "benchmark")
  expect_error(rv_cal_bounded(
    fixture$cache, 0.1, vcov = "CR1"
  ), "cluster")
})
