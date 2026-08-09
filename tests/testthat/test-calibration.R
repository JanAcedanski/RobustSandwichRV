test_that("scalar and grouped benchmarks use numerical-rank normalization", {
  set.seed(30)
  n <- 120
  d <- rnorm(n)
  x <- rnorm(n)
  group <- factor(rep(letters[1:6], length.out = n))
  y <- 0.5 * d + 0.2 * x + model.matrix(~ group)[, -1] %*%
    seq(0.1, 0.5, length.out = 5) + rnorm(n)
  model <- lm(as.numeric(y) ~ d + x + group)
  scalar <- benchmark_from_model(model, "d", "x", vcov = "HC1")
  grouped <- benchmark_from_model(model, "d", "group", vcov = "HC1")
  expect_equal(scalar$group_rank, 1L)
  expect_equal(grouped$group_rank, 5L)
  expect_equal(grouped$delta_L_signed,
               0.5 * log(grouped$L_full / grouped$L_reduced) / 5,
               tolerance = 1e-13)
})

test_that("HC0/HC1 and CR0/CR1 alignment calibrations agree", {
  fixture <- make_rv_fixture()
  hc0 <- benchmark_from_model(fixture$model, "d", "x1", vcov = "HC0")
  hc1 <- benchmark_from_model(fixture$model, "d", "x1", vcov = "HC1")
  cr0 <- benchmark_from_model(
    fixture$model, "d", "x1", vcov = "CR0", cluster = fixture$cluster
  )
  cr1 <- benchmark_from_model(
    fixture$model, "d", "x1", vcov = "CR1", cluster = fixture$cluster
  )
  expect_equal(hc0$delta_L_signed, hc1$delta_L_signed, tolerance = 1e-13)
  expect_equal(cr0$delta_L_signed, cr1$delta_L_signed, tolerance = 1e-13)
})

test_that("symmetric and directional bands are constructed exactly", {
  symmetric <- bounded_benchmark_band(1.2, 0.1, "symmetric")
  expect_equal(symmetric$omega_lower, 1.2 * exp(-0.1))
  expect_equal(symmetric$omega_upper, 1.2 * exp(0.1))
  positive <- bounded_benchmark_band(
    1.2, 0.1, "observed", delta_L_signed = 0.1
  )
  negative <- bounded_benchmark_band(
    1.2, 0.1, "observed", delta_L_signed = -0.1
  )
  expect_equal(positive$omega_lower, 1.2)
  expect_equal(negative$omega_upper, 1.2)
  expect_error(bounded_benchmark_band(
    1.2, 0.1, "observed", delta_L_signed = -0.2
  ))
})

test_that("point and bounded solvers return independently refitted witnesses", {
  fixture <- make_rv_fixture(n = 80, seed = 31)
  cache <- fixture$cache
  origin <- RobustSandwichRV:::.rv_omega_origin(cache, "HC1")
  point <- conditional_omega_rv(
    cache, origin, vcov = "HC1", n_starts = 8,
    max_iterations = 20, max_radius = 0.5
  )
  expect_true(is.finite(point$RV_lower))
  expect_true(point$direct_refit_valid)
  expect_equal(point$best_witness_omega, origin, tolerance = 1e-6)

  bounded <- rv_cal_bounded(
    cache, 0, vcov = "HC1", n_starts = 8,
    max_iterations = 20, max_radius = 0.5
  )
  expect_equal(bounded$omega_lower, origin, tolerance = 1e-12)
  expect_equal(bounded$omega_upper, origin, tolerance = 1e-12)
  expect_equal(bounded$RV_lower, point$RV_lower, tolerance = 1e-10)
  expect_true(is.finite(bounded$RV_upper))
  expect_true(bounded$direct_refit_valid)
  expect_true(bounded$witness_literal_decision_change)
  if (identical(bounded$certificate_type, "adverse_endpoint_equality")) {
    expect_true(bounded$theorem_exact)
  }
  expect_true(bounded$numerically_closed)
})

test_that("CR1 bounded calibration and nested witness propagation work", {
  fixture <- make_rv_fixture(n = 84, seed = 33)
  benchmarks <- list(0.01, 0.04)
  results <- rv_cal_bounded_sequence(
    fixture$cache, benchmarks, vcov = "CR1", cluster = fixture$cluster,
    n_starts = 8, max_iterations = 20, max_radius = 0.6
  )
  expect_length(results, 2)
  expect_lte(results[[2]]$RV_lower, results[[1]]$RV_lower + 1e-12)
  expect_true(is.finite(results[[1]]$RV_upper))
  expect_lte(results[[2]]$RV_upper, results[[1]]$RV_upper + 1e-8)
})

test_that("failed search is a lower-bound status, never infeasible", {
  fixture <- make_rv_fixture(n = 80, seed = 35)
  result <- rv_cal_bounded(
    fixture$cache, 5, vcov = "HC1", n_starts = 8,
    max_iterations = 5, max_radius = 1e-7,
    run_witness_search = FALSE
  )
  expect_equal(result$status, "LOWER_BOUND_ONLY")
  expect_false(grepl("INFEASIBLE", result$status, fixed = TRUE))
})
