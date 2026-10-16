test_that("fixed HAC agrees with independent sandwich covariance", {
  dat <- make_rv_fixture(n = 80, seed = 501)
  z <- sin(seq_len(dat$n))
  full <- lm(y ~ d + x1 + x2 + z, data = cbind(model.frame(dat$model), z = z))
  for (adjust in c(FALSE, TRUE)) for (kernel in c("bartlett", "QS")) {
    rule <- fixed_hac_rule(dat$n, kernel, if (kernel == "bartlett") 4 else 5,
                            adjust = adjust)
    external <- if (kernel == "bartlett")
      sandwich::NeweyWest(full, lag = 4, prewhite = FALSE, adjust = adjust) else
      sandwich::kernHAC(full, kernel = "Quadratic Spectral", bw = 5,
                         prewhite = FALSE, adjust = adjust)
    fast <- evaluate_confounder(dat$cache, z, vcov = rule)
    refit <- evaluate_confounder_reference(dat$cache, z, vcov = rule)
    expect_equal(fast$variance_robust, external["d", "d"], tolerance = 2e-10)
    expect_equal(refit$variance_robust, external["d", "d"], tolerance = 2e-10)
    expect_equal(fast$t_robust, refit$t_robust, tolerance = 2e-10)
    comparison <- validate_sandwich_reference(dat$model, "d", vcov = rule)
    expect_equal(diff(comparison$standard_error), 0, tolerance = 1e-10)
  }
  expect_equal(evaluate_confounder(dat$cache, z, vcov = "HAC_Bartlett_L4")$variance_robust,
    evaluate_confounder(dat$cache, z, vcov = fixed_hac_rule(80))$variance_robust)
  expect_equal(evaluate_confounder(dat$cache, z, vcov = "HAC_QS_BW5")$variance_robust,
    evaluate_confounder(dat$cache, z, vcov = fixed_hac_rule(80, "QS", 5))$variance_robust)
})

test_that("raw two-way covariance uses separate component corrections", {
  dat <- make_rv_fixture(n = 96, seed = 502)
  c1 <- rep(1:12, each = 8); c2 <- rep(1:8, 12)
  z <- cos(seq_len(dat$n))
  full <- lm(y ~ d + x1 + x2 + z, data = cbind(model.frame(dat$model), z = z))
  for (variant in c("CR0", "CR1")) {
    rule <- fixed_twoway_cluster_rule(c1, c2, variant)
    V <- sandwich::vcovCL(full, cluster = data.frame(c1, c2),
      type = if (variant == "CR1") "HC1" else "HC0", cadjust = variant == "CR1",
      multi0 = FALSE, fix = FALSE)
    expect_gt(V["d", "d"], 0)
    fast <- evaluate_confounder(dat$cache, z, vcov = rule)
    refit <- evaluate_confounder_reference(dat$cache, z, vcov = rule)
    expect_equal(fast$variance_robust, V["d", "d"], tolerance = 1e-10)
    expect_equal(refit$variance_robust, V["d", "d"], tolerance = 1e-10)
    alias <- paste0("TW", variant, "_raw")
    expect_equal(evaluate_confounder(dat$cache, z, vcov = alias,
      cluster = data.frame(c1, c2))$variance_robust, fast$variance_robust)
    ctx <- RobustSandwichRV:::.rv_loss_context(dat$model, "d", rule, NULL, .05, NULL, 0)
    expect_equal(ctx$critical, qt(.975, min(length(unique(c1)), length(unique(c2))) - 1))
    expect_false(rule$psd)
    expect_false(rule$metadata$psd_repair)
    comparison <- validate_sandwich_reference(dat$model, "d", vcov = rule)
    expect_equal(diff(comparison$standard_error), 0, tolerance = 1e-10)
  }
})

test_that("HAC ordering and two-way formulas preserve retained sample alignment", {
  dat <- make_rv_fixture(n = 96, seed = 502)
  set.seed(29); permutation <- sample(96)
  rule <- fixed_hac_rule(96, "QS", 5, order = permutation)
  scores <- rnorm(96)
  default <- fixed_hac_rule(96, "QS", 5)
  expect_equal(sum(scores * rule$apply(scores)),
               sum(scores[permutation] * default$apply(scores[permutation])), tolerance = 1e-12)
  c1 <- rep(1:12, each = 8); c2 <- rep(1:8, 12)
  fit <- dat$model
  a <- rv_sand(fit, "d", vcov = "TWCR1_raw", cluster = ~ c1 + c2, run_search = FALSE)
  b <- rv_sand(fit, "d", vcov = "TWCR1_raw", cluster = data.frame(c1, c2), run_search = FALSE)
  expect_equal(a$RV_lower, b$RV_lower)
  expect_error(fixed_hac_rule(48, bandwidth = 48), "smaller")
  expect_error(fixed_hac_rule(48, order = rep(1, 48)), "permutation")
  expect_error(fixed_twoway_cluster_rule(rep(1, 48), c2), "two clusters")
  expect_error(rv_sand(dat$model, "d", vcov = fixed_hac_rule(50)), "observation count")
  expect_error(rv_sand(dat$model, "d", vcov = "TWCR1_raw", cluster = c1), "two-column")
})

test_that("signed quadratic gradients, roots and envelopes handle indefinite operators", {
  dat <- make_rv_fixture(n = 48, seed = 504)
  rules <- list(fixed_hac_rule(48), fixed_hac_rule(48, "QS", 5),
    fixed_twoway_cluster_rule(rep(1:8, each = 6), rep(1:6, 8)))
  for (rule in rules) {
    ctx <- RobustSandwichRV:::.rv_loss_context(dat$model, "d", rule, NULL, .05, 2, 0)
    set.seed(73); v <- rnorm(48); h <- 1e-6
    state <- function(v) RobustSandwichRV:::.rv_loss_direction_objective(ctx, .2, .15, -1, v)
    a <- state(v)
    numerical <- vapply(seq_along(v), function(i) {
      delta <- numeric(48); delta[i] <- h
      (state(v + delta)$value - state(v - delta)$value) / (2 * h)
    }, numeric(1))
    expect_equal(a$gradient, numerical, tolerance = 1e-6)
    env <- RobustSandwichRV:::.rv_loss_envelope(ctx)
    for (r in c(0, .001, .1, .6)) for (s in c(-1, 1)) {
      value <- RobustSandwichRV:::.rv_loss_direction_objective(ctx, r, r / 2, s, v)$value
      expect_lte(value, env$bound(r)^2 + 1e-10)
    }
  }
  # Deterministic indefinite example with a negative-to-positive path.
  ctx <- RobustSandwichRV:::.rv_loss_context(dat$model, "d", rules[[3]], NULL, .05, 2, 0)
  ctl <- RobustSandwichRV:::.rv_loss_controls(8, 80, 29, 5e-8, .99, 6, TRUE)
  interval <- RobustSandwichRV:::.rv_loss_direction_interval(ctx, .2, .2, 1, ctl)
  expect_lt(interval$quadratic_low, 0)
  expect_gt(interval$quadratic_high, 0)
  root <- RobustSandwichRV:::.rv_loss_omega_root(ctx, .2, .2, 1, interval,
    sqrt(interval$quadratic_high / 2), 29)
  expect_true(root$found)
  expect_equal(RobustSandwichRV:::.rv_loss_direction_objective(ctx, .2, .2, 1, root$u)$value,
               interval$quadratic_high / 2, tolerance = 1e-9)
  z <- single_parameterization(dat$cache, sqrt(.2), sqrt(.2), interval$u_low)
  expect_error(evaluate_confounder_reference(dat$cache, z, vcov = rules[[3]]), "nonpositive")
})

test_that("all new estimators work with both current calibrated RVs", {
  dat <- make_rv_fixture(n = 72, seed = 505, target_t = 4)
  specifications <- list(fixed_hac_rule(72), fixed_hac_rule(72, "QS", 5),
    fixed_twoway_cluster_rule(rep(1:12, each = 6), rep(1:6, 12)))
  for (rule in specifications) {
    result <- sandwich_rv(dat$model, "d", c("x1", "x2"), vcov = rule,
      critical = 2, n_starts = 4, max_iterations = 40, search_steps = 4)
    expect_equal(nrow(result$table), 4)
    for (value in result$results) {
      expect_true(is.finite(value$RV_upper))
      expect_true(value$direct_refit_valid)
      expect_lte(value$RV_lower, value$RV_upper)
      refit <- evaluate_confounder_reference(dat$cache, value$best_witness_z, vcov = rule)
      expect_lte(abs(refit$t_robust), 2 + 5e-8)
      expect_gt(refit$variance_robust, 0)
    }
    b <- result$benchmark
    expect_true(b$complete)
    if (rule$metadata$family == "HAC") {
      expect_equal(b$records$omega_minus_common_chi,
        b$records$omega_minus * sqrt(dat$cache$df_short / (dat$cache$df_short + 1)))
    }
  }
})
