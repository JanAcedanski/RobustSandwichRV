test_that("single-confounder parameterization realizes partial R-squared", {
  fixture <- make_rv_fixture()
  cache <- fixture$cache
  set.seed(10)
  u <- rnorm(cache$n)
  x <- 0.12
  y <- 0.08
  for (sign in c(-1, 1)) {
    z <- single_parameterization(cache, sqrt(x), sign * sqrt(y), u)
    evaluation <- evaluate_confounder(cache, z, vcov = "HC1")
    expect_equal(evaluation$r2_d, x, tolerance = 1e-10)
    expect_equal(evaluation$r2_y, y, tolerance = 1e-10)
    expect_equal(evaluation$t_classical,
                 RobustSandwichRV:::.rv_classical_ch_t_signed(
                   cache, x, y, sign
                 ), tolerance = 1e-9)
  }
})

test_that("omega adapters reproduce direct HC and CR evaluations", {
  fixture <- make_rv_fixture()
  cache <- fixture$cache
  set.seed(11)
  u <- rnorm(cache$n)
  x <- 0.07
  y <- 0.11
  sign <- 1
  z <- single_parameterization(cache, sqrt(x), sign * sqrt(y), u)
  for (vcov in c("HC0", "HC1", "CR0", "CR1")) {
    evaluator <- conditional_omega_evaluator(
      cache, vcov, cluster = if (grepl("CR", vcov)) fixture$cluster else NULL
    )
    explicit <- evaluate_confounder(
      cache, z, vcov = vcov,
      cluster = if (grepl("CR", vcov)) fixture$cluster else NULL
    )
    normalized_u <- RobustSandwichRV:::.rv_normalize_direction(cache, u)
    expect_equal(evaluator(x, y, sign, normalized_u), explicit$omega,
                 tolerance = 1e-9)
  }
})

test_that("protected formulas classify transition and reverse edge correctly", {
  fixture <- make_rv_fixture()
  cache <- fixture$cache
  base_noise <- cache$e_short
  edge_f <- 2
  edge_y <- base_noise + edge_f * norm(base_noise, "2") /
    norm(cache$d, "2") * cache$d
  edge_cache <- prepare_model(edge_y, cache$d, cache$X)
  edge_short <- RobustSandwichRV:::.rv_short_quantities(edge_cache)
  expect_equal(edge_short$f, edge_f, tolerance = 1e-10)

  kappa <- 2.5
  omega <- kappa * sqrt(edge_short$nu - 1)
  reverse <- RobustSandwichRV:::.rv_conditional_protected(
    edge_cache, omega, 1, FALSE
  )
  expect_equal(reverse$branch, "asymmetric_increase")
  expect_lt(reverse$x, reverse$y)

  equality_f <- 0.5
  equality_y <- base_noise + equality_f * norm(base_noise, "2") /
    norm(cache$d, "2") * cache$d
  equality_cache <- prepare_model(equality_y, cache$d, cache$X)
  equality_short <- RobustSandwichRV:::.rv_short_quantities(equality_cache)
  equality_omega <- (1 / equality_f) * sqrt(equality_short$nu - 1)
  equality <- RobustSandwichRV:::.rv_conditional_protected(
    equality_cache, equality_omega, 1, FALSE
  )
  expect_equal(equality$branch, "equal_strength_increase")
  expect_equal(equality$x, equality$y, tolerance = 1e-13)
})

test_that("conventional RV agrees with sensemakr", {
  skip_if_not_installed("sensemakr")
  fixture <- make_rv_fixture()
  comparison <- compare_sensemakr(fixture$model, "d")
  expect_equal(comparison$RV[1], comparison$RV[2], tolerance = 1e-10)
})

