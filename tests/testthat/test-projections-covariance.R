test_that("QR residualization and FWL identities are exact", {
  fixture <- make_rv_fixture()
  cache <- fixture$cache
  expect_lt(max(abs(crossprod(cache$Qx, cache$d))), 1e-12)
  expect_lt(max(abs(crossprod(cache$Qx, cache$e_short))), 1e-12)
  expect_lt(abs(sum(cache$d * cache$e_short)), 1e-12)

  set.seed(2)
  z <- rnorm(cache$n)
  for (vcov in c("classical", "HC0", "HC1")) {
    fast <- evaluate_confounder(cache, z, vcov = vcov)
    reference <- evaluate_confounder_reference(cache, z, vcov = vcov)
    expect_equal(fast$tau, reference$tau, tolerance = 1e-11)
    expect_equal(fast$residuals, reference$residuals, tolerance = 1e-10)
    expect_equal(fast$variance_robust, reference$variance_robust,
                 tolerance = 1e-10)
    expect_equal(fast$t_robust, reference$t_robust, tolerance = 1e-10)
    expect_equal(fast$omega, fast$se_robust / fast$se_classical,
                 tolerance = 1e-13)
    expect_equal(fast$t_robust, fast$t_classical / fast$omega,
                 tolerance = 1e-12)
  }
})

test_that("CR0 and CR1 agree with explicit cluster bread-meat-bread", {
  fixture <- make_rv_fixture()
  set.seed(3)
  z <- rnorm(fixture$n)
  for (vcov in c("CR0", "CR1")) {
    fast <- evaluate_confounder(
      fixture$cache, z, vcov = vcov, cluster = fixture$cluster
    )
    reference <- evaluate_confounder_reference(
      fixture$cache, z, vcov = vcov, cluster = fixture$cluster
    )
    expect_equal(fast$tau, reference$tau, tolerance = 1e-11)
    expect_equal(fast$variance_robust, reference$variance_robust,
                 tolerance = 1e-10)
    expect_equal(fast$omega, reference$omega, tolerance = 1e-10)
  }
})

test_that("short-model covariance matches sandwich", {
  skip_if_not_installed("sandwich")
  fixture <- make_rv_fixture()
  for (vcov in c("HC0", "HC1")) {
    comparison <- validate_sandwich_reference(
      fixture$model, "d", vcov = vcov
    )
    expect_equal(comparison$standard_error[1],
                 comparison$standard_error[2], tolerance = 1e-10)
  }
  for (vcov in c("CR0", "CR1")) {
    comparison <- validate_sandwich_reference(
      fixture$model, "d", vcov = vcov, cluster = fixture$cluster
    )
    expect_equal(comparison$standard_error[1],
                 comparison$standard_error[2], tolerance = 1e-10)
  }
})

test_that("scale and sign changes of a scalar confounder are invariant", {
  fixture <- make_rv_fixture()
  set.seed(4)
  z <- rnorm(fixture$n)
  first <- evaluate_confounder(fixture$cache, z, vcov = "HC1")
  second <- evaluate_confounder(fixture$cache, -7 * z, vcov = "HC1")
  expect_equal(first$tau, second$tau, tolerance = 1e-12)
  expect_equal(first$variance_robust, second$variance_robust,
               tolerance = 1e-12)
  expect_equal(first$r2_d, second$r2_d, tolerance = 1e-12)
  expect_equal(first$r2_y, second$r2_y, tolerance = 1e-12)
})

