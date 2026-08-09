.rv_H <- function(a) {
  (sqrt(a^4 + 4 * a^2) - a^2) / 2
}

.rv_short_conventional <- function(cache, null = 0) {
  dd <- .rv_dot(cache$d, cache$d)
  se <- sqrt((.rv_dot(cache$e_short, cache$e_short) / cache$df_short) / dd)
  if (!(se > 0)) {
    stop("the short model has zero conventional standard error", call. = FALSE)
  }
  list(se = se, t = (cache$tau_short - null) / se)
}

.rv_short_quantities <- function(cache, null = 0) {
  short <- .rv_short_conventional(cache, null)
  list(
    nu = cache$df_short,
    t_short_classical = short$t,
    f = abs(short$t) / sqrt(cache$df_short)
  )
}

#' Construct a scalar confounder from partial correlations
#'
#' @param cache An `rv_cache`.
#' @param rho_d Signed treatment-side partial correlation.
#' @param rho_y Signed outcome-side partial correlation.
#' @param u Residual direction. It is projected onto the complement of
#'   `span(X, d, e_short)` and normalized.
#' @return A unit-norm residualized confounder vector.
#' @export
single_parameterization <- function(cache, rho_d, rho_y, u) {
  if (!inherits(cache, "rv_cache")) {
    stop("cache must be an rv_cache", call. = FALSE)
  }
  if (!.rv_is_scalar_number(rho_d) || abs(rho_d) >= 1) {
    stop("abs(rho_d) must be strictly smaller than one", call. = FALSE)
  }
  if (!.rv_is_scalar_number(rho_y) || abs(rho_y) > 1) {
    stop("abs(rho_y) must not exceed one", call. = FALSE)
  }
  u <- .rv_normalize_direction(cache, u)
  d_hat <- cache$d / .rv_norm(cache$d)
  e_hat <- cache$e_short / .rv_norm(cache$e_short)
  z <- rho_d * d_hat + sqrt(1 - rho_d^2) *
    (rho_y * e_hat + sqrt(max(0, 1 - rho_y^2)) * u)
  z / .rv_norm(z)
}

#' Exact normalized score directions
#'
#' @param x Treatment-side partial R-squared.
#' @param y Outcome-side partial R-squared.
#' @param bias_sign Either `-1` or `1`.
#' @inheritParams single_parameterization
#' @return A list containing normalized treatment and outcome residuals.
#' @export
alignment_directions <- function(cache, x, y, bias_sign, u) {
  .rv_assert_probability(x, "x", upper_open = TRUE)
  .rv_assert_probability(y, "y", upper_open = TRUE)
  if (!(bias_sign %in% c(-1, 1))) {
    stop("bias_sign must be -1 or 1", call. = FALSE)
  }
  u <- .rv_normalize_direction(cache, u)
  d_hat <- cache$d / .rv_norm(cache$d)
  e_hat <- cache$e_short / .rv_norm(cache$e_short)
  q_hat <- sqrt(1 - x) * d_hat -
    bias_sign * sqrt(x * y) * e_hat - sqrt(x * (1 - y)) * u
  r_hat <- sqrt(1 - y) * e_hat - bias_sign * sqrt(y) * u
  list(q_hat = q_hat, r_hat = r_hat, u = u)
}

#' Construct the estimator-specific attainable omega evaluator
#'
#' The returned function has signature `function(x, y, bias_sign, u)` and is
#' the only estimator-specific input needed by the common point and bounded
#' calibration algorithms.
#'
#' @param cache An `rv_cache`.
#' @param vcov HC0, HC1, CR0, or CR1.
#' @param cluster Optional one-way cluster labels.
#' @return A function computing `se_R / se_classical`.
#' @export
conditional_omega_evaluator <- function(cache, vcov = "HC1", cluster = NULL) {
  vcov <- match.arg(vcov, c("HC0", "HC1", "CR0", "CR1"))
  cluster <- .rv_resolve_cluster(cluster, cache)
  df_full <- cache$df_short - 1L
  if (vcov %in% c("HC0", "HC1")) {
    scale2 <- if (vcov == "HC0") df_full else cache$n
    return(function(x, y, bias_sign, u) {
      directions <- alignment_directions(cache, x, y, bias_sign, u)
      score <- directions$q_hat * directions$r_hat
      sqrt(scale2 * sum(score^2))
    })
  }
  info <- .rv_cluster_codes(cluster, cache$n)
  scale2 <- if (vcov == "CR0") {
    df_full
  } else {
    info$G / (info$G - 1) * (cache$n - 1)
  }
  function(x, y, bias_sign, u) {
    directions <- alignment_directions(cache, x, y, bias_sign, u)
    score <- directions$q_hat * directions$r_hat
    totals <- rowsum(score, info$codes, reorder = FALSE)
    sqrt(scale2 * sum(totals^2))
  }
}

.rv_classical_ch_t_signed <- function(cache, x, y, bias_sign, null = 0) {
  .rv_assert_probability(x, "x", upper_open = TRUE)
  .rv_assert_probability(y, "y", upper_open = TRUE)
  if (!(bias_sign %in% c(-1, 1))) {
    stop("bias_sign must be -1 or 1", call. = FALSE)
  }
  short <- .rv_short_quantities(cache, null)
  first_term <- sign(short$t_short_classical) * short$f *
    sqrt((1 - x) / (1 - y))
  bias_term <- bias_sign * sqrt(x * y / (1 - y))
  sqrt(short$nu - 1) * (first_term - bias_term)
}

#' Conventional full-model t statistic from CH geometry
#'
#' @inheritParams alignment_directions
#' @param null Null coefficient value.
#' @param adverse If `TRUE`, select the coefficient-reducing sign.
#' @export
classical_ch_t <- function(cache, x, y, null = 0, adverse = TRUE) {
  short <- .rv_short_quantities(cache, null)
  coefficient_sign <- if (short$t_short_classical < 0) -1 else 1
  bias_sign <- if (adverse) coefficient_sign else -coefficient_sign
  .rv_classical_ch_t_signed(cache, x, y, bias_sign, null)
}

.rv_protected_regime <- function(branch) {
  if (branch %in% c("equal_strength_decrease", "equal_strength_increase")) {
    "equal_strength"
  } else if (branch %in% c("asymmetric_decrease", "asymmetric_increase")) {
    "interior"
  } else {
    "zero_strength"
  }
}

.rv_conditional_equal <- function(cache, omega, critical, reject_short,
                                  null = 0) {
  short <- .rv_short_quantities(cache, null)
  kappa <- critical * omega / sqrt(short$nu - 1)
  a <- if (reject_short) max(short$f - kappa, 0) else max(kappa - short$f, 0)
  .rv_H(a)
}

.rv_conditional_protected <- function(cache, omega, critical, reject_short,
                                      null = 0) {
  short <- .rv_short_quantities(cache, null)
  kappa <- critical * omega / sqrt(short$nu - 1)
  f <- short$f
  fk <- f * kappa
  transition_tol <- 64 * .Machine$double.eps * max(1, abs(fk))
  asymmetric <- fk > 1 + transition_tol
  if (reject_short) {
    if (kappa >= f) {
      return(list(rv = 0, x = 0, y = 0, branch = "already_nonrejecting"))
    }
    if (asymmetric) {
      r <- (f^2 - kappa^2) / (f^2 + 1)
      return(list(rv = r, x = r, y = r / (f^2 * (1 - r)),
                  branch = "asymmetric_decrease"))
    }
    r <- .rv_H(f - kappa)
    return(list(rv = r, x = r, y = r,
                branch = "equal_strength_decrease"))
  }
  if (kappa <= f) {
    return(list(rv = 0, x = 0, y = 0, branch = "already_rejecting"))
  }
  if (asymmetric) {
    r <- (kappa^2 - f^2) / (kappa^2 + 1)
    return(list(rv = r, x = r / (f^2 + r), y = r,
                branch = "asymmetric_increase"))
  }
  r <- .rv_H(kappa - f)
  list(rv = r, x = r, y = r, branch = "equal_strength_increase")
}

#' Equal-strength conditional robustness value
#'
#' @param cache An `rv_cache`.
#' @param omega Positive target sandwich-to-classical standard-error ratio.
#' @param alpha Two-sided test level.
#' @param null Null coefficient value.
#' @param critical Optional fixed critical value.
#' @param vcov Covariance estimator used for the baseline decision and critical
#'   value convention.
#' @param cluster Optional one-way cluster labels.
#' @return A scalar robustness value.
#' @export
conditional_equal_strength_rv <- function(cache, omega, alpha = 0.05,
                                          null = 0, critical = NULL,
                                          vcov = "HC1", cluster = NULL) {
  if (!.rv_is_scalar_number(omega) || omega <= 0) {
    stop("omega must be positive and finite", call. = FALSE)
  }
  vcov <- .rv_match_vcov(vcov)
  cluster <- .rv_resolve_cluster(cluster, cache)
  critical <- .rv_critical_value(cache, vcov, 1L, alpha, cluster, critical)
  short <- .rv_short_evaluation(cache, vcov, cluster, null)
  .rv_conditional_equal(cache, omega, critical,
                        .rv_reject(short$t_robust, critical), null)
}

#' Analytical protected conditional robustness value
#'
#' Minimizes `max(x, y)` over the complete conventional CH decision surface,
#' before imposing sandwich attainability.
#'
#' @inheritParams conditional_equal_strength_rv
#' @return A list containing the lower bound, protected coordinates, branch,
#'   and baseline decision.
#' @export
conditional_protected_rv <- function(cache, omega, alpha = 0.05,
                                     null = 0, critical = NULL,
                                     vcov = "HC1", cluster = NULL) {
  if (!.rv_is_scalar_number(omega) || omega <= 0) {
    stop("omega must be positive and finite", call. = FALSE)
  }
  vcov <- .rv_match_vcov(vcov)
  cluster <- .rv_resolve_cluster(cluster, cache)
  critical <- .rv_critical_value(cache, vcov, 1L, alpha, cluster, critical)
  short <- .rv_short_evaluation(cache, vcov, cluster, null)
  reject_short <- .rv_reject(short$t_robust, critical)
  out <- .rv_conditional_protected(cache, omega, critical, reject_short, null)
  out$critical_value <- critical
  out$reject_short <- reject_short
  out$protected_regime <- .rv_protected_regime(out$branch)
  out
}

#' Point-estimate robustness value
#'
#' @param cache An `rv_cache`.
#' @param null Target coefficient value.
#' @export
point_estimate_rv <- function(cache, null = 0) {
  short <- .rv_short_quantities(cache, null)
  .rv_H(short$f)
}

#' Conventional Cinelli-Hazlett robustness value
#'
#' Returns the standard loss-of-significance robustness value. A short model
#' that does not reject has conventional RV zero under this definition.
#'
#' @param cache An `rv_cache`.
#' @param alpha Two-sided test level.
#' @param null Null coefficient value.
#' @param critical Optional fixed critical value.
#' @export
conventional_rv <- function(cache, alpha = 0.05, null = 0, critical = NULL) {
  critical <- .rv_critical_value(cache, "classical", 1L, alpha, NULL, critical)
  short <- .rv_short_evaluation(cache, "classical", null = null)
  if (!.rv_reject(short$t_classical, critical)) {
    return(0)
  }
  .rv_conditional_protected(cache, 1, critical, TRUE, null)$rv
}

#' Incorrect robust-t plug-in benchmark
#'
#' This deliberately inserts the robust short-model t statistic into the
#' conventional CH formula. It is provided only as a comparison benchmark.
#'
#' @inheritParams conditional_equal_strength_rv
#' @export
naive_robust_t_rv <- function(cache, vcov = "HC1", cluster = NULL,
                              alpha = 0.05, null = 0, critical = NULL) {
  vcov <- .rv_match_vcov(vcov)
  cluster <- .rv_resolve_cluster(cluster, cache)
  critical <- .rv_critical_value(cache, vcov, 1L, alpha, cluster, critical)
  short <- .rv_short_evaluation(cache, vcov, cluster, null)
  if (!.rv_reject(short$t_robust, critical)) {
    return(0)
  }
  nu <- cache$df_short
  f <- abs(short$t_robust) / sqrt(nu)
  kappa <- critical / sqrt(nu - 1)
  if (f * kappa > 1 + 64 * .Machine$double.eps) {
    (f^2 - kappa^2) / (f^2 + 1)
  } else {
    .rv_H(max(f - kappa, 0))
  }
}

#' Compare the classical implementation with `sensemakr`
#'
#' @param model An `lm` object.
#' @param treatment Treatment coefficient name.
#' @param alpha Two-sided test level.
#' @return A data frame with both conventional robustness values.
#' @export
compare_sensemakr <- function(model, treatment, alpha = 0.05) {
  if (!requireNamespace("sensemakr", quietly = TRUE)) {
    stop("package 'sensemakr' is required for this comparison", call. = FALSE)
  }
  cache <- cache_from_lm(model, treatment)
  external <- sensemakr::robustness_value(
    model = model, covariates = treatment, alpha = alpha
  )
  data.frame(
    method = c("RobustSandwichRV", "sensemakr"),
    RV = c(conventional_rv(cache, alpha = alpha), as.numeric(external)[1L]),
    row.names = NULL
  )
}
