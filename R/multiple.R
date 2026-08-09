#' Reduce an omitted-confounder block to at most two active directions
#'
#' Preserves the treatment residual, full outcome residual, coefficient,
#' group partial R-squared values, and scalar score vector. The returned
#' `original_rank` must be retained for classical, HC1, and CR1 finite-sample
#' corrections.
#'
#' @param cache An `rv_cache`.
#' @param U Omitted-confounder matrix or orthonormal basis.
#' @return A list containing the reduced basis and original numerical rank.
#' @export
reduce_score_subspace <- function(cache, U) {
  U <- .rv_as_numeric_matrix(U, "U")
  if (nrow(U) != cache$n) stop("U has the wrong number of rows", call. = FALSE)
  Ux <- .rv_qr_basis(.rv_residualize(cache$Qx, U))
  k <- ncol(Ux)
  if (k == 0L) {
    return(list(
      basis = Ux, original_rank = 0L, q = cache$d,
      residuals = cache$e_short
    ))
  }
  treatment_fitted <- as.numeric(Ux %*% crossprod(Ux, cache$d))
  q <- cache$d - treatment_fitted
  qq <- .rv_dot(q, q)
  if (!(qq > 100 * .Machine$double.eps * .rv_dot(cache$d, cache$d))) {
    stop("treatment is unidentified after adding the subspace", call. = FALSE)
  }
  coefficients <- as.numeric(crossprod(cache$d, Ux)) /
    .rv_dot(cache$d, cache$d)
  outcome_space <- .rv_qr_basis(Ux - cache$d %o% coefficients)
  e_full <- as.numeric(.rv_residualize(outcome_space, cache$e_short))
  outcome_fitted <- cache$e_short - e_full
  lambda <- .rv_dot(q, outcome_fitted) / qq
  bridge <- outcome_fitted - lambda * cache$d
  candidates <- cbind(treatment_fitted, bridge)
  candidates <- as.matrix(Ux %*% crossprod(Ux, candidates))
  reduced <- .rv_qr_basis(candidates)
  if (ncol(reduced) > 2L) stop("internal two-direction reduction failure")
  list(
    basis = reduced,
    original_rank = k,
    q = q,
    residuals = e_full,
    treatment_fitted = treatment_fitted,
    outcome_fitted = outcome_fitted,
    bridge = bridge,
    lambda = lambda
  )
}

#' Verify the two-direction reduction numerically
#'
#' @inheritParams reduce_score_subspace
#' @param vcov Covariance estimator.
#' @param cluster Optional cluster labels.
#' @param tolerance Numerical comparison tolerance.
#' @export
verify_score_reduction <- function(cache, U, vcov = "HC1", cluster = NULL,
                                   tolerance = 1e-8) {
  Ux <- .rv_qr_basis(.rv_residualize(cache$Qx, U))
  original <- evaluate_subspace(
    cache, Ux, vcov = vcov, cluster = cluster,
    original_rank = ncol(Ux)
  )
  construction <- reduce_score_subspace(cache, Ux)
  reduced <- evaluate_subspace(
    cache, construction$basis, vcov = vcov, cluster = cluster,
    original_rank = construction$original_rank
  )
  errors <- c(
    q = max(abs(original$q - reduced$q)),
    residuals = max(abs(original$residuals - reduced$residuals)),
    tau = abs(original$tau - reduced$tau),
    variance = abs(original$variance_robust - reduced$variance_robust),
    r2_d = abs(original$r2_d - reduced$r2_d),
    r2_y = abs(original$r2_y - reduced$r2_y)
  )
  list(
    passed = all(errors <= tolerance),
    errors = errors,
    original = original,
    reduced = reduced,
    construction = construction
  )
}
