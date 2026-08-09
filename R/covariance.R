.rv_cluster_scalar_meat <- function(q, e, cluster) {
  info <- .rv_cluster_codes(cluster, length(q))
  totals <- rowsum(q * e, info$codes, reorder = FALSE)
  list(meat = sum(totals^2), G = info$G)
}

.rv_scalar_variance <- function(vcov, q, e, df_full, cluster = NULL) {
  vcov <- .rv_match_vcov(vcov)
  n <- length(q)
  qq <- .rv_dot(q, q)
  rss <- .rv_dot(e, e)
  if (!(qq > 0) || df_full <= 0L) {
    stop("the full treatment coefficient is not identified", call. = FALSE)
  }
  classical <- (rss / df_full) / qq
  if (vcov == "classical") {
    return(list(classical = classical, robust = classical, correction = 1,
                meat = NA_real_, n_groups = n))
  }
  if (vcov %in% c("HC0", "HC1")) {
    meat <- sum((q * e)^2)
    correction <- if (vcov == "HC1") n / df_full else 1
    return(list(
      classical = classical,
      robust = correction * meat / qq^2,
      correction = correction,
      meat = meat,
      n_groups = n
    ))
  }
  cluster_result <- .rv_cluster_scalar_meat(q, e, cluster)
  correction <- if (vcov == "CR1") {
    cluster_result$G / (cluster_result$G - 1) * (n - 1) / df_full
  } else {
    1
  }
  list(
    classical = classical,
    robust = correction * cluster_result$meat / qq^2,
    correction = correction,
    meat = cluster_result$meat,
    n_groups = cluster_result$G
  )
}

.rv_assemble_evaluation <- function(cache, U, q, e_full, tau,
                                    original_rank, vcov, cluster, null) {
  df_full <- cache$n - cache$rank_x - 1L - as.integer(original_rank)
  variances <- .rv_scalar_variance(vcov, q, e_full, df_full, cluster)
  if (variances$classical < 0 || variances$robust < 0) {
    stop("a negative variance was produced", call. = FALSE)
  }
  se_classical <- sqrt(variances$classical)
  se_robust <- sqrt(variances$robust)
  if (!(se_classical > 0) || !(se_robust > 0)) {
    stop("the requested full-model standard error is zero", call. = FALSE)
  }
  dd <- .rv_dot(cache$d, cache$d)
  ee <- .rv_dot(cache$e_short, cache$e_short)
  r2_d <- min(max(1 - .rv_dot(q, q) / dd, 0), 1)
  r2_y <- if (ee > 100 * .Machine$double.eps) {
    min(max(1 - .rv_dot(e_full, e_full) / ee, 0), 1)
  } else {
    0
  }
  structure(
    list(
      tau = as.numeric(tau),
      null = as.numeric(null),
      se_classical = se_classical,
      se_robust = se_robust,
      variance_classical = variances$classical,
      variance_robust = variances$robust,
      omega = se_robust / se_classical,
      t_classical = (tau - null) / se_classical,
      t_robust = (tau - null) / se_robust,
      q = as.numeric(q),
      residuals = as.numeric(e_full),
      confounder_basis = U,
      r2_d = r2_d,
      r2_y = r2_y,
      geometric_rank = ncol(U),
      original_rank = as.integer(original_rank),
      df_residual = df_full,
      vcov = vcov,
      finite_sample_correction = variances$correction,
      meat = variances$meat,
      n_groups = variances$n_groups
    ),
    class = "rv_confounder_evaluation"
  )
}

#' Evaluate a residualized confounder subspace
#'
#' @param cache An [prepare_model()] cache.
#' @param U Matrix whose columns span the omitted-confounder subspace.
#' @param vcov One of `"classical"`, `"HC0"`, `"HC1"`, `"CR0"`, or
#'   `"CR1"`.
#' @param cluster Optional one-way cluster labels.
#' @param original_rank Original omitted-block rank used in degrees-of-freedom
#'   corrections. This may exceed the geometric rank of a reduced
#'   two-direction representation.
#' @param null Null coefficient value.
#' @return An `rv_confounder_evaluation`.
#' @export
evaluate_subspace <- function(cache, U, vcov = "HC1", cluster = NULL,
                              original_rank = ncol(as.matrix(U)), null = 0) {
  if (!inherits(cache, "rv_cache")) {
    stop("cache must be an rv_cache", call. = FALSE)
  }
  vcov <- .rv_match_vcov(vcov)
  cluster <- .rv_resolve_cluster(cluster, cache)
  U <- .rv_as_numeric_matrix(U, "U")
  if (nrow(U) != cache$n) {
    stop("U has the wrong number of rows", call. = FALSE)
  }
  Ux <- .rv_qr_basis(.rv_residualize(cache$Qx, U))
  geometric_rank <- ncol(Ux)
  if (original_rank < geometric_rank) {
    stop("original_rank cannot be smaller than the represented subspace rank",
         call. = FALSE)
  }
  if (geometric_rank == 0L) {
    projection_d <- numeric(cache$n)
  } else {
    projection_d <- as.numeric(Ux %*% crossprod(Ux, cache$d))
  }
  q <- cache$d - projection_d
  qq <- .rv_dot(q, q)
  if (!(qq > 100 * .Machine$double.eps * .rv_dot(cache$d, cache$d))) {
    stop("the confounder subspace makes treatment unidentified",
         call. = FALSE)
  }
  if (geometric_rank == 0L) {
    outcome_space <- matrix(numeric(), nrow = cache$n, ncol = 0L)
  } else {
    coefficients <- as.numeric(crossprod(cache$d, Ux)) /
      .rv_dot(cache$d, cache$d)
    MdU <- Ux - cache$d %o% coefficients
    outcome_space <- .rv_qr_basis(MdU)
  }
  e_full <- as.numeric(.rv_residualize(outcome_space, cache$e_short))
  tau <- .rv_dot(q, cache$y_x) / qq
  .rv_assemble_evaluation(cache, Ux, q, e_full, tau, original_rank,
                          vcov, cluster, null)
}

.rv_signed_partial_correlations <- function(cache, z) {
  zz <- .rv_dot(z, z)
  dd <- .rv_dot(cache$d, cache$d)
  a <- .rv_dot(cache$d, z) / dd
  w <- z - a * cache$d
  ww <- .rv_dot(w, w)
  ee <- .rv_dot(cache$e_short, cache$e_short)
  rho_d <- .rv_dot(z, cache$d) / sqrt(zz * dd)
  rho_y <- if (ww > 0 && ee > 0) {
    .rv_dot(w, cache$e_short) / sqrt(ww * ee)
  } else {
    0
  }
  c(rho_d = rho_d, rho_y = rho_y)
}

.rv_evaluate_single_fast <- function(cache, z, vcov, cluster, null) {
  z <- as.numeric(z)
  if (length(z) != cache$n || any(!is.finite(z))) {
    stop("z must be a finite vector with one entry per observation",
         call. = FALSE)
  }
  z_x <- as.numeric(.rv_residualize(cache$Qx, z))
  zz <- .rv_dot(z_x, z_x)
  if (!(zz > 100 * .Machine$double.eps * max(.rv_dot(z, z), 1))) {
    stop("z lies in span(X) or is numerically zero", call. = FALSE)
  }
  a <- .rv_dot(z_x, cache$d)
  q <- cache$d - z_x * a / zz
  qq <- .rv_dot(q, q)
  dd <- .rv_dot(cache$d, cache$d)
  if (!(qq > 100 * .Machine$double.eps * dd)) {
    stop("z makes treatment unidentified", call. = FALSE)
  }
  w <- z_x - cache$d * a / dd
  ww <- .rv_dot(w, w)
  if (!(ww > 100 * .Machine$double.eps * zz)) {
    stop("z has no outcome-side direction after conditioning on treatment",
         call. = FALSE)
  }
  e_full <- cache$e_short - w * .rv_dot(w, cache$e_short) / ww
  tau <- .rv_dot(q, cache$y_x) / qq
  U <- .rv_qr_basis(z_x)
  result <- .rv_assemble_evaluation(cache, U, q, e_full, tau, 1L, vcov,
                                    cluster, null)
  correlations <- .rv_signed_partial_correlations(cache, z_x)
  result$rho_d <- correlations[["rho_d"]]
  result$rho_y <- correlations[["rho_y"]]
  result
}

#' Fast FWL evaluation of an omitted confounder
#'
#' The function accepts either `evaluate_confounder(cache, z, ...)` or
#' `evaluate_confounder(y, d, X, z, ...)`.
#'
#' @param y An `rv_cache` or outcome vector.
#' @param d Treatment vector, or `z` when `y` is an `rv_cache`.
#' @param X Control matrix.
#' @param z Confounder vector or matrix.
#' @inheritParams evaluate_subspace
#' @export
evaluate_confounder <- function(y, d = NULL, X = NULL, z = NULL,
                                vcov = "HC1", cluster = NULL, null = 0) {
  if (inherits(y, "rv_cache")) {
    cache <- y
    candidate <- z %||% d
  } else {
    cache <- prepare_model(y, d, X)
    candidate <- z
  }
  if (is.null(candidate)) {
    stop("z is required", call. = FALSE)
  }
  vcov <- .rv_match_vcov(vcov)
  cluster <- .rv_resolve_cluster(cluster, cache)
  if (is.null(dim(candidate))) {
    .rv_evaluate_single_fast(cache, candidate, vcov, cluster, null)
  } else {
    result <- evaluate_subspace(cache, candidate, vcov = vcov,
                                cluster = cluster,
                                original_rank = ncol(as.matrix(candidate)),
                                null = null)
    result$rho_d <- NA_real_
    result$rho_y <- NA_real_
    result
  }
}

.rv_cluster_matrix_meat <- function(W, e, cluster) {
  info <- .rv_cluster_codes(cluster, nrow(W))
  scores <- rowsum(W * e, info$codes, reorder = FALSE)
  list(meat = crossprod(scores), G = info$G)
}

.rv_evaluate_reference_cache <- function(cache, z, vcov, cluster, null) {
  Z <- .rv_as_numeric_matrix(z, "z")
  if (nrow(Z) != cache$n) {
    stop("z has the wrong number of rows", call. = FALSE)
  }
  U <- .rv_qr_basis(.rv_residualize(cache$Qx, Z))
  k <- ncol(U)
  if (k == 0L) {
    stop("z is absorbed by X", call. = FALSE)
  }
  W <- cbind(cache$Qx, cache$d, U)
  fit <- stats::lm.fit(W, cache$y)
  if (fit$rank != ncol(W)) {
    stop("the full regression is rank deficient", call. = FALSE)
  }
  treatment_index <- cache$rank_x + 1L
  tau <- unname(fit$coefficients[treatment_index])
  e <- as.numeric(fit$residuals)
  bread <- solve(crossprod(W))
  df_full <- cache$n - ncol(W)
  variance_classical <- .rv_dot(e, e) / df_full *
    bread[treatment_index, treatment_index]
  if (vcov == "classical") {
    V <- .rv_dot(e, e) / df_full * bread
    correction <- 1
    n_groups <- cache$n
    meat_scalar <- NA_real_
  } else if (vcov %in% c("HC0", "HC1")) {
    meat_matrix <- crossprod(W, W * e^2)
    correction <- if (vcov == "HC1") cache$n / df_full else 1
    V <- correction * bread %*% meat_matrix %*% bread
    n_groups <- cache$n
    q <- as.numeric(.rv_residualize(U, cache$d))
    meat_scalar <- sum((q * e)^2)
  } else {
    cluster_meat <- .rv_cluster_matrix_meat(W, e, cluster)
    correction <- if (vcov == "CR1") {
      cluster_meat$G / (cluster_meat$G - 1) * (cache$n - 1) / df_full
    } else {
      1
    }
    V <- correction * bread %*% cluster_meat$meat %*% bread
    n_groups <- cluster_meat$G
    q <- as.numeric(.rv_residualize(U, cache$d))
    meat_scalar <- .rv_cluster_scalar_meat(q, e, cluster)$meat
  }
  variance_robust <- V[treatment_index, treatment_index]
  se_classical <- sqrt(variance_classical)
  se_robust <- sqrt(variance_robust)
  q <- as.numeric(.rv_residualize(U, cache$d))
  ee <- .rv_dot(cache$e_short, cache$e_short)
  result <- structure(list(
    tau = tau,
    null = null,
    se_classical = se_classical,
    se_robust = se_robust,
    variance_classical = variance_classical,
    variance_robust = variance_robust,
    omega = se_robust / se_classical,
    t_classical = (tau - null) / se_classical,
    t_robust = (tau - null) / se_robust,
    q = q,
    residuals = e,
    confounder_basis = U,
    r2_d = min(max(1 - .rv_dot(q, q) / .rv_dot(cache$d, cache$d), 0), 1),
    r2_y = if (ee > 0) min(max(1 - .rv_dot(e, e) / ee, 0), 1) else 0,
    geometric_rank = k,
    original_rank = k,
    df_residual = df_full,
    vcov = vcov,
    finite_sample_correction = correction,
    meat = meat_scalar,
    n_groups = n_groups
  ), class = "rv_confounder_evaluation")
  if (k == 1L) {
    correlations <- .rv_signed_partial_correlations(cache, U[, 1L])
    result$rho_d <- correlations[["rho_d"]]
    result$rho_y <- correlations[["rho_y"]]
  } else {
    result$rho_d <- NA_real_
    result$rho_y <- NA_real_
  }
  result
}

#' Independent explicit full-regression evaluation
#'
#' This reference calculation fits the complete regression and forms the full
#' bread-meat-bread covariance matrix. It is deliberately separate from the
#' optimized scalar FWL evaluator.
#'
#' @inheritParams evaluate_confounder
#' @export
evaluate_confounder_reference <- function(y, d = NULL, X = NULL, z = NULL,
                                          vcov = "HC1", cluster = NULL,
                                          null = 0) {
  if (inherits(y, "rv_cache")) {
    cache <- y
    candidate <- z %||% d
  } else {
    cache <- prepare_model(y, d, X)
    candidate <- z
  }
  if (is.null(candidate)) {
    stop("z is required", call. = FALSE)
  }
  vcov <- .rv_match_vcov(vcov)
  cluster <- .rv_resolve_cluster(cluster, cache)
  .rv_evaluate_reference_cache(cache, candidate, vcov, cluster, null)
}

.rv_short_evaluation <- function(cache, vcov = "HC1", cluster = NULL,
                                 null = 0) {
  evaluate_subspace(cache, matrix(numeric(), nrow = cache$n, ncol = 0L),
                    vcov = vcov, cluster = cluster, original_rank = 0L,
                    null = null)
}

#' @export
print.rv_confounder_evaluation <- function(x, ...) {
  cat("Omitted-confounder evaluation\n")
  cat("  Covariance:", x$vcov, "\n")
  cat("  Partial R2 (treatment, outcome):",
      format(x$r2_d, digits = 6), ",", format(x$r2_y, digits = 6), "\n")
  cat("  Coefficient:", format(x$tau, digits = 7), "\n")
  cat("  Classical / robust SE:", format(x$se_classical, digits = 7), "/",
      format(x$se_robust, digits = 7), "\n")
  cat("  Robust t:", format(x$t_robust, digits = 7), "\n")
  invisible(x)
}
