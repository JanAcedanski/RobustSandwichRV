#' @export
print.rv_benchmark <- function(x, ...) {
  cat("Observed-control benchmark calibration\n")
  cat("  Name:", x$benchmark_name, "\n")
  cat("  Type / rank:", x$benchmark_type, "/", x$group_rank, "\n")
  cat("  Signed per-rank delta_L:",
      format(x$delta_L_signed, digits = 7), "\n")
  cat("  Bounded magnitude b:", format(x$b, digits = 7), "\n")
  invisible(x)
}

#' Print a robustness-value result
#'
#' @param x A `robust_rv_result`.
#' @param ... Unused.
#' @export
print.robust_rv_result <- function(x, ...) {
  title <- if (identical(x$calibration_mode, "bounded")) {
    "Benchmark-bounded robust RV"
  } else {
    "Point-conditioned robust RV"
  }
  cat(title, "\n")
  cat("  Covariance:", x$vcov, "\n")
  cat("  Baseline decision:", x$baseline_decision, "\n")
  cat("  Status:", x$status, "/", x$substatus, "\n")
  if (!is.null(x$certificate_type)) {
    cat("  Certificate:", x$certificate_type, "\n")
  }
  cat("  Certified lower bound:", format(x$RV_lower, digits = 7), "\n")
  cat("  Certified upper bound:",
      if (is.finite(x$RV_upper)) format(x$RV_upper, digits = 7) else "not found",
      "\n")
  if (isTRUE(x$exact)) {
    cat("  Reported RV:", format(x$rv, digits = 7), "\n")
  }
  if (!is.null(x$best_witness_omega) && is.finite(x$best_witness_omega)) {
    cat("  Witness omega / robust t:",
        format(x$best_witness_omega, digits = 7), "/",
        format(x$best_witness_t, digits = 7), "\n")
  }
  cat("  Runtime:", format(x$runtime, digits = 4), "seconds\n")
  invisible(x)
}

#' Summarize a robustness-value result
#'
#' @param object A `robust_rv_result`.
#' @param ... Unused.
#' @export
summary.robust_rv_result <- function(object, ...) {
  table <- data.frame(
    calibration_mode = object$calibration_mode,
    vcov = object$vcov,
    status = object$status,
    substatus = object$substatus,
    baseline_decision = object$baseline_decision,
    RV_lower = object$RV_lower,
    RV_upper = object$RV_upper,
    RV_gap = object$RV_gap,
    exact = object$exact,
    theorem_exact = object$theorem_exact %||% NA,
    numerically_closed = object$numerically_closed %||% NA,
    certificate_type = object$certificate_type %||% NA_character_,
    omega_lower = object$omega_lower %||% object$omega_target,
    omega_upper = object$omega_upper %||% object$omega_target,
    witness_omega = object$best_witness_omega %||% NA_real_,
    witness_t = object$best_witness_t %||% NA_real_,
    runtime = object$runtime,
    objective_evals = object$objective_evals,
    stringsAsFactors = FALSE
  )
  structure(list(call = object$model_call, coefficients = table,
                 result = object), class = "summary.robust_rv_result")
}

#' @export
print.summary.robust_rv_result <- function(x, ...) {
  cat("Certified robustness-value summary\n")
  print(x$coefficients, row.names = FALSE)
  invisible(x)
}

#' @export
as.data.frame.robust_rv_result <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  summary(x)$coefficients
}

#' Tidy one robustness-value result
#'
#' @param x A `robust_rv_result`.
#' @return A one-row data frame.
#' @export
tidy_rv <- function(x) {
  if (!inherits(x, "robust_rv_result")) {
    stop("x must be a robust_rv_result", call. = FALSE)
  }
  as.data.frame(x)
}

#' Plot the analytical decision geometry
#'
#' Shows the fixed-omega CH decision region at the adverse endpoint (bounded
#' mode) or target omega (point mode), together with the protected solution and
#' the best verified witness.
#'
#' @param x A `robust_rv_result`.
#' @param cache The `rv_cache` used to compute `x`.
#' @param maximum Largest displayed partial R-squared.
#' @param grid_size Grid resolution.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @export
plot.robust_rv_result <- function(x, cache = NULL, maximum = NULL,
                                  grid_size = 100L, ...) {
  if (is.null(cache) || !inherits(cache, "rv_cache")) {
    stop("plotting decision geometry requires the original rv_cache",
         call. = FALSE)
  }
  target <- x$omega_adverse %||% x$omega_target
  if (!is.finite(target)) stop("result has no omega target", call. = FALSE)
  radius <- max(c(x$RV_upper, x$RV_lower, 0.05), na.rm = TRUE)
  maximum <- maximum %||% min(0.95, max(0.1, 1.35 * radius))
  values <- seq(0, maximum, length.out = as.integer(grid_size))
  changed <- matrix(FALSE, nrow = length(values), ncol = length(values))
  reject_short <- identical(x$baseline_decision, "reject")
  for (i in seq_along(values)) {
    for (j in seq_along(values)) {
      changed[i, j] <- any(vapply(c(-1, 1), function(sign) {
        t <- .rv_classical_ch_t_signed(cache, values[[i]], values[[j]], sign)
        if (reject_short) abs(t / target) <= x$critical_value else
          abs(t / target) >= x$critical_value
      }, logical(1)))
    }
  }
  graphics::image(values, values, changed,
                  col = c("#dbeafe", "#fecaca"),
                  xlab = expression(R[D]^2), ylab = expression(R[Y]^2), ...)
  graphics::abline(0, 1, lty = 2, col = "grey35")
  graphics::points(x$x_prot, x$y_prot, pch = 21, bg = "gold", cex = 1.2)
  witness_x <- x$best_witness_x %||% x$x_star
  witness_y <- x$best_witness_y %||% x$y_star
  if (!is.null(witness_x) && is.finite(witness_x) &&
      !is.null(witness_y) && is.finite(witness_y)) {
    graphics::points(witness_x, witness_y, pch = 23, bg = "black", cex = 1.2)
  }
  graphics::legend("topright",
                   legend = c("decision preserved", "decision changed",
                              "protected optimum", "verified witness"),
                   fill = c("#dbeafe", "#fecaca", "gold", "black"),
                   bty = "n", cex = 0.8)
  invisible(x)
}
