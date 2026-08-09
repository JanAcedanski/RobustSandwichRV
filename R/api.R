#' Calibrated robustness value
#'
#' This is the cache-level dispatcher. Bounded calibration is the recommended
#' default; point conditioning remains available without reinterpretation.
#'
#' @param cache An `rv_cache`.
#' @param benchmark Benchmark object or `b` in bounded mode; omega target in
#'   point mode when `point_omega` is not supplied.
#' @param calibration `"bounded"` or `"point"`.
#' @param point_omega Explicit point-conditioned omega target.
#' @param point_method Point solver method.
#' @param ... Arguments passed to [rv_cal_bounded()] or
#'   [conditional_omega_rv()].
#' @export
calibrated_rv <- function(cache, benchmark,
                          calibration = c("bounded", "point"),
                          point_omega = NULL,
                          point_method = "protected_first", ...) {
  calibration <- match.arg(calibration)
  if (calibration == "bounded") {
    return(rv_cal_bounded(cache, benchmark, ...))
  }
  target <- point_omega
  if (is.null(target)) {
    if (.rv_is_scalar_number(benchmark)) {
      target <- benchmark
    } else {
      target <- benchmark$omega_target %||% benchmark$omega_benchmark
    }
  }
  if (!.rv_is_scalar_number(target) || target <= 0) {
    stop("point calibration requires a positive omega target", call. = FALSE)
  }
  conditional_omega_rv(cache, target, method = point_method, ...)
}

#' Main model-facing robustness-value workflow
#'
#' @param model An unweighted `lm` object or an `rv_cache`.
#' @param treatment Name of the scalar treatment model-matrix column. Not used
#'   when `model` is already an `rv_cache`.
#' @param vcov HC0, HC1, CR0, or CR1.
#' @param calibration `"bounded"` (recommended) or `"point"`.
#' @param benchmark In bounded mode, an observed control term, vector of terms,
#'   benchmark object, or direct nonnegative `b`. In point mode, an omega target
#'   or object containing `omega_target`.
#' @param b Optional direct bounded benchmark magnitude.
#' @param omega_target Optional point-conditioned ratio.
#' @param cluster Optional cluster formula or labels.
#' @param ... Numerical solver options.
#' @return A `robust_rv_result`.
#' @export
robust_rv <- function(model, treatment = NULL, vcov = "HC1",
                      calibration = c("bounded", "point"),
                      benchmark = NULL, b = NULL, omega_target = NULL,
                      cluster = NULL, ...) {
  calibration <- match.arg(calibration)
  vcov <- match.arg(vcov, c("HC0", "HC1", "CR0", "CR1"))
  if (inherits(model, "rv_cache")) {
    cache <- model
    original_model <- NULL
  } else {
    cache <- cache_from_lm(model, treatment)
    original_model <- model
  }
  if (calibration == "bounded") {
    if (!is.null(b)) {
      calibration_object <- b
    } else if (inherits(benchmark, "rv_benchmark") ||
               (.rv_is_scalar_number(benchmark) && benchmark >= 0)) {
      calibration_object <- benchmark
    } else if (!is.null(original_model) && !is.null(benchmark)) {
      calibration_object <- benchmark_from_model(
        original_model, treatment, benchmark, vcov = vcov, cluster = cluster
      )
    } else {
      stop("bounded calibration requires b, an rv_benchmark, or an observed-control benchmark",
           call. = FALSE)
    }
    result <- rv_cal_bounded(
      cache, calibration_object, vcov = vcov, cluster = cluster, ...
    )
    result$model_call <- if (is.null(original_model)) NULL else original_model$call
    result$treatment <- treatment %||% cache$treatment
    result$benchmark_object <- calibration_object
    return(result)
  }
  target <- omega_target
  if (is.null(target)) {
    if (.rv_is_scalar_number(benchmark)) {
      target <- benchmark
    } else if (is.list(benchmark)) {
      target <- benchmark$omega_target %||% benchmark$omega_benchmark
    }
  }
  if (!.rv_is_scalar_number(target) || target <= 0) {
    stop("point calibration requires omega_target", call. = FALSE)
  }
  result <- conditional_omega_rv(
    cache, target, vcov = vcov, cluster = cluster, ...
  )
  result$model_call <- if (is.null(original_model)) NULL else original_model$call
  result$treatment <- treatment %||% cache$treatment
  result
}

#' @rdname robust_rv
#' @export
observed_control_rv <- robust_rv

#' Validate a short-model covariance against `sandwich`
#'
#' @param model An `lm` model.
#' @param treatment Treatment coefficient name.
#' @param vcov HC0, HC1, CR0, or CR1.
#' @param cluster Cluster labels or formula for CR estimators.
#' @return A data frame comparing standard errors and t statistics.
#' @export
validate_sandwich_reference <- function(model, treatment, vcov = "HC1",
                                        cluster = NULL) {
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("package 'sandwich' is required for this validation", call. = FALSE)
  }
  vcov <- match.arg(vcov, c("HC0", "HC1", "CR0", "CR1"))
  cache <- cache_from_lm(model, treatment)
  labels <- .rv_resolve_cluster(cluster, cache)
  internal <- .rv_short_evaluation(cache, vcov, labels)
  if (vcov %in% c("HC0", "HC1")) {
    external_vcov <- sandwich::vcovHC(model, type = vcov)
  } else {
    external_vcov <- sandwich::vcovCL(
      model, cluster = labels,
      type = if (vcov == "CR1") "HC1" else "HC0",
      cadjust = vcov == "CR1"
    )
  }
  external_se <- sqrt(external_vcov[treatment, treatment])
  estimate <- stats::coef(model)[[treatment]]
  data.frame(
    implementation = c("RobustSandwichRV", "sandwich"),
    standard_error = c(internal$se_robust, external_se),
    t_statistic = c(internal$t_robust, estimate / external_se),
    row.names = NULL
  )
}
