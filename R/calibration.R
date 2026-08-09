#' Construct a rank-normalized benchmark calibration
#'
#' @param delta_L Signed per-synthetic-index change in log square-root score
#'   alignment.
#' @param benchmark_type `"scalar"` or `"group"`.
#' @param benchmark_name Human-readable benchmark name.
#' @param group_rank Numerical rank of the observed benchmark block.
#' @param rank_normalized Whether `delta_L` is already divided by group rank.
#' @param L_full,L_reduced Optional score-alignment values. When supplied,
#'   `delta_L` is computed as `log(L_full / L_reduced) / (2 * group_rank)`.
#' @return An `rv_benchmark` object.
#' @export
benchmark_calibration <- function(delta_L = NULL,
                                  benchmark_type = c("scalar", "group"),
                                  benchmark_name = "benchmark",
                                  group_rank = 1L,
                                  rank_normalized = TRUE,
                                  L_full = NULL, L_reduced = NULL) {
  benchmark_type <- match.arg(benchmark_type)
  group_rank <- as.integer(group_rank)
  if (group_rank < 1L) {
    stop("group_rank must be positive", call. = FALSE)
  }
  if (!is.null(L_full) || !is.null(L_reduced)) {
    if (!.rv_is_scalar_number(L_full) || !.rv_is_scalar_number(L_reduced) ||
        L_full <= 0 || L_reduced <= 0) {
      stop("L_full and L_reduced must be positive and finite", call. = FALSE)
    }
    delta_L <- 0.5 * log(L_full / L_reduced) / group_rank
    rank_normalized <- TRUE
  }
  if (!.rv_is_scalar_number(delta_L)) {
    stop("a finite signed delta_L is required", call. = FALSE)
  }
  if (!isTRUE(rank_normalized)) {
    delta_L <- delta_L / group_rank
  }
  structure(list(
    benchmark_type = benchmark_type,
    benchmark_name = as.character(benchmark_name),
    group_rank = group_rank,
    delta_L_signed = as.numeric(delta_L),
    b = abs(as.numeric(delta_L)),
    b_absolute = abs(as.numeric(delta_L)),
    rank_normalized = TRUE,
    L_full = L_full,
    L_reduced = L_reduced
  ), class = "rv_benchmark")
}

.rv_score_alignment <- function(cache, vcov, cluster = NULL) {
  d_hat <- cache$d / .rv_norm(cache$d)
  e_hat <- cache$e_short / .rv_norm(cache$e_short)
  score <- d_hat * e_hat
  if (vcov %in% c("HC0", "HC1")) {
    return(sum(score^2))
  }
  info <- .rv_cluster_codes(cluster, cache$n)
  totals <- rowsum(score, info$codes, reorder = FALSE)
  sum(totals^2)
}

.rv_omega_origin <- function(cache, vcov, cluster = NULL) {
  L <- .rv_score_alignment(cache, vcov, cluster)
  if (!(L > 0)) {
    stop("the zero-strength sandwich origin is degenerate", call. = FALSE)
  }
  scale2 <- switch(vcov,
    HC0 = cache$df_short - 1L,
    HC1 = cache$n,
    CR0 = cache$df_short - 1L,
    CR1 = {
      G <- .rv_cluster_codes(cluster, cache$n)$G
      G / (G - 1) * (cache$n - 1)
    },
    stop("unsupported covariance estimator", call. = FALSE)
  )
  sqrt(scale2 * L)
}

.rv_benchmark_column_indices <- function(extracted, benchmark) {
  benchmark <- as.character(benchmark)
  if (length(benchmark) < 1L) {
    stop("benchmark must identify at least one observed control", call. = FALSE)
  }
  mm_names <- colnames(extracted$model_matrix)
  selected <- integer()
  for (name in benchmark) {
    column_hit <- which(mm_names == name)
    term_hit <- which(extracted$term_labels == name)
    term_columns <- if (length(term_hit)) {
      which(extracted$assign %in% term_hit)
    } else {
      integer()
    }
    hits <- union(column_hit, term_columns)
    if (length(hits) == 0L) {
      stop("benchmark '", name, "' is neither a model term nor a matrix column",
           call. = FALSE)
    }
    selected <- union(selected, hits)
  }
  if (extracted$treatment_index %in% selected) {
    stop("the treatment cannot be used as its own benchmark", call. = FALSE)
  }
  selected
}

#' Calibrate an observed scalar control or control group
#'
#' The observed group is deleted from the control matrix, its numerical rank
#' is measured, and the score-alignment change is divided by that rank before
#' it is transported to the one-index sensitivity problem.
#'
#' @param model An unweighted `lm` model.
#' @param treatment Name of the treatment model-matrix column.
#' @param benchmark One or more model term or model-matrix column names.
#' @param vcov HC0, HC1, CR0, or CR1.
#' @param cluster Optional cluster formula or labels.
#' @return An `rv_benchmark` with full calibration metadata.
#' @export
benchmark_from_model <- function(model, treatment, benchmark,
                                 vcov = "HC1", cluster = NULL) {
  vcov <- match.arg(vcov, c("HC0", "HC1", "CR0", "CR1"))
  extracted <- .rv_extract_lm(model, treatment)
  full_cache <- cache_from_lm(model, treatment)
  cluster <- .rv_resolve_cluster(cluster, full_cache)
  selected_mm <- .rv_benchmark_column_indices(extracted, benchmark)
  control_mm <- setdiff(seq_len(ncol(extracted$model_matrix)),
                        extracted$treatment_index)
  selected_control <- match(selected_mm, control_mm)
  if (anyNA(selected_control)) {
    stop("internal benchmark column mapping failed", call. = FALSE)
  }
  X_reduced <- extracted$X[, -selected_control, drop = FALSE]
  reduced_cache <- prepare_model(extracted$y, extracted$d, X_reduced)
  group_rank <- full_cache$rank_x - reduced_cache$rank_x
  if (group_rank < 1L) {
    stop("the selected benchmark has zero numerical rank", call. = FALSE)
  }
  L_full <- .rv_score_alignment(full_cache, vcov, cluster)
  L_reduced <- .rv_score_alignment(reduced_cache, vcov, cluster)
  calibration <- benchmark_calibration(
    benchmark_type = if (group_rank == 1L && length(selected_mm) == 1L)
      "scalar" else "group",
    benchmark_name = paste(benchmark, collapse = " + "),
    group_rank = group_rank,
    L_full = L_full,
    L_reduced = L_reduced
  )
  calibration$selected_columns <- colnames(extracted$model_matrix)[selected_mm]
  calibration$omega_origin <- .rv_omega_origin(full_cache, vcov, cluster)
  calibration$omega_target <- calibration$omega_origin *
    exp(calibration$delta_L_signed)
  calibration$vcov <- vcov
  calibration
}

#' Construct a benchmark sandwich-ratio band
#'
#' @param omega_origin Positive zero-strength synthetic-index ratio.
#' @param b Nonnegative log-scale benchmark magnitude.
#' @param bound_direction `"symmetric"` or `"observed"`.
#' @param delta_L_signed Signed observed change, required for the observed
#'   direction mode.
#' @export
bounded_benchmark_band <- function(omega_origin, b,
                                   bound_direction = c("symmetric", "observed"),
                                   delta_L_signed = NA_real_) {
  bound_direction <- match.arg(bound_direction)
  if (!.rv_is_scalar_number(omega_origin) || omega_origin <= 0) {
    stop("omega_origin must be positive and finite", call. = FALSE)
  }
  if (!.rv_is_scalar_number(b) || b < 0) {
    stop("b must be finite and nonnegative", call. = FALSE)
  }
  if (bound_direction == "symmetric") {
    lower <- omega_origin * exp(-b)
    upper <- omega_origin * exp(b)
  } else {
    if (!.rv_is_scalar_number(delta_L_signed) ||
        abs(abs(delta_L_signed) - b) > max(1e-12, 1e-10 * b)) {
      stop("observed direction requires b = abs(delta_L_signed)",
           call. = FALSE)
    }
    if (delta_L_signed >= 0) {
      lower <- omega_origin
      upper <- omega_origin * exp(b)
    } else {
      lower <- omega_origin * exp(-b)
      upper <- omega_origin
    }
  }
  list(
    omega_origin = omega_origin,
    omega_lower = lower,
    omega_upper = upper,
    bound_direction = bound_direction
  )
}

.rv_coerce_benchmark <- function(benchmark) {
  if (inherits(benchmark, "rv_benchmark")) return(benchmark)
  if (.rv_is_scalar_number(benchmark) && benchmark >= 0) {
    return(structure(list(
      benchmark_type = "direct",
      benchmark_name = "direct b",
      group_rank = 1L,
      delta_L_signed = NA_real_,
      b = as.numeric(benchmark),
      b_absolute = as.numeric(benchmark),
      rank_normalized = TRUE
    ), class = "rv_benchmark"))
  }
  if (is.list(benchmark)) {
    delta <- benchmark$delta_L_signed %||% benchmark$delta_L
    b <- benchmark$b_absolute %||% benchmark$b
    rank <- as.integer(benchmark$group_rank %||% 1L)
    if (.rv_is_scalar_number(delta)) {
      return(benchmark_calibration(
        delta_L = delta,
        benchmark_type = benchmark$benchmark_type %||%
          if (rank == 1L) "scalar" else "group",
        benchmark_name = benchmark$benchmark_name %||% "benchmark",
        group_rank = rank,
        rank_normalized = benchmark$rank_normalized %||% TRUE
      ))
    }
    if (.rv_is_scalar_number(b) && b >= 0) {
      return(structure(list(
        benchmark_type = benchmark$benchmark_type %||% "direct",
        benchmark_name = benchmark$benchmark_name %||% "benchmark",
        group_rank = rank,
        delta_L_signed = NA_real_, b = b, b_absolute = b,
        rank_normalized = TRUE
      ), class = "rv_benchmark"))
    }
  }
  stop("benchmark must be a nonnegative b or an rv_benchmark", call. = FALSE)
}

.rv_validate_bounded_z <- function(cache, z, vcov, cluster, null,
                                   short_t, critical, omega_lower, omega_upper,
                                   band_tolerance = 5e-7,
                                   decision_tolerance = 5e-7) {
  attained <- tryCatch(
    evaluate_confounder_reference(cache, z, vcov = vcov,
                                  cluster = cluster, null = null),
    error = function(e) e
  )
  if (inherits(attained, "error")) {
    return(list(valid = FALSE, closure_valid = FALSE,
                literal_change = FALSE, error = conditionMessage(attained)))
  }
  band_violation <- max(log(omega_lower / attained$omega),
                        log(attained$omega / omega_upper), 0)
  reject_short <- .rv_reject(short_t, critical)
  decision_violation <- .rv_decision_violation(
    reject_short, attained$t_robust, critical
  )
  literal <- .rv_decision_change(short_t, attained$t_robust, critical)
  closure <- band_violation <= band_tolerance &&
    decision_violation <= decision_tolerance
  list(
    valid = closure && literal,
    closure_valid = closure,
    literal_change = literal,
    attained = attained,
    radius = max(attained$r2_d, attained$r2_y),
    band_violation = band_violation,
    decision_violation = decision_violation,
    z = as.numeric(z),
    error = NULL
  )
}

.rv_fixed_xy_bounded_witness <- function(cache, x, y, bias_sign,
                                         omega_evaluator, vcov, cluster,
                                         null, short_t, critical,
                                         omega_lower, omega_upper,
                                         n_starts, seed, max_iterations,
                                         root_tolerance, warm = list()) {
  interval <- .rv_direction_interval(
    cache, x, y, bias_sign, omega_evaluator,
    n_starts = n_starts, seed = seed,
    max_iterations = max_iterations, warm = warm
  )
  t_classical <- .rv_classical_ch_t_signed(cache, x, y, bias_sign, null)
  threshold <- abs(t_classical) / critical
  reject_short <- .rv_reject(short_t, critical)
  safety <- max(1e-9, 10 * root_tolerance)
  if (reject_short) {
    feasible_lower <- max(omega_lower, threshold * exp(safety))
    feasible_upper <- omega_upper
  } else {
    feasible_lower <- omega_lower
    feasible_upper <- min(omega_upper, threshold * exp(-safety))
  }
  intersection_lower <- max(feasible_lower, interval$omega_low)
  intersection_upper <- min(feasible_upper, interval$omega_high)
  if (!(intersection_lower < intersection_upper)) {
    return(list(valid = FALSE, evaluations = interval$evaluations,
                interval = interval, reason = "empty_interval_intersection"))
  }
  target <- sqrt(intersection_lower * intersection_upper)
  root <- .rv_omega_root(
    cache, x, y, bias_sign, interval, target, omega_evaluator,
    seed = seed + 300007L, root_tolerance = root_tolerance
  )
  evaluations <- interval$evaluations + root$evaluations
  if (!isTRUE(root$found)) {
    return(list(valid = FALSE, evaluations = evaluations,
                interval = interval, root = root,
                reason = "bounded_root_failure"))
  }
  z <- single_parameterization(cache, sqrt(x), bias_sign * sqrt(y), root$u)
  validation <- .rv_validate_bounded_z(
    cache, z, vcov, cluster, null, short_t, critical,
    omega_lower, omega_upper
  )
  c(validation, list(
    evaluations = evaluations,
    interval = interval,
    root = root,
    target_omega = target,
    x = x,
    y = y,
    bias_sign = bias_sign,
    u = root$u
  ))
}

.rv_bounded_candidate_at_radius <- function(cache, radius, protected,
                                            omega_evaluator, vcov, cluster,
                                            null, short_t, critical,
                                            omega_lower, omega_upper,
                                            n_starts, seed, max_iterations,
                                            root_tolerance,
                                            n_random = 8L, warm = list()) {
  reject_short <- .rv_reject(short_t, critical)
  signs <- .rv_protected_signs(cache, reject_short, null)
  pairs <- .rv_candidate_pairs(cache, radius, protected$branch,
                               n_random, seed)
  best <- NULL
  evaluations <- 0L
  for (pair_index in seq_len(nrow(pairs))) {
    for (sign_index in seq_along(signs)) {
      witness <- .rv_fixed_xy_bounded_witness(
        cache, pairs[pair_index, 1L], pairs[pair_index, 2L],
        signs[[sign_index]], omega_evaluator, vcov, cluster,
        null, short_t, critical, omega_lower, omega_upper,
        n_starts = n_starts,
        seed = seed + 1013L * pair_index + 10007L * sign_index,
        max_iterations = max_iterations,
        root_tolerance = root_tolerance,
        warm = warm
      )
      evaluations <- evaluations + witness$evaluations
      if (isTRUE(witness$valid) &&
          (is.null(best) || witness$radius < best$radius)) {
        best <- witness
      }
    }
  }
  list(found = !is.null(best), witness = best, evaluations = evaluations)
}

.rv_search_bounded_upper <- function(cache, protected, omega_evaluator,
                                     vcov, cluster, null, short_t, critical,
                                     omega_lower, omega_upper, lower,
                                     n_starts, seed, max_iterations,
                                     root_tolerance, max_radius = 0.9,
                                     radial_points = 18L, warm = list()) {
  radii <- unique(c(
    lower,
    pmin(max_radius, lower + c(1e-7, 1e-5, 1e-3, 1e-2)),
    seq(lower, max_radius, length.out = radial_points)
  ))
  radii <- sort(radii[radii >= 0 & radii < 1])
  previous <- lower
  best <- NULL
  evaluations <- 0L
  for (index in seq_along(radii)) {
    starts <- c(warm, if (is.null(best)) list() else list(best$u))
    candidate <- .rv_bounded_candidate_at_radius(
      cache, radii[[index]], protected, omega_evaluator, vcov, cluster,
      null, short_t, critical, omega_lower, omega_upper,
      n_starts = n_starts, seed = seed + 100003L * index,
      max_iterations = max_iterations, root_tolerance = root_tolerance,
      n_random = max(6L, n_starts %/% 2L), warm = starts
    )
    evaluations <- evaluations + candidate$evaluations
    if (candidate$found) {
      best <- candidate$witness
      lo <- previous
      hi <- radii[[index]]
      for (step in seq_len(10L)) {
        mid <- (lo + hi) / 2
        refined <- .rv_bounded_candidate_at_radius(
          cache, mid, protected, omega_evaluator, vcov, cluster,
          null, short_t, critical, omega_lower, omega_upper,
          n_starts = n_starts,
          seed = seed + 800011L + step,
          max_iterations = max_iterations,
          root_tolerance = root_tolerance,
          n_random = 4L, warm = list(best$u)
        )
        evaluations <- evaluations + refined$evaluations
        if (refined$found) {
          hi <- mid
          if (refined$witness$radius < best$radius) best <- refined$witness
        } else {
          lo <- mid
        }
      }
      break
    }
    previous <- radii[[index]]
  }
  list(found = !is.null(best), witness = best, evaluations = evaluations)
}

.rv_endpoint_perturbation <- function(cache, certificate, protected,
                                      omega_target, omega_evaluator,
                                      vcov, cluster, null, short_t, critical,
                                      seed, max_iterations, root_tolerance) {
  if (isTRUE(certificate$witness$literal_change)) {
    return(certificate$witness)
  }
  for (index in seq_along(c(2e-7, 2e-6, 2e-5, 1e-4, 5e-4, 2e-3))) {
    delta <- c(2e-7, 2e-6, 2e-5, 1e-4, 5e-4, 2e-3)[[index]]
    radius <- min(protected$rv + delta, 1 - 1e-7)
    pair <- .rv_branch_pair(cache, protected$branch, radius)
    if (any(!is.finite(pair)) || any(pair < 0) || any(pair >= 1)) next
    witness <- .rv_fixed_xy_point_witness(
      cache, pair[[1L]], pair[[2L]], certificate$witness$bias_sign,
      omega_target, omega_evaluator, vcov, cluster, null, short_t,
      critical, n_starts = 24L,
      seed = seed + 100003L * index,
      max_iterations = max_iterations,
      root_tolerance = root_tolerance,
      warm = list(certificate$witness$u),
      require_literal = TRUE
    )
    if (isTRUE(witness$valid)) return(witness)
  }
  NULL
}

.rv_validate_supplied_witnesses <- function(cache, witnesses, vcov, cluster,
                                            null, short_t, critical,
                                            omega_lower, omega_upper) {
  output <- list()
  for (candidate in witnesses) {
    z <- if (is.numeric(candidate)) candidate else candidate$z %||%
      candidate$best_witness_z
    if (is.null(z)) next
    validation <- .rv_validate_bounded_z(
      cache, z, vcov, cluster, null, short_t, critical,
      omega_lower, omega_upper
    )
    if (isTRUE(validation$valid)) output[[length(output) + 1L]] <- validation
  }
  output
}

#' Benchmark-bounded calibrated robustness value
#'
#' Computes the analytical adverse-endpoint lower bound and searches for an
#' explicit, independently refitted decision-changing confounder whose
#' sandwich ratio lies anywhere inside the benchmark band.
#'
#' @param cache An `rv_cache`.
#' @param benchmark A nonnegative `b` or an [benchmark_calibration()] object.
#' @param vcov HC0, HC1, CR0, or CR1.
#' @param cluster Optional one-way cluster labels.
#' @param bound_direction Symmetric or observed-direction band.
#' @param alpha Two-sided test level.
#' @param null Null coefficient value.
#' @param critical Optional fixed critical value.
#' @param n_starts Number of sphere starts.
#' @param seed Local random seed.
#' @param max_iterations Maximum local sphere iterations.
#' @param root_tolerance Log-omega root tolerance.
#' @param exact_tolerance Numerical interval-closing tolerance.
#' @param max_radius Largest searched partial-R-squared radius.
#' @param supplied_witnesses Optional earlier witness vectors for nested bands.
#' @param run_witness_search If `FALSE`, stop after the endpoint certificate.
#' @return A `robust_rv_result` with certified lower and upper endpoints.
#' @export
rv_cal_bounded <- function(cache, benchmark, vcov = "HC1", cluster = NULL,
                           bound_direction = c("symmetric", "observed"),
                           alpha = 0.05, null = 0, critical = NULL,
                           n_starts = 16L, seed = 1234L,
                           max_iterations = 60L,
                           root_tolerance = 1e-9,
                           exact_tolerance = 2e-4,
                           max_radius = 0.9,
                           supplied_witnesses = list(),
                           run_witness_search = TRUE) {
  started <- .rv_timer()
  vcov <- match.arg(vcov, c("HC0", "HC1", "CR0", "CR1"))
  bound_direction <- match.arg(bound_direction)
  cluster <- .rv_resolve_cluster(cluster, cache)
  calibration <- .rv_coerce_benchmark(benchmark)
  critical <- .rv_critical_value(cache, vcov, 1L, alpha, cluster, critical)
  short <- .rv_short_evaluation(cache, vcov, cluster, null)
  reject_short <- .rv_reject(short$t_robust, critical)
  omega_evaluator <- conditional_omega_evaluator(cache, vcov, cluster)
  omega_origin <- .rv_omega_origin(cache, vcov, cluster)
  band <- bounded_benchmark_band(
    omega_origin, calibration$b,
    bound_direction = bound_direction,
    delta_L_signed = calibration$delta_L_signed
  )
  omega_adverse <- if (reject_short) band$omega_upper else band$omega_lower
  RV_eq <- .rv_conditional_equal(cache, omega_adverse, critical,
                                 reject_short, null)
  protected <- .rv_conditional_protected(cache, omega_adverse, critical,
                                         reject_short, null)
  endpoint <- .rv_protected_certificate(
    cache, protected, omega_adverse, omega_evaluator, vcov, cluster,
    null, short$t_robust, critical,
    n_starts = n_starts, seed = seed,
    max_iterations = max_iterations,
    root_tolerance = root_tolerance,
    require_literal = FALSE
  )
  base <- list(
    calibration_mode = "bounded",
    bound_direction = bound_direction,
    benchmark_type = calibration$benchmark_type,
    benchmark_name = calibration$benchmark_name,
    group_rank = calibration$group_rank,
    delta_L_signed = calibration$delta_L_signed,
    b = calibration$b,
    omega_origin = omega_origin,
    omega_lower = band$omega_lower,
    omega_upper = band$omega_upper,
    omega_adverse = omega_adverse,
    baseline_decision = if (reject_short) "reject" else "do_not_reject",
    short_t_robust = short$t_robust,
    critical_value = critical,
    vcov = vcov,
    RV_eq_adverse = RV_eq,
    RV_prot_adverse = protected$rv,
    protected_regime = .rv_protected_regime(protected$branch),
    protected_branch = protected$branch,
    x_prot = protected$x,
    y_prot = protected$y,
    endpoint_certificate_attempted = TRUE,
    endpoint_certificate_found = isTRUE(endpoint$found)
  )
  endpoint_witness <- NULL
  if (isTRUE(endpoint$found)) {
    endpoint_witness <- if (reject_short) {
      .rv_endpoint_perturbation(
        cache, endpoint, protected, omega_adverse, omega_evaluator,
        vcov, cluster, null, short$t_robust, critical,
        seed + 500009L, max_iterations, root_tolerance
      )
    } else if (isTRUE(endpoint$witness$literal_change)) {
      endpoint$witness
    } else {
      NULL
    }
  }
  if (!is.null(endpoint_witness)) {
    return(structure(c(base, list(
      rv = protected$rv,
      RV_cal_bounded = protected$rv,
      RV_lower = protected$rv,
      RV_upper = protected$rv,
      RV_gap = 0,
      best_witness_x = endpoint_witness$attained$r2_d,
      best_witness_y = endpoint_witness$attained$r2_y,
      best_witness_s = endpoint_witness$bias_sign,
      best_witness_u = endpoint_witness$u,
      best_witness_z = endpoint_witness$z,
      best_witness_omega = endpoint_witness$attained$omega,
      best_witness_t = endpoint_witness$attained$t_robust,
      direct_refit_valid = TRUE,
      witness_literal_decision_change = TRUE,
      theorem_exact = TRUE,
      numerically_closed = TRUE,
      certificate_type = "adverse_endpoint_equality",
      exact = TRUE,
      status = "EXACT_BOUNDED_RV",
      substatus = if (.rv_protected_regime(protected$branch) == "equal_strength")
        "ADVERSE_ENDPOINT_PROTECTED_DIAGONAL" else
        "ADVERSE_ENDPOINT_PROTECTED_INTERIOR",
      runtime = .rv_timer() - started,
      objective_evals = endpoint$evaluations,
      witness_count = 1L
    )), class = "robust_rv_result"))
  }
  verified <- .rv_validate_supplied_witnesses(
    cache, supplied_witnesses, vcov, cluster, null, short$t_robust,
    critical, band$omega_lower, band$omega_upper
  )
  search <- list(found = FALSE, evaluations = 0L)
  if (isTRUE(run_witness_search)) {
    warm <- lapply(verified, function(item) {
      correlations <- .rv_signed_partial_correlations(
        cache, item$attained$confounder_basis[, 1L]
      )
      z_hat <- item$attained$confounder_basis[, 1L]
      rd <- correlations[["rho_d"]]
      ry <- correlations[["rho_y"]]
      denominator <- sqrt(max(1 - rd^2, 0)) * sqrt(max(1 - ry^2, 0))
      if (denominator <= 1e-10) return(NULL)
      d_hat <- cache$d / .rv_norm(cache$d)
      e_hat <- cache$e_short / .rv_norm(cache$e_short)
      tryCatch(.rv_normalize_direction(
        cache, (z_hat - rd * d_hat - sqrt(1 - rd^2) * ry * e_hat) /
          denominator
      ), error = function(e) NULL)
    })
    warm <- Filter(Negate(is.null), warm)
    search <- .rv_search_bounded_upper(
      cache, protected, omega_evaluator, vcov, cluster, null,
      short$t_robust, critical, band$omega_lower, band$omega_upper,
      lower = protected$rv, n_starts = n_starts,
      seed = seed + 900001L, max_iterations = max_iterations,
      root_tolerance = root_tolerance, max_radius = max_radius,
      warm = warm
    )
  }
  candidates <- verified
  if (isTRUE(search$found)) candidates[[length(candidates) + 1L]] <- search$witness
  if (length(candidates) == 0L) {
    return(structure(c(base, list(
      rv = NA_real_, RV_cal_bounded = NA_real_,
      RV_lower = protected$rv, RV_upper = NA_real_, RV_gap = NA_real_,
      best_witness_z = NULL, direct_refit_valid = FALSE,
      witness_literal_decision_change = FALSE, exact = FALSE,
      theorem_exact = FALSE, numerically_closed = FALSE,
      certificate_type = "analytical_lower_bound_only",
      status = "LOWER_BOUND_ONLY",
      substatus = if (run_witness_search) "BOUNDED_WITNESS_NOT_FOUND" else
        "ENDPOINT_CERTIFICATE_FAILED",
      runtime = .rv_timer() - started,
      objective_evals = endpoint$evaluations + search$evaluations,
      witness_count = 0L
    )), class = "robust_rv_result"))
  }
  radii <- vapply(candidates, function(x) x$radius, numeric(1))
  best <- candidates[[which.min(radii)]]
  upper <- best$radius
  if (upper < protected$rv - exact_tolerance) {
    stop("verified upper bound lies below the analytical lower bound",
         call. = FALSE)
  }
  gap <- upper - protected$rv
  closed <- gap <= exact_tolerance
  structure(c(base, list(
    rv = if (closed) protected$rv else NA_real_,
    RV_cal_bounded = if (closed) protected$rv else NA_real_,
    RV_lower = protected$rv,
    RV_upper = upper,
    RV_gap = gap,
    best_witness_x = best$attained$r2_d,
    best_witness_y = best$attained$r2_y,
    best_witness_s = best$bias_sign %||% NA_integer_,
    best_witness_u = best$u %||% NULL,
    best_witness_z = best$z,
    best_witness_omega = best$attained$omega,
    best_witness_t = best$attained$t_robust,
    direct_refit_valid = TRUE,
    witness_literal_decision_change = TRUE,
    theorem_exact = FALSE,
    numerically_closed = closed,
    certificate_type = if (closed) "verified_interval_within_tolerance" else
      "verified_interval",
    exact = closed,
    status = if (closed) "EXACT_BOUNDED_RV" else "CERTIFIED_BOUNDS",
    substatus = if (closed) "BOUNDS_CLOSED_NUMERICALLY" else
      "BOUNDED_WITNESS_FOUND",
    runtime = .rv_timer() - started,
    objective_evals = endpoint$evaluations + search$evaluations,
    witness_count = length(candidates)
  )), class = "robust_rv_result")
}

#' Evaluate nested benchmark bands with witness reuse
#'
#' @param cache An `rv_cache`.
#' @param benchmarks A list of benchmark objects or magnitudes.
#' @param ... Arguments passed to [rv_cal_bounded()].
#' @export
rv_cal_bounded_sequence <- function(cache, benchmarks, ...) {
  calibrations <- lapply(benchmarks, .rv_coerce_benchmark)
  order <- order(vapply(calibrations, function(x) x$b, numeric(1)))
  output <- vector("list", length(benchmarks))
  witnesses <- list()
  previous_upper <- Inf
  for (index in order) {
    result <- rv_cal_bounded(cache, calibrations[[index]],
                             supplied_witnesses = witnesses, ...)
    if (!is.null(result$best_witness_z) &&
        isTRUE(result$witness_literal_decision_change)) {
      witnesses[[length(witnesses) + 1L]] <- result$best_witness_z
      if (is.finite(result$RV_upper)) {
        if (result$RV_upper > previous_upper + 1e-8) {
          stop("a nested-band upper bound worsened despite witness reuse",
               call. = FALSE)
        }
        previous_upper <- min(previous_upper, result$RV_upper)
      }
    }
    output[[index]] <- result
  }
  output
}
