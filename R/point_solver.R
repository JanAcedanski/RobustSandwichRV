.rv_direction_screen <- function(cache, x, y, bias_sign, omega_evaluator,
                                 n_starts, seed, warm = list()) {
  directions <- .rv_structured_directions(cache, n_starts, seed, warm)
  omega <- vapply(directions, function(u) {
    omega_evaluator(x, y, bias_sign, u)
  }, numeric(1))
  if (any(!is.finite(omega)) || any(omega <= 0)) {
    stop("omega evaluator returned a nonpositive or non-finite value",
         call. = FALSE)
  }
  low <- which.min(omega)
  high <- which.max(omega)
  list(
    u_low = directions[[low]],
    u_high = directions[[high]],
    omega_low = omega[[low]],
    omega_high = omega[[high]],
    evaluations = length(omega),
    starts = directions
  )
}

.rv_sphere_extreme <- function(cache, x, y, bias_sign, omega_evaluator,
                               start, sense = c("min", "max"), seed = 1234L,
                               max_iterations = 60L,
                               step_tolerance = 1e-6) {
  sense <- match.arg(sense)
  .rv_with_seed(seed, {
    u <- .rv_normalize_direction(cache, start)
    omega <- omega_evaluator(x, y, bias_sign, u)
    objective <- log(omega)
    evaluations <- 1L
    step <- 0.45
    converged <- FALSE
    iterations <- 0L
    is_better <- if (sense == "min") {
      function(new, old) new < old - 1e-12
    } else {
      function(new, old) new > old + 1e-12
    }
    for (iteration in seq_len(min(as.integer(max_iterations), 80L))) {
      iterations <- iteration
      candidate_u <- u
      candidate_omega <- omega
      candidate_objective <- objective
      improved <- FALSE
      for (probe_index in seq_len(6L)) {
        probe <- tryCatch(.rv_normalize_direction(cache, stats::rnorm(cache$n)),
                          error = function(e) NULL)
        if (is.null(probe)) next
        tangent <- probe - u * .rv_dot(u, probe)
        tangent_norm <- .rv_norm(tangent)
        if (!(tangent_norm > 1e-12)) next
        tangent <- tangent / tangent_norm
        for (orientation in c(-1, 1)) {
          proposed <- cos(step) * u + orientation * sin(step) * tangent
          proposed <- .rv_normalize_direction(cache, proposed)
          proposed_omega <- omega_evaluator(x, y, bias_sign, proposed)
          proposed_objective <- log(proposed_omega)
          evaluations <- evaluations + 1L
          if (is_better(proposed_objective, candidate_objective)) {
            candidate_u <- proposed
            candidate_omega <- proposed_omega
            candidate_objective <- proposed_objective
            improved <- TRUE
          }
        }
      }
      if (improved) {
        u <- candidate_u
        omega <- candidate_omega
        objective <- candidate_objective
        step <- min(0.65, 1.08 * step)
      } else {
        step <- step / 2
      }
      if (step < step_tolerance) {
        converged <- TRUE
        break
      }
    }
    list(u = u, omega = omega, evaluations = evaluations,
         iterations = iterations, converged = converged)
  })
}

.rv_direction_interval <- function(cache, x, y, bias_sign, omega_evaluator,
                                   n_starts = 16L, seed = 1234L,
                                   max_iterations = 60L, warm = list()) {
  screen <- .rv_direction_screen(cache, x, y, bias_sign, omega_evaluator,
                                 n_starts, seed, warm)
  low <- .rv_sphere_extreme(
    cache, x, y, bias_sign, omega_evaluator, screen$u_low,
    sense = "min", seed = seed + 100003L,
    max_iterations = max_iterations
  )
  high <- .rv_sphere_extreme(
    cache, x, y, bias_sign, omega_evaluator, screen$u_high,
    sense = "max", seed = seed + 200003L,
    max_iterations = max_iterations
  )
  list(
    u_low = low$u,
    u_high = high$u,
    omega_low = low$omega,
    omega_high = high$omega,
    evaluations = screen$evaluations + low$evaluations + high$evaluations,
    screen = screen,
    local_low = low,
    local_high = high
  )
}

.rv_slerp_basic <- function(u0, u1, t) {
  cosine <- min(max(.rv_dot(u0, u1), -1), 1)
  if (cosine > 1 - 1e-10) {
    value <- (1 - t) * u0 + t * u1
    return(value / .rv_norm(value))
  }
  angle <- acos(cosine)
  denominator <- sin(angle)
  value <- sin((1 - t) * angle) / denominator * u0 +
    sin(t * angle) / denominator * u1
  value / .rv_norm(value)
}

.rv_geodesic_path <- function(cache, u0, u1, seed = 1234L) {
  u0 <- .rv_normalize_direction(cache, u0)
  u1 <- .rv_normalize_direction(cache, u1)
  cosine <- .rv_dot(u0, u1)
  if (cosine > -1 + 1e-8) {
    return(function(t) .rv_slerp_basic(u0, u1, t))
  }
  candidates <- .rv_structured_directions(cache, 12L, seed)
  bridge <- NULL
  for (candidate in candidates) {
    tangent <- candidate - u0 * .rv_dot(u0, candidate)
    if (.rv_norm(tangent) > 1e-8) {
      bridge <- .rv_normalize_direction(cache, tangent)
      break
    }
  }
  if (is.null(bridge)) return(NULL)
  function(t) {
    if (t <= 0.5) {
      .rv_slerp_basic(u0, bridge, 2 * t)
    } else {
      .rv_slerp_basic(bridge, u1, 2 * t - 1)
    }
  }
}

.rv_omega_root <- function(cache, x, y, bias_sign, interval, omega_target,
                           omega_evaluator, seed = 1234L,
                           root_tolerance = 1e-9,
                           max_iterations = 100L) {
  if (abs(log(interval$omega_low / omega_target)) <= root_tolerance) {
    return(list(found = TRUE, u = interval$u_low,
                omega = interval$omega_low, error = 0,
                evaluations = 0L, iterations = 0L))
  }
  if (abs(log(interval$omega_high / omega_target)) <= root_tolerance) {
    return(list(found = TRUE, u = interval$u_high,
                omega = interval$omega_high, error = 0,
                evaluations = 0L, iterations = 0L))
  }
  if (!(interval$omega_low < omega_target &&
        omega_target < interval$omega_high)) {
    return(list(found = FALSE, reason = "target_not_bracketed",
                evaluations = 0L, iterations = 0L))
  }
  path <- .rv_geodesic_path(cache, interval$u_low, interval$u_high, seed)
  if (is.null(path)) {
    return(list(found = FALSE, reason = "antipodal_bridge_failure",
                evaluations = 0L, iterations = 0L))
  }
  value <- function(t) {
    log(omega_evaluator(x, y, bias_sign, path(t)) / omega_target)
  }
  lo <- 0
  hi <- 1
  f_lo <- value(lo)
  f_hi <- value(hi)
  evaluations <- 2L
  if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo > 0 || f_hi < 0) {
    return(list(found = FALSE, reason = "invalid_root_bracket",
                evaluations = evaluations, iterations = 0L))
  }
  iterations <- 0L
  for (iteration in seq_len(as.integer(max_iterations))) {
    iterations <- iteration
    mid <- (lo + hi) / 2
    f_mid <- value(mid)
    evaluations <- evaluations + 1L
    if (!is.finite(f_mid)) {
      return(list(found = FALSE, reason = "nonfinite_root_value",
                  evaluations = evaluations, iterations = iterations))
    }
    if (abs(f_mid) <= root_tolerance) {
      lo <- hi <- mid
      break
    }
    if (f_mid < 0) lo <- mid else hi <- mid
  }
  u <- .rv_normalize_direction(cache, path((lo + hi) / 2))
  omega <- omega_evaluator(x, y, bias_sign, u)
  error <- abs(log(omega / omega_target))
  list(found = is.finite(error) && error <= 5 * root_tolerance,
       u = u, omega = omega, error = error,
       evaluations = evaluations, iterations = iterations,
       reason = if (error <= 5 * root_tolerance) "ok" else "root_tolerance_failure")
}

.rv_validate_point_z <- function(cache, z, vcov, cluster, null,
                                 short_t, critical, omega_target,
                                 x_target, y_target,
                                 omega_tolerance = 5e-7,
                                 r2_tolerance = 5e-7,
                                 decision_tolerance = 5e-7,
                                 require_literal = FALSE) {
  attained <- tryCatch(
    evaluate_confounder_reference(cache, z, vcov = vcov,
                                  cluster = cluster, null = null),
    error = function(e) e
  )
  if (inherits(attained, "error")) {
    return(list(valid = FALSE, error = conditionMessage(attained)))
  }
  omega_error <- abs(log(attained$omega / omega_target))
  x_error <- abs(attained$r2_d - x_target)
  y_error <- abs(attained$r2_y - y_target)
  reject_short <- .rv_reject(short_t, critical)
  decision_violation <- .rv_decision_violation(
    reject_short, attained$t_robust, critical
  )
  literal <- .rv_decision_change(short_t, attained$t_robust, critical)
  valid <- omega_error <= omega_tolerance &&
    x_error <= r2_tolerance && y_error <= r2_tolerance &&
    decision_violation <= decision_tolerance &&
    (!require_literal || literal)
  list(
    valid = valid,
    attained = attained,
    z = as.numeric(z),
    radius = max(attained$r2_d, attained$r2_y),
    omega_error = omega_error,
    x_error = x_error,
    y_error = y_error,
    decision_violation = decision_violation,
    literal_change = literal,
    error = NULL
  )
}

.rv_fixed_xy_point_witness <- function(cache, x, y, bias_sign, omega_target,
                                       omega_evaluator, vcov, cluster, null,
                                       short_t, critical, n_starts, seed,
                                       max_iterations, root_tolerance,
                                       warm = list(), require_literal = FALSE) {
  interval <- .rv_direction_interval(
    cache, x, y, bias_sign, omega_evaluator,
    n_starts = n_starts, seed = seed,
    max_iterations = max_iterations, warm = warm
  )
  root <- .rv_omega_root(
    cache, x, y, bias_sign, interval, omega_target, omega_evaluator,
    seed = seed + 300007L, root_tolerance = root_tolerance
  )
  evaluations <- interval$evaluations + root$evaluations
  if (!isTRUE(root$found)) {
    return(list(valid = FALSE, interval = interval, root = root,
                evaluations = evaluations))
  }
  z <- single_parameterization(
    cache, sqrt(x), bias_sign * sqrt(y), root$u
  )
  validation <- .rv_validate_point_z(
    cache, z, vcov, cluster, null, short_t, critical, omega_target,
    x, y, require_literal = require_literal,
    omega_tolerance = max(5e-7, 10 * root_tolerance)
  )
  c(validation, list(
    interval = interval,
    root = root,
    u = root$u,
    bias_sign = bias_sign,
    x = x,
    y = y,
    evaluations = evaluations
  ))
}

.rv_protected_signs <- function(cache, reject_short, null = 0) {
  short <- .rv_short_quantities(cache, null)
  coefficient_sign <- if (short$t_short_classical < 0) -1 else 1
  expected <- if (reject_short) coefficient_sign else -coefficient_sign
  c(expected, -expected)
}

.rv_protected_certificate <- function(cache, protected, omega_target,
                                      omega_evaluator, vcov, cluster, null,
                                      short_t, critical, n_starts, seed,
                                      max_iterations, root_tolerance,
                                      require_literal = FALSE,
                                      diagonal_only = FALSE,
                                      warm = list()) {
  regime <- .rv_protected_regime(protected$branch)
  if (diagonal_only && regime != "equal_strength") {
    return(list(found = FALSE, attempts = list(), evaluations = 0L,
                reason = "protected_optimizer_is_not_diagonal"))
  }
  reject_short <- .rv_reject(short_t, critical)
  attempts <- list()
  evaluations <- 0L
  signs <- .rv_protected_signs(cache, reject_short, null)
  for (index in seq_along(signs)) {
    witness <- .rv_fixed_xy_point_witness(
      cache, protected$x, protected$y, signs[[index]], omega_target,
      omega_evaluator, vcov, cluster, null, short_t, critical,
      n_starts = n_starts, seed = seed + 100003L * index,
      max_iterations = max_iterations, root_tolerance = root_tolerance,
      warm = warm, require_literal = require_literal
    )
    attempts[[index]] <- witness
    evaluations <- evaluations + witness$evaluations
    if (isTRUE(witness$valid)) {
      return(list(found = TRUE, witness = witness, attempts = attempts,
                  evaluations = evaluations, regime = regime))
    }
  }
  list(found = FALSE, attempts = attempts, evaluations = evaluations,
       regime = regime, reason = "protected_target_not_certified")
}

.rv_branch_pair <- function(cache, branch, radius) {
  f <- .rv_short_quantities(cache)$f
  if (branch == "asymmetric_decrease") {
    c(radius, radius / (f^2 * (1 - radius)))
  } else if (branch == "asymmetric_increase") {
    c(radius / (f^2 + radius), radius)
  } else {
    c(radius, radius)
  }
}

.rv_candidate_pairs <- function(cache, radius, branch = NULL,
                                n_random = 8L, seed = 1234L) {
  radius <- min(max(radius, 0), 1 - 1e-8)
  pairs <- rbind(
    c(radius, radius),
    c(radius, radius / 2),
    c(radius / 2, radius),
    c(radius, radius / 10),
    c(radius / 10, radius)
  )
  if (!is.null(branch)) {
    branch_pair <- .rv_branch_pair(cache, branch, radius)
    if (all(is.finite(branch_pair)) && all(branch_pair >= 0) &&
        all(branch_pair < 1) && max(branch_pair) <= radius + 1e-10) {
      pairs <- rbind(branch_pair, pairs)
    }
  }
  random_pairs <- .rv_with_seed(seed, {
    out <- matrix(0, nrow = as.integer(n_random), ncol = 2L)
    for (i in seq_len(nrow(out))) {
      if (stats::runif(1) < 0.5) {
        out[i, ] <- c(radius, stats::runif(1, 0, radius))
      } else {
        out[i, ] <- c(stats::runif(1, 0, radius), radius)
      }
    }
    out
  })
  pairs <- unique(round(rbind(pairs, random_pairs), 14L))
  pairs[pairs[, 1L] < 1 & pairs[, 2L] < 1, , drop = FALSE]
}

.rv_point_candidate_at_radius <- function(cache, radius, protected,
                                          omega_target, omega_evaluator,
                                          vcov, cluster, null, short_t,
                                          critical, n_starts, seed,
                                          max_iterations, root_tolerance,
                                          n_random = 8L, warm = list()) {
  reject_short <- .rv_reject(short_t, critical)
  signs <- .rv_protected_signs(cache, reject_short, null)
  pairs <- .rv_candidate_pairs(cache, radius, protected$branch,
                               n_random, seed)
  evaluations <- 0L
  best <- NULL
  for (pair_index in seq_len(nrow(pairs))) {
    x <- pairs[pair_index, 1L]
    y <- pairs[pair_index, 2L]
    for (sign_index in seq_along(signs)) {
      bias_sign <- signs[[sign_index]]
      t_classical <- .rv_classical_ch_t_signed(cache, x, y, bias_sign, null)
      t_at_target <- t_classical / omega_target
      closure_ok <- if (reject_short) {
        abs(t_at_target) <= critical + 1e-8
      } else {
        abs(t_at_target) >= critical - 1e-8
      }
      if (!closure_ok) next
      witness <- .rv_fixed_xy_point_witness(
        cache, x, y, bias_sign, omega_target, omega_evaluator,
        vcov, cluster, null, short_t, critical,
        n_starts = n_starts,
        seed = seed + 1009L * pair_index + 9176L * sign_index,
        max_iterations = max_iterations,
        root_tolerance = root_tolerance,
        warm = warm,
        require_literal = FALSE
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

.rv_search_point_upper <- function(cache, protected, omega_target,
                                   omega_evaluator, vcov, cluster, null,
                                   short_t, critical, lower, n_starts,
                                   seed, max_iterations, root_tolerance,
                                   max_radius = 0.9, radial_points = 18L) {
  lower_search <- max(lower, 0)
  radii <- unique(c(
    lower_search,
    pmin(max_radius, lower_search + c(1e-7, 1e-5, 1e-3, 1e-2)),
    seq(lower_search, max_radius, length.out = radial_points)
  ))
  radii <- sort(radii[radii < 1])
  previous <- lower_search
  best <- NULL
  evaluations <- 0L
  for (index in seq_along(radii)) {
    candidate <- .rv_point_candidate_at_radius(
      cache, radii[[index]], protected, omega_target, omega_evaluator,
      vcov, cluster, null, short_t, critical, n_starts,
      seed + 100003L * index, max_iterations, root_tolerance,
      n_random = max(6L, n_starts %/% 2L),
      warm = if (is.null(best)) list() else list(best$u)
    )
    evaluations <- evaluations + candidate$evaluations
    if (candidate$found) {
      best <- candidate$witness
      lo <- previous
      hi <- radii[[index]]
      for (step in seq_len(10L)) {
        mid <- (lo + hi) / 2
        refined <- .rv_point_candidate_at_radius(
          cache, mid, protected, omega_target, omega_evaluator,
          vcov, cluster, null, short_t, critical, n_starts,
          seed + 700001L + step, max_iterations, root_tolerance,
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

.rv_point_result <- function(base, certificate, fallback, started,
                             exact_tolerance) {
  if (isTRUE(certificate$found)) {
    w <- certificate$witness
    return(structure(c(base, list(
      rv = base$RV_prot,
      RV_att = base$RV_prot,
      RV_lower = base$RV_prot,
      RV_upper = base$RV_prot,
      RV_gap = 0,
      exact = TRUE,
      status = "EXACT_POINT_RV",
      substatus = if (certificate$regime == "equal_strength")
        "PROTECTED_DIAGONAL" else "PROTECTED_INTERIOR",
      x_star = w$x,
      y_star = w$y,
      sign_star = w$bias_sign,
      best_witness_z = w$z,
      best_witness_u = w$u,
      best_witness_omega = w$attained$omega,
      best_witness_t = w$attained$t_robust,
      direct_refit_valid = TRUE,
      literal_decision_change = w$literal_change,
      objective_evals = certificate$evaluations,
      runtime = .rv_timer() - started
    )), class = "robust_rv_result"))
  }
  if (!isTRUE(fallback$found)) {
    return(structure(c(base, list(
      rv = NA_real_, RV_att = NA_real_, RV_lower = base$RV_prot,
      RV_upper = NA_real_, RV_gap = NA_real_, exact = FALSE,
      status = "LOWER_BOUND_ONLY", substatus = "POINT_WITNESS_NOT_FOUND",
      best_witness_z = NULL, direct_refit_valid = FALSE,
      objective_evals = certificate$evaluations + fallback$evaluations,
      runtime = .rv_timer() - started
    )), class = "robust_rv_result"))
  }
  w <- fallback$witness
  gap <- w$radius - base$RV_prot
  closed <- gap <= exact_tolerance
  structure(c(base, list(
    rv = if (closed) base$RV_prot else NA_real_,
    RV_att = if (closed) base$RV_prot else NA_real_,
    RV_lower = base$RV_prot,
    RV_upper = w$radius,
    RV_gap = gap,
    exact = closed,
    status = if (closed) "POINT_BOUNDS_CLOSED" else "CERTIFIED_POINT_BOUNDS",
    substatus = "POINT_WITNESS_FOUND",
    x_star = w$x,
    y_star = w$y,
    sign_star = w$bias_sign,
    best_witness_z = w$z,
    best_witness_u = w$u,
    best_witness_omega = w$attained$omega,
    best_witness_t = w$attained$t_robust,
    direct_refit_valid = TRUE,
    literal_decision_change = w$literal_change,
    objective_evals = certificate$evaluations + fallback$evaluations,
    runtime = .rv_timer() - started
  )), class = "robust_rv_result")
}

#' Point-conditioned attainable robustness value
#'
#' @param cache An `rv_cache`.
#' @param omega_target Target sandwich-to-classical standard-error ratio.
#' @param method `"protected_first"`, `"diagonal_first"`, or `"full"`.
#' @param vcov HC0, HC1, CR0, or CR1.
#' @param cluster Optional one-way cluster labels.
#' @param alpha Two-sided test level.
#' @param null Null coefficient value.
#' @param critical Optional fixed critical value.
#' @param n_starts Number of structured/random sphere starts.
#' @param seed Random seed used locally without changing `.Random.seed`.
#' @param max_iterations Maximum local sphere iterations.
#' @param root_tolerance Log-omega root tolerance.
#' @param exact_tolerance Gap used to label a numerically closed interval.
#' @param max_radius Largest partial-R-squared radius searched by the fallback.
#' @return A `robust_rv_result` containing a theorem certificate or certified
#'   numerical bounds.
#' @export
conditional_omega_rv <- function(cache, omega_target,
                                 method = c("protected_first", "diagonal_first", "full"),
                                 vcov = "HC1", cluster = NULL,
                                 alpha = 0.05, null = 0, critical = NULL,
                                 n_starts = 16L, seed = 1234L,
                                 max_iterations = 60L,
                                 root_tolerance = 1e-9,
                                 exact_tolerance = 2e-4,
                                 max_radius = 0.9) {
  started <- .rv_timer()
  method <- match.arg(method)
  vcov <- match.arg(vcov, c("HC0", "HC1", "CR0", "CR1"))
  cluster <- .rv_resolve_cluster(cluster, cache)
  if (!.rv_is_scalar_number(omega_target) || omega_target <= 0) {
    stop("omega_target must be positive and finite", call. = FALSE)
  }
  critical <- .rv_critical_value(cache, vcov, 1L, alpha, cluster, critical)
  short <- .rv_short_evaluation(cache, vcov, cluster, null)
  reject_short <- .rv_reject(short$t_robust, critical)
  RV_eq <- .rv_conditional_equal(cache, omega_target, critical,
                                 reject_short, null)
  protected <- .rv_conditional_protected(cache, omega_target, critical,
                                         reject_short, null)
  omega_evaluator <- conditional_omega_evaluator(cache, vcov, cluster)
  base <- list(
    calibration_mode = "point",
    method = method,
    vcov = vcov,
    omega_target = omega_target,
    critical_value = critical,
    baseline_decision = if (reject_short) "reject" else "do_not_reject",
    short_t_robust = short$t_robust,
    RV_eq = RV_eq,
    RV_prot = protected$rv,
    protected_branch = protected$branch,
    protected_regime = .rv_protected_regime(protected$branch),
    x_prot = protected$x,
    y_prot = protected$y
  )
  if (method == "full") {
    certificate <- list(found = FALSE, evaluations = 0L,
                        reason = "full_search_requested")
  } else {
    certificate <- .rv_protected_certificate(
      cache, protected, omega_target, omega_evaluator, vcov, cluster,
      null, short$t_robust, critical,
      n_starts = n_starts, seed = seed,
      max_iterations = max_iterations,
      root_tolerance = root_tolerance,
      diagonal_only = method == "diagonal_first"
    )
  }
  if (isTRUE(certificate$found)) {
    return(.rv_point_result(base, certificate,
                            list(found = FALSE, evaluations = 0L),
                            started, exact_tolerance))
  }
  fallback <- .rv_search_point_upper(
    cache, protected, omega_target, omega_evaluator, vcov, cluster,
    null, short$t_robust, critical, protected$rv,
    n_starts = if (method == "full") max(24L, 2L * n_starts) else n_starts,
    seed = seed + 900001L, max_iterations = max_iterations,
    root_tolerance = root_tolerance, max_radius = max_radius,
    radial_points = if (method == "full") 24L else 18L
  )
  .rv_point_result(base, certificate, fallback, started, exact_tolerance)
}

#' @rdname conditional_omega_rv
#' @export
rv_att_point <- conditional_omega_rv
