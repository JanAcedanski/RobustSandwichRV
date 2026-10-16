# Loss-of-significance estimands in the current manuscript. The older point
# and symmetric-band interfaces remain available under their original names.

.rv_loss_context <- function(model, treatment, vcov, cluster, alpha, critical,
                             null) {
  if (!.rv_is_scalar_number(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must lie in (0,1)", call. = FALSE)
  cache <- if (inherits(model, "rv_cache")) model else
    cache_from_lm(model, treatment)
  if (cache$rank_x != ncol(cache$X))
    stop("the maintained design must have full column rank", call. = FALSE)
  if (.rv_norm(cache$e_short) <= 100 * .Machine$double.eps * max(.rv_norm(cache$y), 1))
    stop("the observed model has zero residual variation", call. = FALSE)
  if (.rv_norm(.rv_residualize(cache$Qx, rep(1, cache$n))) > 1e-7 * sqrt(cache$n))
    stop("the maintained controls must include an intercept", call. = FALSE)
  if (!.rv_is_scalar_number(null)) stop("null must be finite", call. = FALSE)
  cluster <- .rv_resolve_cluster(cluster, cache)
  rule <- .rv_resolve_rule(vcov, cache$n, cluster)
  critical <- .rv_critical_value(cache, rule, 1L, alpha, NULL, critical)
  short <- .rv_short_evaluation(cache, rule, NULL, null)
  list(cache = cache, vcov = rule$name, rule = rule, cluster = NULL, critical = critical,
       short = short, null = null, alpha = alpha,
       treatment = treatment %||% cache$treatment)
}

.rv_loss_chi <- function(ctx, df) {
  df * .rv_rule_correction(ctx$rule, df)
}

# A continuous upper envelope valid for EVERY scalar index with x,y <= r.
# Expand qhat*ehat, use the triangle inequality in the score-operator norm,
# ||d*u|| <= ||d||_inf, ||e*u|| <= ||e||_inf, ||u*u|| <= 1, and
# ||qhat*ehat||^2 <= 1/2 for orthogonal unit residuals. Rule constructors
# supply valid spectral bounds. Indefinite rules use the signed expansion below.
.rv_loss_envelope <- function(ctx) {
  d <- ctx$cache$d / .rv_norm(ctx$cache$d)
  e <- ctx$cache$e_short / .rv_norm(ctx$cache$e_short)
  chi <- .rv_loss_chi(ctx, ctx$cache$df_short - 1L)
  lambda <- ctx$rule$norm_bound
  a0 <- .rv_omega_origin(ctx$cache, ctx$rule)
  multiplier <- sqrt(chi * lambda)
  global <- sqrt(chi * ctx$rule$positive_bound / 2)
  c1 <- max(abs(d)) + max(abs(e))
  c2 <- sqrt(sum(e^4)) + 1
  if (ctx$rule$psd) return(list(a0 = a0, global = global,
    bound = function(r) min(global, a0 + multiplier * (c1 * sqrt(r) + c2 * r))))
  # For indefinite operators a seminorm triangle inequality is invalid.
  # Write scores=a*d*e+delta with 0<=a<=1. Expansion of the signed form gives
  # a0^2 + chi*(2*||A(d*e)||*||delta|| + ||A||*||delta||^2).
  Ad <- .rv_norm(ctx$rule$apply(d * e))
  list(a0 = a0, global = global, bound = function(r) {
    delta <- c1 * sqrt(r) + c2 * r
    min(global, sqrt(a0^2 + chi * (2 * Ad * delta + lambda * delta^2)))
  })
}

.rv_loss_lower <- function(ctx, cap) {
  lo <- 0; hi <- 1
  for (i in seq_len(70L)) {
    mid <- (lo + hi) / 2
    p <- .rv_conditional_protected(ctx$cache, cap(mid), ctx$critical, TRUE, ctx$null)
    if (p$rv > mid) lo <- mid else hi <- mid
  }
  max(0, lo - 1e-12)
}

.rv_loss_controls <- function(n_starts, max_iterations, seed, tolerance,
                              max_radius, search_steps, run_search) {
  integer_arg <- function(x, name, minimum, maximum = 1e8) {
    if (!.rv_is_scalar_number(x) || x != floor(x) || x < minimum || x > maximum)
      stop(name, " is outside its supported integer range", call. = FALSE)
    as.integer(x)
  }
  if (!.rv_is_scalar_number(tolerance) || tolerance <= 0 || tolerance > 1e-3)
    stop("tolerance must be positive and at most 1e-3", call. = FALSE)
  if (!.rv_is_scalar_number(max_radius) || max_radius <= 0 || max_radius >= 1)
    stop("max_radius must lie in (0,1)", call. = FALSE)
  if (!is.logical(run_search) || length(run_search) != 1L || is.na(run_search))
    stop("run_search must be TRUE or FALSE", call. = FALSE)
  list(n_starts = integer_arg(n_starts, "n_starts", 2),
       max_iterations = integer_arg(max_iterations, "max_iterations", 1),
       seed = integer_arg(seed, "seed", 0), tolerance = tolerance,
       max_radius = max_radius,
       search_steps = integer_arg(search_steps, "search_steps", 1, 60),
       run_search = run_search)
}

# Gradient with respect to a direction in the residual sphere. Projection and
# normalization are differentiated as well; optim's failed searches never
# contribute a lower bound. Only independently refitted vectors give uppers.
.rv_loss_direction_objective <- function(ctx, x, y, s, v) {
  cache <- ctx$cache
  projected <- .rv_project_complement(cache, v)
  nv <- .rv_norm(projected)
  if (nv < 1e-14) stop("degenerate direction", call. = FALSE)
  u <- projected / nv
  d <- cache$d / .rv_norm(cache$d); e <- cache$e_short / .rv_norm(cache$e_short)
  rd <- sqrt(x); ry <- s * sqrt(y); k <- sqrt(1 - y)
  q <- sqrt(1 - x) * d - rd * ry * e - rd * k * u
  residual <- k * e - ry * u
  score <- q * residual
  applied <- ctx$rule$apply(score)
  meat <- sum(score * applied)
  if (!is.finite(meat)) stop("nonfinite score quadratic form", call. = FALSE)
  chi <- .rv_loss_chi(ctx, cache$df_short - 1L)
  gradient <- 2 * chi * (-rd * k * residual - ry * q) * applied
  gradient <- .rv_project_complement(cache, gradient)
  gradient <- (gradient - u * sum(u * gradient)) / nv
  list(value = chi * meat, gradient = gradient, u = u,
       omega = if (meat > 0) sqrt(chi * meat) else NA_real_)
}

.rv_loss_direction_interval <- function(ctx, x, y, s, control, warm = list()) {
  starts <- .rv_structured_directions(ctx$cache, control$n_starts, control$seed, warm)
  values <- vapply(starts, function(u)
    .rv_loss_direction_objective(ctx, x, y, s, u)$value, numeric(1))
  low <- list(u = starts[[which.min(values)]], value = min(values))
  high <- list(u = starts[[which.max(values)]], value = max(values))
  evaluations <- length(starts)
  failures <- 0L
  # Use several screened starts in each direction, including supplied witnesses.
  for (sense in c(-1, 1)) {
    order <- order(sense * values)
    for (i in head(order, max(2L, ceiling(control$n_starts / 4)))) {
      last_v <- NULL; last_state <- NULL
      state <- function(v) {
        if (is.null(last_v) || !identical(v, last_v)) {
          last_state <<- .rv_loss_direction_objective(ctx, x, y, s, v)
          last_v <<- v
          evaluations <<- evaluations + 1L
        }
        last_state
      }
      fit <- tryCatch(stats::optim(starts[[i]],
        function(v) sense * state(v)$value,
        function(v) sense * state(v)$gradient, method = "L-BFGS-B",
        control = list(maxit = control$max_iterations, factr = 1e5, pgtol = 1e-7)),
        error = function(e) NULL)
      if (is.null(fit)) { failures <- failures + 1L; next }
      candidate <- state(fit$par)
      if (candidate$value < low$value) low <- candidate
      if (candidate$value > high$value) high <- candidate
    }
  }
  list(u_low = low$u, u_high = high$u, quadratic_low = low$value,
       quadratic_high = high$value, omega_low = sqrt(max(0, low$value)),
       omega_high = sqrt(max(0, high$value)), evaluations = evaluations, failures = failures)
}

.rv_loss_check <- function(ctx, z, cap, tolerance) {
  attained <- tryCatch(evaluate_confounder_reference(ctx$cache, z,
    vcov = ctx$rule, null = ctx$null), error = function(e) NULL)
  if (is.null(attained)) return(list(valid = FALSE))
  radius <- max(attained$r2_d, attained$r2_y)
  actual_cap <- cap(attained$r2_d, attained$r2_y)
  violation <- max(log(attained$omega / actual_cap), 0)
  valid <- is.finite(attained$variance_robust) && attained$variance_robust > 0 &&
    attained$r2_d < 1 && attained$r2_y < 1 &&
    violation <= tolerance && abs(attained$t_robust) <= ctx$critical + tolerance
  list(valid = valid, z = as.numeric(z), radius = radius, attained = attained,
       cap = actual_cap, cap_violation = violation,
       literal_change = abs(attained$t_robust) < ctx$critical)
}

# A positive target can be bracketed across negative raw two-way variances.
# The signed quadratic form is continuous along the path; log(omega) is not.
.rv_loss_omega_root <- function(ctx, x, y, s, interval, target, seed) {
  target2 <- target^2
  tolerance <- 1e-11 * target2
  for (end in c("low", "high")) {
    if (abs(interval[[paste0("quadratic_", end)]] - target2) <= tolerance)
      return(list(found = TRUE, u = interval[[paste0("u_", end)]], evaluations = 0L))
  }
  if (interval$quadratic_low > target2 || interval$quadratic_high < target2)
    return(list(found = FALSE, evaluations = 0L))
  path <- .rv_geodesic_path(ctx$cache, interval$u_low, interval$u_high, seed)
  if (is.null(path)) return(list(found = FALSE, evaluations = 0L))
  evaluate <- function(u) {
    residuals <- alignment_directions(ctx$cache, x, y, s, u)
    score <- residuals$q_hat * residuals$r_hat
    .rv_loss_chi(ctx, ctx$cache$df_short - 1L) * sum(score * ctx$rule$apply(score))
  }
  lo <- 0; hi <- 1; evaluations <- 0L
  for (i in seq_len(100L)) {
    mid <- (lo + hi) / 2; u <- path(mid)
    value <- evaluate(u); evaluations <- evaluations + 1L
    if (!is.finite(value)) return(list(found = FALSE, evaluations = evaluations))
    if (abs(value - target2) <= tolerance) break
    if (value < target2) lo <- mid else hi <- mid
  }
  list(found = value > 0 && abs(value - target2) <= 10 * tolerance,
       u = u, evaluations = evaluations)
}

.rv_loss_xy_witness <- function(ctx, x, y, s, cap, control, endpoint = FALSE,
                                warm = list()) {
  interval <- .rv_loss_direction_interval(ctx, x, y, s, control, warm)
  if (interval$quadratic_high <= 0)
    return(list(valid = FALSE, evaluations = interval$evaluations))
  threshold <- abs(.rv_classical_ch_t_signed(ctx$cache, x, y, s, ctx$null)) / ctx$critical
  upper <- min(cap(x, y), interval$omega_high)
  lower <- max(threshold, interval$omega_low)
  if (!endpoint) {
    if (lower > upper * (1 + 1e-12))
      return(list(valid = FALSE, evaluations = interval$evaluations))
    directions <- list(interval$u_high, interval$u_low)
    omegas <- c(interval$omega_high, interval$omega_low)
    for (i in seq_along(directions)) {
      if (omegas[i] < threshold * (1 - 1e-12) ||
          omegas[i] > cap(x, y) * (1 + 1e-12)) next
      u <- directions[[i]]
      z <- single_parameterization(ctx$cache, sqrt(x), s * sqrt(y), u)
      check <- .rv_loss_check(ctx, z, cap, control$tolerance)
      if (isTRUE(check$valid))
        return(c(check, list(u = u, evaluations = interval$evaluations)))
    }
  }
  if (endpoint) {
    target <- cap(x, y)
    if (target < interval$omega_low * (1 - 1e-10) ||
        target > interval$omega_high * (1 + 1e-10))
      return(list(valid = FALSE, evaluations = interval$evaluations))
  } else {
    if (lower > upper * (1 + 1e-12))
      return(list(valid = FALSE, evaluations = interval$evaluations))
    target <- if (lower == 0) upper / 2 else sqrt(lower * upper)
  }
  root <- .rv_loss_omega_root(ctx, x, y, s, interval, target, control$seed)
  if (!isTRUE(root$found))
    return(list(valid = FALSE, evaluations = interval$evaluations + root$evaluations))
  z <- single_parameterization(ctx$cache, sqrt(x), s * sqrt(y), root$u)
  check <- .rv_loss_check(ctx, z, cap, control$tolerance)
  c(check, list(u = root$u, evaluations = interval$evaluations + root$evaluations))
}

.rv_loss_result <- function(ctx, mode, lower, best, elapsed, evaluations,
                            control, envelope, cap, status = NULL,
                            benchmark = NULL, profile = NULL, search_error = "") {
  upper <- if (is.null(best)) NA_real_ else best$radius
  if (is.finite(upper) && lower > upper + control$tolerance)
    stop("analytical lower bound exceeds independently verified upper bound", call. = FALSE)
  if (is.finite(upper)) lower <- min(lower, upper)
  closed <- is.finite(upper) && upper - lower <= control$tolerance
  if (is.null(status)) status <- if (closed) "POINT_CERTIFICATE" else
    if (is.finite(upper)) "CERTIFIED_BOUNDS" else "LOWER_BOUND_ONLY"
  benchmark_status <- if (is.null(benchmark)) NA_character_ else
    if (!benchmark$complete || (mode == "profile" && !benchmark$profile_available)) "unavailable" else
    if (mode == "profile") {
      if (benchmark$L == 0 && benchmark$beta == 0) "zero" else "positive"
    } else if (benchmark$delta_positive == 0) "zero" else "positive"
  structure(list(calibration_mode = mode, measure = switch(mode,
    sand = "RV_sand", increment = "RV_cal_delta", profile = "RV_cal_profile", cap = "RV_cal_cap"),
    rv = if (closed) lower else NA_real_, RV_lower = lower, RV_upper = upper,
    RV_gap = if (is.finite(upper)) upper - lower else NA_real_, exact = closed,
    theorem_exact = FALSE, numerically_closed = closed,
    certificate_type = if (closed) "verified_interval_within_tolerance" else
      if (is.finite(upper)) "verified_interval" else "analytical_lower_bound_only",
    vcov = ctx$vcov, covariance_rule = ctx$rule,
    critical_value = ctx$critical, alpha = ctx$alpha, null = ctx$null,
    baseline_decision = if (abs(ctx$short$t_robust) > ctx$critical) "reject" else "do_not_reject",
    short_t_robust = ctx$short$t_robust, omega_origin = envelope$a0,
    omega_lower = 0, omega_upper = cap(1, 1), omega_adverse = cap(1, 1),
    status = status, substatus = mode, runtime = elapsed, objective_evals = evaluations,
    best_witness_z = if (is.null(best)) NULL else best$z,
    witness_refit = if (is.null(best)) NULL else best$attained,
    best_witness_x = if (is.null(best)) NA_real_ else best$attained$r2_d,
    best_witness_y = if (is.null(best)) NA_real_ else best$attained$r2_y,
    best_witness_omega = if (is.null(best)) NA_real_ else best$attained$omega,
    best_witness_t = if (is.null(best)) NA_real_ else best$attained$t_robust,
    best_witness_cap = if (is.null(best)) NA_real_ else best$cap,
    direct_refit_valid = !is.null(best), witness_literal_decision_change =
      !is.null(best) && best$literal_change,
    treatment = ctx$treatment, benchmark = benchmark, benchmark_status = benchmark_status,
    profile = profile,
    control = control, search_error = search_error,
    search_max_radius = control$max_radius,
    model_call = if (is.null(ctx$cache$model)) NULL else ctx$cache$model$call),
    class = "robust_rv_result")
}

.rv_loss_solve <- function(ctx, mode, control, M = NULL, benchmark = NULL,
                           k = 1, supplied_witnesses = list()) {
  started <- .rv_timer(); evaluations <- 0L
  envelope <- .rv_loss_envelope(ctx)
  profile <- if (mode == "profile") list(a0 = envelope$a0, L = benchmark$L,
                                         beta = benchmark$beta, k = k) else NULL
  cap <- switch(mode, sand = function(x, y) envelope$global,
    increment = function(x, y) benchmark$M_delta,
    cap = function(x, y) M,
    profile = function(x, y) profile$a0 * exp(profile$k *
      min(profile$L, profile$beta * (sqrt(x) + sqrt(y)))))
  box <- function(r) min(envelope$bound(r), cap(r, r))
  lower <- .rv_loss_lower(ctx, box)
  if (abs(ctx$short$t_robust) <= ctx$critical) {
    out <- .rv_loss_result(ctx, mode, 0, NULL, .rv_timer() - started,
      0L, control, envelope, cap, status = "ALREADY_NONREJECTING",
      benchmark = benchmark, profile = profile)
    out$rv <- out$RV_upper <- out$RV_gap <- 0
    out$exact <- out$numerically_closed <- TRUE
    out$certificate_type <- "observed_nonrejection_convention"
    return(out)
  }
  best <- NULL; warm <- list(); search_error <- character()
  retain <- function(check) {
    if (isTRUE(check$valid) && (is.null(best) || check$radius < best$radius)) {
      best <<- check
      if (!is.null(check$u)) warm <<- c(list(check$u), head(warm, 3L))
    }
  }
  for (candidate in supplied_witnesses) {
    z <- if (is.numeric(candidate)) candidate else candidate$best_witness_z %||% candidate$z
    if (!is.null(z)) retain(.rv_loss_check(ctx, z, cap, control$tolerance))
  }
  signs <- .rv_protected_signs(ctx$cache, TRUE, ctx$null)
  attempt <- function(x, y, s, endpoint = FALSE) {
    if (x < 0 || y < 0 || x >= 1 || y >= 1) return(invisible(NULL))
    check <- tryCatch(.rv_loss_xy_witness(ctx, x, y, s, cap, control,
      endpoint = endpoint, warm = warm), error = function(e) {
        search_error <<- unique(c(search_error, conditionMessage(e)))
        list(valid = FALSE, evaluations = 0L)
      })
    evaluations <<- evaluations + check$evaluations
    retain(check)
    invisible(check)
  }
  # A zero-strength index has a fixed ratio. This also handles mechanical
  # HC0/CR0 dimension corrections without substituting the observed ratio.
  if (envelope$a0 <= cap(0, 0) * (1 + control$tolerance)) {
    u <- .rv_structured_directions(ctx$cache, 2L, control$seed)[[1L]]
    z <- single_parameterization(ctx$cache, 0, 0, u)
    retain(.rv_loss_check(ctx, z, cap, control$tolerance))
  }
  if (mode != "sand") {
    endpoint_cap <- cap(lower, lower)
    protected <- .rv_conditional_protected(ctx$cache, endpoint_cap,
      ctx$critical, TRUE, ctx$null)
    for (s in signs) {
      attempt(protected$x, protected$y, s, endpoint = TRUE)
      if (!is.null(best) && best$radius - lower <= control$tolerance) break
    }
  }
  if (control$run_search && (is.null(best) || best$radius - lower > control$tolerance)) {
    # Always keep a directly checked coefficient-nullifying index if its
    # actual SE ratio satisfies this calibration. It gives a useful initial U.
    point <- point_estimate_rv(ctx$cache, ctx$null)
    if (point < 1) {
      for (u in head(.rv_structured_directions(ctx$cache, control$n_starts,
                                              control$seed), control$n_starts)) {
        z <- single_parameterization(ctx$cache, sqrt(point), signs[[1L]] * sqrt(point), u)
        retain(.rv_loss_check(ctx, z, cap, control$tolerance))
      }
    }
    lo_search <- lower
    hi_search <- if (is.null(best)) control$max_radius else min(best$radius, control$max_radius)
    # A failed query only changes where we search. It never changes `lower`.
    for (step in seq_len(control$search_steps)) {
      if (hi_search < lo_search || hi_search - lo_search <= control$tolerance) break
      radius <- if (is.null(best) && step == 1L) hi_search else (lo_search + hi_search) / 2
      f <- .rv_short_quantities(ctx$cache, ctx$null)$f
      asymmetric_y <- min(radius, radius / (f^2 * (1 - radius)))
      pairs <- unique(rbind(c(radius, radius), c(radius, 0), c(0, radius),
        c(radius, asymmetric_y), c(radius, radius / 2), c(radius / 2, radius)))
      old_best <- if (is.null(best)) Inf else best$radius
      for (j in seq_len(nrow(pairs))) {
        for (s in signs) {
          attempt(pairs[j, 1L], pairs[j, 2L], s)
          if (!is.null(best) && best$radius <= radius + control$tolerance) break
        }
        if (!is.null(best) && best$radius <= radius + control$tolerance) break
      }
      if (!is.null(best) && best$radius <= radius + control$tolerance) {
        hi_search <- min(radius, best$radius)
      } else lo_search <- radius
      if (!is.null(best) && best$radius - lower <= control$tolerance) break
      if (is.infinite(old_best) && is.null(best) && step == 1L) break
    }
  }
  .rv_loss_result(ctx, mode, lower, best, .rv_timer() - started, evaluations,
    control, envelope, cap, benchmark = benchmark, profile = profile,
    search_error = paste(search_error, collapse = "; "))
}

#' Observed-control benchmarks for increment and profile calibration
#'
#' Keeps the estimation rows and covariance rule fixed when deleting each
#' prescribed control unit. Empty or incomplete pools are unavailable. Negative
#' increments are retained. Profile calibration requires scalar control units.
#'
#' @param model An unweighted, full-rank lm model including an intercept.
#' @param treatment Target model-matrix column name.
#' @param benchmarks Character vector of separate control terms, or a named
#'   list whose elements specify whole blocks. NULL selects scalar terms.
#' @param vcov A built-in estimator name or a fixed HAC/two-way rule.
#' @param cluster One-way labels, two-column labels, or a cluster formula.
#' @return An rv_calibration with the complete deletion table and both benchmarks.
#' @export
rv_benchmarks <- function(model, treatment, benchmarks = NULL, vcov = "HC1",
                           cluster = NULL) {
  ctx <- .rv_loss_context(model, treatment, vcov, cluster, 0.05, NULL, 0)
  extracted <- .rv_extract_lm(model, treatment)
  if (is.null(benchmarks)) {
    candidates <- extracted$term_labels
    benchmarks <- Filter(function(term) {
      columns <- which(extracted$assign == match(term, extracted$term_labels))
      length(columns) == 1L && !(extracted$treatment_index %in% columns) &&
        diff(range(extracted$model_matrix[, columns])) > 0
    }, candidates)
  }
  if (is.character(benchmarks)) units <- as.list(benchmarks) else units <- benchmarks
  if (!is.list(units)) stop("benchmarks must be a character vector or list", call. = FALSE)
  labels <- names(units)
  if (is.null(labels)) labels <- vapply(units, paste, collapse = " + ", FUN.VALUE = character(1))
  if (any(!nzchar(labels)) || anyDuplicated(labels))
    stop("benchmark units need unique, nonempty names", call. = FALSE)
  full <- ctx$short; records <- list(); seen <- integer()
  for (i in seq_along(units)) {
    selected <- .rv_benchmark_column_indices(extracted, units[[i]])
    if (any(extracted$assign[selected] == 0L))
      stop("the intercept cannot be a benchmark", call. = FALSE)
    if (length(intersect(seen, selected)))
      stop("benchmark units must not overlap", call. = FALSE)
    seen <- union(seen, selected)
    controls <- setdiff(seq_len(ncol(extracted$model_matrix)), extracted$treatment_index)
    reduced_X <- extracted$X[, -match(selected, controls), drop = FALSE]
    failure <- ""
    row <- tryCatch({
      reduced <- prepare_model(extracted$y, extracted$d, reduced_X)
      if (.rv_norm(.rv_residualize(reduced$Qx, rep(1, reduced$n))) > 1e-7 * sqrt(reduced$n))
        stop("benchmark deletion removes the intercept")
      minus <- .rv_short_evaluation(reduced, ctx$rule, NULL, 0)
      if (!is.finite(minus$omega) || minus$omega <= 100 * .Machine$double.eps)
        stop("reduced target score variance is numerically zero")
      rank <- ctx$cache$rank_x - reduced$rank_x
      if (rank != length(selected)) stop("benchmark is not identified")
      x <- 1 - sum(ctx$cache$d^2) / sum(reduced$d^2)
      y <- 1 - sum(ctx$cache$e_short^2) / sum(reduced$e_short^2)
      x <- max(0, min(1, x)); y <- max(0, min(1, y))
      common <- minus$omega * sqrt(.rv_loss_chi(ctx, ctx$cache$df_short) /
                                    .rv_loss_chi(ctx, reduced$df_short))
      ell <- log(full$omega / common)
      h <- sqrt(x) + sqrt(y)
      slope_valid <- h > 1e-12 || abs(ell) <= 1e-10
      list(rank = rank, omega_minus = minus$omega, omega_minus_common_chi = common,
           delta = full$omega - minus$omega, ell = ell, x = x, y = y,
           slope = if (h > 1e-12) max(ell, 0) / h else if (slope_valid) 0 else NA_real_,
           slope_valid = slope_valid)
    }, error = function(e) { failure <<- conditionMessage(e); NULL })
    records[[i]] <- data.frame(benchmark = labels[[i]],
      rank = if (is.null(row)) length(selected) else row$rank,
      omega_full = full$omega, omega_minus = if (is.null(row)) NA_real_ else row$omega_minus,
      omega_minus_common_chi = if (is.null(row)) NA_real_ else row$omega_minus_common_chi,
      delta_omega = if (is.null(row)) NA_real_ else row$delta,
      log_increment = if (is.null(row)) NA_real_ else row$ell,
      r2_d = if (is.null(row)) NA_real_ else row$x,
      r2_y = if (is.null(row)) NA_real_ else row$y,
      slope = if (is.null(row)) NA_real_ else row$slope,
      status = if (is.null(row)) "INVALID_REFIT" else "AVAILABLE",
      profile_valid = !is.null(row) && row$rank == 1L && row$slope_valid,
      error = failure, stringsAsFactors = FALSE)
  }
  table <- if (length(records)) do.call(rbind, records) else data.frame()
  complete <- length(records) > 0L && all(table$status == "AVAILABLE")
  scalar <- complete && all(table$rank == 1L)
  mixed <- complete && any(table$rank == 1L) && any(table$rank > 1L)
  if (mixed) stop("use separate scalar and block benchmark pools", call. = FALSE)
  delta <- if (complete) max(0, table$delta_omega) else NA_real_
  profile_ok <- scalar && all(table$profile_valid)
  L <- if (profile_ok) max(0, table$log_increment) else NA_real_
  beta <- if (profile_ok) max(table$slope) else NA_real_
  status <- if (!length(records)) "EMPTY_POOL" else if (!complete) "INCOMPLETE_POOL" else
    if (delta == 0) "COMPLETE_ZERO_INCREMENT" else "COMPLETE_POSITIVE_INCREMENT"
  structure(list(records = table, status = status, complete = complete,
    scalar = scalar, profile_available = profile_ok,
    omega_full = full$omega, M_delta = if (complete) full$omega + delta else NA_real_,
    delta_positive = delta, L = L, beta = beta,
    binding_delta = if (complete) table$benchmark[table$delta_omega == delta] else character(),
    binding_level = if (profile_ok) table$benchmark[table$log_increment == L] else character(),
    binding_slope = if (profile_ok) table$benchmark[table$slope == beta] else character(),
    treatment = treatment, vcov = ctx$vcov, covariance_rule = ctx$rule), class = "rv_calibration")
}

#' Sandwich robustness value for loss of significance
#'
#' Searches actual added scalar indices on the Cinelli-Hazlett strength scale.
#' An analytical continuous envelope supplies the lower bound; every upper
#' bound is independently refitted. Unresolved minima remain intervals.
#'
#' @param model An unweighted lm model, or an rv_cache.
#' @param treatment Target model-matrix column name (for lm).
#' @param vcov HC0/HC1, CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw,
#'   TWCR1_raw, or a fixed covariance rule.
#' @param cluster One-way labels, two-column labels, or a cluster formula.
#' @param alpha Two-sided significance level.
#' @param critical Optional positive cutoff, held fixed in the search.
#' @param null Null coefficient value.
#' @param n_starts Number of screened residual directions.
#' @param max_iterations Maximum iterations per local direction optimization.
#' @param seed Reproducible local seed; the caller's random state is preserved.
#' @param tolerance Numerical refit and interval-closing tolerance.
#' @param max_radius Largest searched strength, strictly below one.
#' @param search_steps Number of radial search queries.
#' @param run_search Whether to run the fallback search after endpoint checks.
#' @param supplied_witnesses Previously found vectors, rechecked for this problem.
#' @return A robust_rv_result with RV_lower, RV_upper, rv (only if resolved),
#'   a witness vector, its refit, and diagnostic status. Bounds are computational,
#'   not confidence intervals. Local failure never proves an empty loss set.
#' @export
rv_sand <- function(model, treatment = NULL, vcov = "HC1", cluster = NULL,
                    alpha = 0.05, critical = NULL, null = 0, n_starts = 8L,
                    max_iterations = 80L, seed = 1234L, tolerance = 5e-8,
                    max_radius = 1 - 1e-6, search_steps = 12L,
                    run_search = TRUE, supplied_witnesses = list()) {
  ctx <- .rv_loss_context(model, treatment, vcov, cluster, alpha, critical, null)
  control <- .rv_loss_controls(n_starts, max_iterations, seed, tolerance,
                               max_radius, search_steps, run_search)
  .rv_loss_solve(ctx, "sand", control, supplied_witnesses = supplied_witnesses)
}

#' Increment or profile calibrated sandwich robustness value
#'
#' Increment calibration uses M = omega_full + max(0, max(delta_omega)).
#' Profile calibration uses a0*exp(k*min(L,beta*(sqrt(x)+sqrt(y)))).
#' Both impose only 0 < omega <= cap; no lower floor at one is imposed.
#' The same rows, covariance tuning and cutoff are held fixed throughout.
#'
#' @inheritParams rv_sand
#' @param calibration Increment, profile, or a directly supplied fixed cap.
#' @param benchmarks Observed-control pool, as in rv_benchmarks.
#' @param k Positive profile multiplier; one is the reference specification.
#' @param M Positive fixed ratio cap when calibration = "cap".
#' @param ... Numerical arguments passed to rv_sand's solver.
#' @return A robust_rv_result, including the complete benchmark table. An empty
#'   or incomplete pool gives CALIBRATION_UNAVAILABLE with missing bounds.
#' @export
rv_cal <- function(model, treatment = NULL, calibration = c("increment", "profile", "cap"),
                   benchmarks = NULL, vcov = "HC1", cluster = NULL, k = 1,
                   M = NULL, alpha = 0.05, critical = NULL, null = 0, ...) {
  calibration <- match.arg(calibration)
  ctx <- .rv_loss_context(model, treatment, vcov, cluster, alpha, critical, null)
  if (!.rv_is_scalar_number(k) || k <= 0) stop("k must be positive and finite", call. = FALSE)
  defaults <- list(n_starts = 8L, max_iterations = 80L, seed = 1234L,
    tolerance = 5e-8, max_radius = 1 - 1e-6, search_steps = 12L, run_search = TRUE)
  dots <- list(...)
  unknown <- setdiff(names(dots), c(names(defaults), "supplied_witnesses"))
  if (length(unknown) || anyDuplicated(names(dots)) ||
      (length(dots) && is.null(names(dots))))
    stop("unknown or unnamed solver arguments: ", paste(unknown, collapse = ", "), call. = FALSE)
  controls <- defaults
  controls[intersect(names(dots), names(defaults))] <- dots[intersect(names(dots), names(defaults))]
  control <- do.call(.rv_loss_controls, controls)
  benchmark <- NULL
  if (calibration == "cap") {
    if (!.rv_is_scalar_number(M) || M <= 0) stop("M must be positive and finite", call. = FALSE)
  } else {
    benchmark_model <- if (inherits(model, "rv_cache")) model$model else model
    if (is.null(benchmark_model)) stop("observed-control calibration requires an lm model", call. = FALSE)
    benchmark <- rv_benchmarks(benchmark_model, ctx$treatment, benchmarks, ctx$rule)
    available <- benchmark$complete && (calibration != "profile" || benchmark$profile_available)
    if (!available) {
      envelope <- .rv_loss_envelope(ctx)
      result <- .rv_loss_result(ctx, calibration, NA_real_, NULL, 0, 0L,
        control, envelope, function(x, y) NA_real_, status = "CALIBRATION_UNAVAILABLE",
        benchmark = benchmark)
      result$substatus <- if (calibration == "profile" && benchmark$complete)
        "PROFILE_REQUIRES_VALID_SCALAR_POOL" else benchmark$status
      return(result)
    }
    if (calibration == "profile" && (!is.finite(k * benchmark$L) ||
        log(.rv_loss_envelope(ctx)$a0) + k * benchmark$L > log(.Machine$double.xmax)))
      stop("profile cap overflows; use a smaller k", call. = FALSE)
  }
  .rv_loss_solve(ctx, calibration, control, M, benchmark, k,
                 supplied_witnesses = dots$supplied_witnesses %||% list())
}

#' Practitioner summary of conventional, unrestricted and calibrated RVs
#'
#' @inheritParams rv_cal
#' @param ... Numerical solver arguments common to all sandwich calculations.
#' @return A sandwich_rv_analysis with a comparison table, complete results,
#'   benchmark records, and conventional references. RV_CH in the table uses
#'   the SAME cutoff as the sandwich test; RV_CH_standard is also retained.
#' @export
sandwich_rv <- function(model, treatment, benchmarks = NULL, vcov = "HC1",
                        cluster = NULL, k = 1, alpha = 0.05, critical = NULL,
                        null = 0, ...) {
  ctx <- .rv_loss_context(model, treatment, vcov, cluster, alpha, critical, null)
  delta <- rv_cal(model, treatment, "increment", benchmarks, ctx$rule,
                   alpha = alpha, critical = ctx$critical, null = null, ...)
  profile <- rv_cal(model, treatment, "profile", benchmarks, ctx$rule,
                     k = k, alpha = alpha, critical = ctx$critical, null = null, ...)
  dots <- list(...)
  witnesses <- c(dots$supplied_witnesses %||% list(),
    Filter(Negate(is.null), list(delta$best_witness_z, profile$best_witness_z)))
  dots$supplied_witnesses <- witnesses
  sand <- do.call(rv_sand, c(list(model = model, treatment = treatment, vcov = ctx$rule,
    cluster = ctx$cluster, alpha = alpha, critical = ctx$critical, null = null), dots))
  results <- list(RV_sand = sand, RV_cal_delta = delta, RV_cal_profile = profile)
  ch <- conventional_rv(ctx$cache, alpha, null, critical = ctx$critical)
  table <- data.frame(measure = "RV_CH", estimate = ch, lower = ch, upper = ch,
    status = "ANALYTICAL", benchmark = NA_character_, critical = ctx$critical, stringsAsFactors = FALSE)
  for (result in results) table <- rbind(table, data.frame(measure = result$measure,
    estimate = result$rv, lower = result$RV_lower, upper = result$RV_upper,
    status = result$status, benchmark = result$benchmark_status,
    critical = ctx$critical, stringsAsFactors = FALSE))
  structure(list(table = table, results = results, benchmark = delta$benchmark,
    k = k, RV_CH = ch, RV_CH_standard = conventional_rv(ctx$cache, alpha, null),
    call = match.call()), class = "sandwich_rv_analysis")
}

#' @export
print.sandwich_rv_analysis <- function(x, ...) {
  cat("Robustness values for loss of significance (partial R-squared scale)\n")
  print(x$table, row.names = FALSE, digits = 6)
  cat("Profile multiplier k =", x$k, "\n")
  cat("Bounds are computational intervals; see status for unavailable or unresolved values.\n")
  invisible(x)
}

#' @export
summary.sandwich_rv_analysis <- function(object, ...) object$table

#' @export
as.data.frame.sandwich_rv_analysis <- function(x, row.names = NULL,
                                             optional = FALSE, ...) x$table
