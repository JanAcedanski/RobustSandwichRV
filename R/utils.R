.rv_supported_vcov <- c("classical", "HC0", "HC1", "CR0", "CR1")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.rv_match_vcov <- function(vcov) {
  vcov <- match.arg(vcov, .rv_supported_vcov)
  vcov
}

.rv_is_scalar_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

.rv_assert_probability <- function(x, name, upper_open = FALSE) {
  upper_ok <- if (upper_open) x < 1 else x <= 1
  if (!.rv_is_scalar_number(x) || x < 0 || !upper_ok) {
    stop(name, " must lie in [0, 1", if (upper_open) ")" else "]",
         call. = FALSE)
  }
  invisible(x)
}

.rv_as_numeric_matrix <- function(x, name) {
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " contains non-finite values", call. = FALSE)
  }
  x
}

.rv_qr_basis <- function(A, tol = NULL) {
  A <- .rv_as_numeric_matrix(A, "A")
  n <- nrow(A)
  if (ncol(A) == 0L) {
    return(matrix(numeric(), nrow = n, ncol = 0L))
  }
  if (is.null(tol)) {
    scale <- max(abs(A), 1)
    tol <- max(dim(A)) * .Machine$double.eps * scale
  }
  fit <- qr(A, tol = tol, LAPACK = FALSE)
  rank <- fit$rank
  if (rank == 0L) {
    return(matrix(numeric(), nrow = n, ncol = 0L))
  }
  qr.Q(fit, complete = FALSE)[, seq_len(rank), drop = FALSE]
}

.rv_residualize <- function(Q, A) {
  if (ncol(Q) == 0L) {
    return(A + 0)
  }
  A - Q %*% crossprod(Q, A)
}

.rv_dot <- function(x, y) {
  sum(x * y)
}

.rv_norm <- function(x) {
  sqrt(sum(x * x))
}

.rv_project_complement <- function(cache, u) {
  u <- as.numeric(u)
  if (length(u) != cache$n) {
    stop("u has the wrong length", call. = FALSE)
  }
  v <- as.numeric(.rv_residualize(cache$Qx, u))
  d_hat <- cache$d / .rv_norm(cache$d)
  e_norm <- .rv_norm(cache$e_short)
  if (!(e_norm > 0)) {
    stop("the short regression has zero residual variation", call. = FALSE)
  }
  e_hat <- cache$e_short / e_norm
  v <- v - d_hat * .rv_dot(d_hat, v)
  v <- v - e_hat * .rv_dot(e_hat, v)
  v
}

.rv_normalize_direction <- function(cache, u) {
  v <- .rv_project_complement(cache, u)
  nv <- .rv_norm(v)
  threshold <- 100 * .Machine$double.eps * max(.rv_norm(u), 1)
  if (!(nv > threshold)) {
    stop("no usable direction orthogonal to X, d, and the short residual",
         call. = FALSE)
  }
  v / nv
}

.rv_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.rv_structured_directions <- function(cache, n_starts = 12L, seed = 1234L,
                                      warm = list()) {
  directions <- list()
  add <- function(candidate) {
    value <- tryCatch(.rv_normalize_direction(cache, candidate),
                      error = function(e) NULL)
    if (!is.null(value)) {
      directions[[length(directions) + 1L]] <<- value
    }
  }
  for (candidate in warm) {
    add(candidate)
  }
  d <- cache$d
  e <- cache$e_short
  for (candidate in list(d * e, d^2 * e, d * e^2, sign(d * e),
                         abs(d) * e, d * abs(e))) {
    add(candidate)
  }
  .rv_with_seed(seed, {
    attempts <- 0L
    target <- max(as.integer(n_starts), 8L)
    while (length(directions) < target && attempts < 20L * target) {
      add(stats::rnorm(cache$n))
      attempts <- attempts + 1L
    }
  })
  if (length(directions) == 0L) {
    stop("the residual-direction sphere is empty", call. = FALSE)
  }
  directions
}

.rv_cluster_codes <- function(cluster, n) {
  if (is.null(cluster)) {
    stop("cluster labels are required for CR0 or CR1", call. = FALSE)
  }
  if (length(cluster) != n) {
    stop("cluster must have one label per retained observation", call. = FALSE)
  }
  if (anyNA(cluster)) {
    stop("cluster labels cannot be missing", call. = FALSE)
  }
  factor_cluster <- factor(cluster, exclude = NULL)
  codes <- as.integer(factor_cluster)
  G <- nlevels(factor_cluster)
  if (G <= 1L) {
    stop("cluster-robust covariance requires at least two clusters",
         call. = FALSE)
  }
  list(codes = codes, G = G, labels = factor_cluster)
}

.rv_reject <- function(t, critical) {
  abs(t) >= critical
}

.rv_decision_change <- function(t_short, t_full, critical) {
  xor(.rv_reject(t_short, critical), .rv_reject(t_full, critical))
}

.rv_decision_violation <- function(reject_short, t_full, critical) {
  if (reject_short) {
    max(abs(t_full) - critical, 0)
  } else {
    max(critical - abs(t_full), 0)
  }
}

.rv_critical_value <- function(cache, vcov, block_rank = 1L,
                               alpha = 0.05, cluster = NULL,
                               critical = NULL) {
  if (!is.null(critical)) {
    if (!.rv_is_scalar_number(critical) || critical <= 0) {
      stop("critical must be positive and finite", call. = FALSE)
    }
    return(as.numeric(critical))
  }
  if (!.rv_is_scalar_number(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie in (0, 1)", call. = FALSE)
  }
  vcov <- .rv_match_vcov(vcov)
  if (vcov %in% c("CR0", "CR1")) {
    G <- .rv_cluster_codes(cluster, cache$n)$G
    return(stats::qt(1 - alpha / 2, df = G - 1L))
  }
  df_full <- cache$df_short - as.integer(block_rank)
  if (df_full <= 0L) {
    stop("the candidate model has no residual degrees of freedom",
         call. = FALSE)
  }
  stats::qt(1 - alpha / 2, df = df_full)
}

.rv_timer <- function() {
  proc.time()[["elapsed"]]
}

