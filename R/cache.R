#' Prepare a linear-regression sensitivity cache
#'
#' Computes the Frisch-Waugh-Lovell quantities used throughout the package.
#' `X` should contain the intercept and all observed controls. Rank-deficient
#' columns are absorbed by a rank-revealing QR decomposition.
#'
#' @param y Numeric outcome vector.
#' @param d Numeric treatment vector.
#' @param X Numeric control matrix, normally including an intercept.
#' @return An object of class `rv_cache`.
#' @export
prepare_model <- function(y, d, X) {
  y <- as.numeric(y)
  d_raw <- as.numeric(d)
  X <- .rv_as_numeric_matrix(X, "X")
  n <- length(y)
  if (length(d_raw) != n || nrow(X) != n) {
    stop("y, d, and X must have the same number of observations",
         call. = FALSE)
  }
  if (n <= 2L) {
    stop("at least three observations are required", call. = FALSE)
  }
  if (any(!is.finite(y)) || any(!is.finite(d_raw))) {
    stop("y and d must contain only finite values", call. = FALSE)
  }
  Qx <- .rv_qr_basis(X)
  rank_x <- ncol(Qx)
  if (rank_x < ncol(X)) {
    warning("X is rank deficient; redundant columns were absorbed",
            call. = FALSE)
  }
  d <- as.numeric(.rv_residualize(Qx, d_raw))
  dd <- .rv_dot(d, d)
  if (!(dd > 100 * .Machine$double.eps * max(.rv_dot(d_raw, d_raw), 1))) {
    stop("d has no independent variation after conditioning on X",
         call. = FALSE)
  }
  y_x <- as.numeric(.rv_residualize(Qx, y))
  tau_short <- .rv_dot(d, y_x) / dd
  e_short <- y_x - tau_short * d
  df_short <- n - rank_x - 1L
  if (df_short <= 1L) {
    stop("one synthetic confounder would leave no residual degrees of freedom",
         call. = FALSE)
  }
  structure(
    list(
      y = y,
      d_raw = d_raw,
      X = X,
      Qx = Qx,
      d = d,
      y_x = y_x,
      e_short = e_short,
      tau_short = tau_short,
      rank_x = rank_x,
      df_short = df_short,
      n = n
    ),
    class = "rv_cache"
  )
}

.rv_extract_lm <- function(model, treatment) {
  if (!inherits(model, "lm")) {
    stop("model must inherit from 'lm'", call. = FALSE)
  }
  if (!is.null(model$weights)) {
    stop("weighted least squares is not supported in version 0.1.0",
         call. = FALSE)
  }
  if (!is.null(model$offset)) {
    stop("models with an offset are not supported in version 0.1.0",
         call. = FALSE)
  }
  mf <- stats::model.frame(model)
  mm <- stats::model.matrix(model)
  y <- stats::model.response(mf)
  if (!is.numeric(y)) {
    stop("the outcome must be numeric", call. = FALSE)
  }
  treatment <- as.character(treatment)
  if (length(treatment) != 1L) {
    stop("treatment must identify one model-matrix column", call. = FALSE)
  }
  index <- which(colnames(mm) == treatment)
  if (length(index) != 1L) {
    stop("treatment must exactly match one model-matrix column; available columns: ",
         paste(colnames(mm), collapse = ", "), call. = FALSE)
  }
  X <- mm[, -index, drop = FALSE]
  list(
    y = as.numeric(y),
    d = as.numeric(mm[, index]),
    X = X,
    model_matrix = mm,
    treatment_index = index,
    treatment = treatment,
    model_frame = mf,
    assign = attr(mm, "assign"),
    term_labels = attr(stats::terms(model), "term.labels")
  )
}

#' Create a sensitivity cache from an `lm` model
#'
#' @param model An unweighted [stats::lm()] object.
#' @param treatment Name of the single treatment model-matrix column.
#' @return An `rv_cache` with the original model metadata attached.
#' @export
cache_from_lm <- function(model, treatment) {
  extracted <- .rv_extract_lm(model, treatment)
  cache <- prepare_model(extracted$y, extracted$d, extracted$X)
  cache$model <- model
  cache$treatment <- extracted$treatment
  cache$model_frame <- extracted$model_frame
  cache$model_matrix <- extracted$model_matrix
  cache$model_assign <- extracted$assign
  cache$term_labels <- extracted$term_labels
  cache$treatment_index <- extracted$treatment_index
  cache
}

.rv_resolve_cluster <- function(cluster, cache) {
  if (is.null(cluster)) {
    return(NULL)
  }
  if (inherits(cluster, "formula")) {
    if (is.null(cache$model_frame)) {
      stop("a cluster formula requires a cache created from lm", call. = FALSE)
    }
    frame <- stats::model.frame(cluster, data = cache$model_frame,
                                na.action = stats::na.pass)
    if (ncol(frame) != 1L) {
      stop("cluster formula must select exactly one variable", call. = FALSE)
    }
    cluster <- frame[[1L]]
  }
  if (length(cluster) != cache$n) {
    stop("cluster must align with the observations retained by lm",
         call. = FALSE)
  }
  cluster
}

#' @export
print.rv_cache <- function(x, ...) {
  cat("Linear-regression robustness cache\n")
  cat("  Observations:", x$n, "\n")
  cat("  Control rank:", x$rank_x, "\n")
  cat("  Short-model residual df:", x$df_short, "\n")
  cat("  Short coefficient:", format(x$tau_short, digits = 7), "\n")
  invisible(x)
}
