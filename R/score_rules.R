# Immutable-by-convention, fixed score operators. The scalar correction is
# evaluated at the actual residual df; order, kernel and cluster weights stay
# fixed across control deletions and added-index candidates.
.rv_new_rule <- function(name, n, apply, correction, psd, norm_bound,
                          positive_bound = norm_bound, critical_df = NULL,
                          metadata = list()) {
  structure(list(name = name, n = n, apply = apply, correction = correction,
    psd = psd, norm_bound = norm_bound, positive_bound = positive_bound,
    critical_df = critical_df %||% function(df) df, metadata = metadata),
    class = "rv_score_rule")
}

.rv_rule_correction <- function(rule, df) {
  if (!.rv_is_scalar_number(df) || df <= 0) stop("residual df must be positive", call. = FALSE)
  value <- rule$correction(df)
  if (!.rv_is_scalar_number(value) || value <= 0)
    stop("the covariance correction must be positive and finite", call. = FALSE)
  value
}

.rv_group_apply <- function(value, codes) {
  totals <- rowsum(as.matrix(value), codes, reorder = TRUE)
  result <- totals[codes, , drop = FALSE]
  if (is.null(dim(value))) as.numeric(result) else result
}

#' Fixed-kernel heteroskedasticity and autocorrelation covariance rule
#'
#' @param n Number of retained observations.
#' @param kernel Bartlett or quadratic spectral (also accepts QS).
#' @param bandwidth For Bartlett, an integer maximum lag; for QS a fixed positive bandwidth.
#' @param order Permutation giving the temporal order of retained rows.
#' @param adjust Apply the HC1-style n/df correction if TRUE.
#' @param critical_df Optional positive cutoff df or function of residual df.
#' @return An rv_score_rule accepted as the vcov argument of the current RV APIs.
#' @export
fixed_hac_rule <- function(n, kernel = "bartlett", bandwidth = 4,
                            order = seq_len(n), adjust = FALSE,
                            critical_df = NULL) {
  if (!.rv_is_scalar_number(n) || n < 2 || n != floor(n))
    stop("n must be an integer greater than one", call. = FALSE)
  n <- as.integer(n)
  if (!is.character(kernel) || length(kernel) != 1L || is.na(kernel))
    stop("kernel must be Bartlett or quadratic_spectral", call. = FALSE)
  kernel <- tolower(kernel)
  if (kernel %in% c("qs", "quadratic spectral")) kernel <- "quadratic_spectral"
  if (!(kernel %in% c("bartlett", "quadratic_spectral")))
    stop("kernel must be Bartlett or quadratic_spectral", call. = FALSE)
  if (!.rv_is_scalar_number(bandwidth) || bandwidth <= 0)
    stop("bandwidth must be positive and finite", call. = FALSE)
  if (kernel == "bartlett" && (bandwidth != floor(bandwidth) || bandwidth >= n))
    stop("Bartlett bandwidth must be an integer maximum lag smaller than n", call. = FALSE)
  if (length(order) != n || anyNA(order) || !identical(sort(as.numeric(order)), as.numeric(seq_len(n))))
    stop("order must be a permutation of 1:n", call. = FALSE)
  if (!is.logical(adjust) || length(adjust) != 1L || is.na(adjust))
    stop("adjust must be TRUE or FALSE", call. = FALSE)
  permutation <- as.integer(order); inverse <- order(permutation)
  lag <- seq_len(n - 1L)
  weights <- if (kernel == "bartlett") pmax(0, 1 - lag / (bandwidth + 1)) else {
    z <- 6 * pi * lag / (5 * bandwidth)
    ifelse(abs(z) < 1e-3, 1 - z^2 / 10 + z^4 / 280 - z^6 / 15120,
           3 / z^2 * (sin(z) / z - cos(z)))
  }
  # Dense BLAS is efficient for the article's moderate time-series sizes.
  # Larger inputs use the identical fixed weights without an n-by-n matrix.
  matrix <- if (n <= 2000L) stats::toeplitz(c(1, weights))[inverse, inverse] else NULL
  apply <- if (!is.null(matrix)) function(value) {
    answer <- matrix %*% value
    if (is.null(dim(value))) as.numeric(answer) else answer
  } else function(value) {
    v <- as.matrix(value)[permutation, , drop = FALSE]; answer <- v
    for (h in which(weights != 0)) {
      i <- seq_len(n - h); j <- i + h
      answer[i, ] <- answer[i, , drop = FALSE] + weights[h] * v[j, , drop = FALSE]
      answer[j, ] <- answer[j, , drop = FALSE] + weights[h] * v[i, , drop = FALSE]
    }
    answer <- answer[inverse, , drop = FALSE]
    if (is.null(dim(value))) as.numeric(answer) else answer
  }
  name <- if (kernel == "bartlett") paste0("HAC_Bartlett_L", bandwidth) else paste0("HAC_QS_BW", bandwidth)
  if (adjust) name <- paste0(name, "_adjusted")
  .rv_new_rule(name, n, apply, if (adjust) function(df) n / df else function(df) 1,
    TRUE, 1 + 2 * sum(abs(weights)), critical_df = critical_df,
    metadata = list(family = "HAC", kernel = kernel, bandwidth = bandwidth,
      order = permutation, weights = weights, adjust = adjust,
      prewhite = FALSE, automatic_bandwidth = FALSE, psd_repair = FALSE))
}

#' Fixed raw two-way cluster covariance rule
#'
#' @param cluster1,cluster2 Labels aligned with the retained regression rows.
#' @param variant CR0 or CR1; CR1 uses a separate G/(G-1) in each component.
#' @param critical_df Optional cutoff df; default min(G1,G2)-1.
#' @return A fixed rv_score_rule for A1 + A2 - A12, without PSD repair.
#' @export
fixed_twoway_cluster_rule <- function(cluster1, cluster2, variant = c("CR1", "CR0"),
                                       critical_df = NULL) {
  variant <- match.arg(variant)
  n <- length(cluster1)
  one <- .rv_cluster_codes(cluster1, n); two <- .rv_cluster_codes(cluster2, n)
  pair <- .rv_cluster_codes(interaction(one$codes, two$codes, drop = TRUE), n)
  groups <- c(one$G, two$G, pair$G)
  multipliers <- if (variant == "CR1") groups / (groups - 1) else rep(1, 3)
  codes <- list(one$codes, two$codes, pair$codes)
  apply <- function(value) multipliers[1] * .rv_group_apply(value, codes[[1]]) +
    multipliers[2] * .rv_group_apply(value, codes[[2]]) -
    multipliers[3] * .rv_group_apply(value, codes[[3]])
  bounds <- multipliers * vapply(codes, function(x) max(tabulate(x)), numeric(1))
  .rv_new_rule(paste0("TW", variant, "_raw"), n, apply,
    if (variant == "CR1") function(df) (n - 1) / df else function(df) 1,
    FALSE, sum(bounds), positive_bound = sum(bounds[1:2]),
    critical_df = critical_df %||% (min(groups[1:2]) - 1),
    metadata = list(family = "two_way_cluster", variant = variant, groups = groups,
      codes = codes, component_multipliers = multipliers, psd_repair = FALSE))
}

.rv_resolve_rule <- function(vcov, n, cluster = NULL) {
  vcov <- .rv_match_vcov(vcov)
  if (inherits(vcov, "rv_score_rule")) {
    if (vcov$n != n) stop("fixed covariance rule has a different observation count", call. = FALSE)
    return(vcov)
  }
  if (vcov == "HAC_Bartlett_L4") return(fixed_hac_rule(n, "bartlett", 4))
  if (vcov == "HAC_QS_BW5") return(fixed_hac_rule(n, "quadratic_spectral", 5))
  if (vcov %in% c("TWCR0_raw", "TWCR1_raw")) {
    if (!(is.data.frame(cluster) || is.matrix(cluster)) || ncol(cluster) != 2L || nrow(cluster) != n)
      stop("two-way clustering requires a two-column cluster data frame or matrix", call. = FALSE)
    return(fixed_twoway_cluster_rule(cluster[, 1L], cluster[, 2L],
      variant = if (vcov == "TWCR1_raw") "CR1" else "CR0"))
  }
  if (vcov %in% c("HC0", "HC1")) return(.rv_new_rule(vcov, n, function(value) value,
    if (vcov == "HC1") function(df) n / df else function(df) 1, TRUE, 1,
    metadata = list(family = "HC")))
  if (vcov %in% c("CR0", "CR1")) {
    info <- .rv_cluster_codes(cluster, n)
    .rv_new_rule(vcov, n, function(value) .rv_group_apply(value, info$codes),
      if (vcov == "CR1") function(df) info$G / (info$G - 1) * (n - 1) / df else function(df) 1,
      TRUE, max(tabulate(info$codes)), critical_df = info$G - 1,
      metadata = list(family = "one_way_cluster", groups = info$G, codes = info$codes))
  } else stop("a sandwich covariance rule is required", call. = FALSE)
}

.rv_rule_df <- function(rule, df) {
  value <- if (is.function(rule$critical_df)) rule$critical_df(df) else rule$critical_df
  if (!.rv_is_scalar_number(value) || value <= 0)
    stop("critical degrees of freedom must be positive and finite", call. = FALSE)
  value
}

#' @export
print.rv_score_rule <- function(x, ...) {
  cat("Fixed sandwich covariance rule:", x$name, "\n")
  cat("  Observations:", x$n, "\n")
  cat("  PSD repair: none\n")
  invisible(x)
}
