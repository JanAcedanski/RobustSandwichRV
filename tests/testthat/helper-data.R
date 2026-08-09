make_rv_fixture <- function(n = 96L, seed = 1701L, target_t = NULL) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  cluster <- rep(seq_len(12L), length.out = n)
  cluster_effect <- rnorm(12L, sd = 0.5)[cluster]
  d <- 0.35 * x1 - 0.2 * x2 + rnorm(n)
  error <- (0.7 + 0.35 * abs(d)) * rnorm(n) + cluster_effect
  y <- 0.55 * d + 0.25 * x1 - 0.1 * x2 + error
  model <- lm(y ~ d + x1 + x2)
  if (!is.null(target_t)) {
    cache0 <- cache_from_lm(model, "d")
    short0 <- RobustSandwichRV::evaluate_subspace(
      cache0, matrix(numeric(), nrow = n, ncol = 0L),
      vcov = "HC1", original_rank = 0L
    )
    shift <- (target_t - short0$t_robust) * short0$se_robust
    y <- y + shift * cache0$d
    model <- lm(y ~ d + x1 + x2)
  }
  list(
    n = n,
    model = model,
    cache = cache_from_lm(model, "d"),
    cluster = cluster,
    x1 = x1,
    x2 = x2,
    d = d,
    y = y
  )
}

