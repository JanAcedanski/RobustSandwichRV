n <- 120L
i <- seq_len(n)
x1 <- sin(0.37 * i) + 0.15 * cos(0.11 * i)
x2 <- cos(0.23 * i) - 0.10 * sin(0.07 * i)
cluster <- (i - 1L) %% 12L + 1L
d <- 0.4 * x1 - 0.3 * x2 + sin(1.17 * i) + 0.1 * cos(0.03 * i^2)
z <- 0.25 * d + 0.35 * x1 + cos(0.73 * i) + 0.08 * sin(0.13 * i^2)
cluster_effect <- 0.35 * sin(0.8 * cluster)
error <- (0.8 + 0.25 * abs(d)) * cos(0.91 * i) + cluster_effect
y <- 0.6 * d + 0.2 * x1 - 0.1 * x2 + error

fixture <- data.frame(y, d, x1, x2, z, cluster)
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
options(digits = 17)
write.csv(
  fixture,
  "inst/extdata/julia_reference_input.csv",
  row.names = FALSE,
  quote = FALSE
)
