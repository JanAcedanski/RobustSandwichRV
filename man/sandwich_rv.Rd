\name{sandwich_rv}
\alias{sandwich_rv}
\alias{print.sandwich_rv_analysis}
\alias{summary.sandwich_rv_analysis}
\alias{as.data.frame.sandwich_rv_analysis}
\title{Compare CH, unrestricted sandwich and both calibrated robustness values}
\usage{
sandwich_rv(model, treatment, benchmarks = NULL, vcov = "HC1",
  cluster = NULL, k = 1, alpha = 0.05, critical = NULL, null = 0, ...)
\method{print}{sandwich_rv_analysis}(x, ...)
\method{summary}{sandwich_rv_analysis}(object, ...)
\method{as.data.frame}{sandwich_rv_analysis}(x, row.names = NULL,
  optional = FALSE, ...)
}
\arguments{
\item{model}{Unweighted \code{lm} with an intercept and full-rank design.}
\item{treatment}{One target model-matrix column.}
\item{benchmarks}{Observed control units; NULL selects scalar terms.}
\item{vcov}{HC0/HC1, one-way CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw, TWCR1_raw, or a fixed covariance rule.}
\item{cluster}{One- or two-way cluster labels aligned with fitted rows, or a formula.}
\item{k}{Positive profile multiplier; default one.}
\item{alpha}{Two-sided significance level.}
\item{critical}{Optional fixed cutoff common to the comparison.}
\item{null}{Null coefficient value.}
\item{...}{Numerical solver arguments from \code{\link{rv_sand}}; unused in methods.}
\item{x,object}{A \code{sandwich_rv_analysis}.}
\item{row.names,optional}{Standard data-frame method arguments.}
}
\description{Returns a practitioner table with RV_CH, RV_sand, RV_cal_delta and RV_cal_profile on the same partial-R-squared scale and test cutoff.}
\details{The calibrated searches are run first so their independently verified indices can also improve the unrestricted upper bound, after rechecking. The two calibrated classes are not generally nested. A point estimate is shown only when bounds close within tolerance. The table's CH row uses the sandwich cutoff; the conventional default-cutoff reference is retained separately as \code{RV_CH_standard}.}
\value{A \code{sandwich_rv_analysis} with \code{table}, three full \code{results}, benchmark diagnostics, \code{RV_CH}, \code{RV_CH_standard}, and k. \code{summary()} and \code{as.data.frame()} return the four-row table. Its benchmark column labels positive, zero and unavailable calibrations separately from computational status.}
\examples{
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
sandwich_rv(fit, "wt", c("hp", "disp"), n_starts = 4,
            max_iterations = 20, search_steps = 3)
}
