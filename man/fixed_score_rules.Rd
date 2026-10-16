\name{fixed_score_rules}
\alias{fixed_hac_rule}
\alias{fixed_twoway_cluster_rule}
\alias{print.rv_score_rule}
\title{Fixed HAC and raw two-way cluster covariance rules}
\usage{
fixed_hac_rule(n, kernel = "bartlett", bandwidth = 4,
  order = seq_len(n), adjust = FALSE, critical_df = NULL)
fixed_twoway_cluster_rule(cluster1, cluster2,
  variant = c("CR1", "CR0"), critical_df = NULL)
\method{print}{rv_score_rule}(x, ...)
}
\arguments{
\item{n}{Number of rows retained by the regression.}
\item{kernel}{\code{"bartlett"} or \code{"quadratic_spectral"}; \code{"QS"} is also accepted.}
\item{bandwidth}{For Bartlett, a positive integer maximum lag smaller than n; for QS, a fixed positive bandwidth.}
\item{order}{Permutation giving the temporal order of the retained rows.}
\item{adjust}{If TRUE, use an n/residual-df correction. FALSE matches the article's unadjusted HAC rules.}
\item{critical_df}{Optional positive cutoff degrees of freedom or function of residual df. Defaults to augmented residual df for HAC and min(G1,G2)-1 for two-way clustering. The resulting test cutoff is held fixed during RV search.}
\item{cluster1,cluster2}{Two cluster dimensions aligned with retained observations, with at least two groups in each.}
\item{variant}{CR0 or CR1.}
\item{x}{A fixed \code{rv_score_rule}.}
\item{...}{Unused.}
}
\description{Construct fixed score operators for use as the \code{vcov} argument in \code{sandwich_rv}, \code{rv_sand}, \code{rv_cal}, \code{rv_benchmarks} and the confounder evaluators.}
\details{HAC freezes row order, bandwidth and kernel weights. Bartlett maximum lag L uses weight 1-h/(L+1) at h<=L. QS uses 3/z^2*(sin(z)/z-cos(z)), z=6*pi*h/(5*bandwidth), at all retained lags. No prewhitening or automatic bandwidth selection is performed. The aliases \code{"HAC_Bartlett_L4"} and \code{"HAC_QS_BW5"} reproduce the article's settings in current fitted-row order.

Raw two-way clustering uses A1+A2-A12. For CR1, each component has its own Gj/(Gj-1) multiplier in the fixed operator, and the common (n-1)/residual-df multiplier changes with model rank. The intersection term is not replaced by HC0. There is no PSD repair. A nonpositive observed target variance is an error; a nonpositive candidate variance is inadmissible; an invalid deletion makes the prescribed calibration pool unavailable.

Alternatively pass \code{vcov="TWCR1_raw"} or \code{"TWCR0_raw"} with a two-column cluster data frame, matrix, or formula. A rule object embeds its order or cluster partitions, so subsequent cluster arguments are unnecessary. Rules are reused unchanged for every augmented and reduced model; only the specified scalar rank correction changes.}
\value{An \code{rv_score_rule} containing fixed weights or partitions, component corrections, critical-df convention and metadata. The complete rule is also stored as \code{result$covariance_rule} for independent refits.}
\seealso{\code{\link{sandwich_rv}}, \code{\link{validate_sandwich_reference}}}
\examples{
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
hac <- fixed_hac_rule(nobs(fit), "QS", bandwidth = 5)
rv_cal(fit, "wt", "profile", c("hp", "disp"), vcov = hac,
       n_starts = 4, max_iterations = 20, run_search = FALSE)
fixed_twoway_cluster_rule(rep(1:8, each = 4), rep(1:4, 8))
}
