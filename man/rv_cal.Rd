\name{rv_cal}
\alias{rv_cal}
\title{Increment and strength-dependent profile calibrated robustness values}
\usage{
rv_cal(model, treatment = NULL,
  calibration = c("increment", "profile", "cap"), benchmarks = NULL,
  vcov = "HC1", cluster = NULL, k = 1, M = NULL,
  alpha = 0.05, critical = NULL, null = 0, ...)
}
\arguments{
\item{model}{An unweighted \code{lm} model, or an \code{rv_cache}. Observed-control calibration requires stored \code{lm} metadata.}
\item{treatment}{Target model-matrix column.}
\item{calibration}{Restoring-increment, strength-dependent profile, or direct fixed cap.}
\item{benchmarks}{Separate control units as in \code{\link{rv_benchmarks}}.}
\item{vcov}{HC0/HC1, one-way CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw, TWCR1_raw, or a fixed covariance rule.}
\item{cluster}{One- or two-way cluster labels aligned with fitted rows, or a formula.}
\item{k}{Positive profile multiplier. The reference is one; two is a sensitivity scenario.}
\item{M}{Positive fixed cap when calibration is \code{"cap"}.}
\item{alpha}{Two-sided test level.}
\item{critical}{Optional fixed test cutoff; see \code{\link{rv_sand}} for defaults.}
\item{null}{Null coefficient value.}
\item{...}{Named numerical arguments from \code{\link{rv_sand}}, including previously found witnesses.}
}
\description{Computes the weakest added scalar index satisfying 0 < omega <= the calibrated cap that removes significance in the maintained full model.}
\details{Increment calibration uses \code{M_delta} from \code{rv_benchmarks}. Profile calibration imposes a0*exp(k*min(L,beta*(sqrt(x)+sqrt(y)))); a0 is the zero-strength augmented-model ratio, with the actual finite-sample correction. Both permit ratios below one. A complete zero profile benchmark gives the constant cap a0 for every k. An empty or invalid pool is never replaced by a zero benchmark.

The profile lower bound uses the crossing of the protected CH bound and the envelope a0*exp(k*min(L,2*beta*sqrt(r))). The unrestricted analytical envelope can strengthen both calibrated lower bounds. Endpoint and fallback candidates are tested against their actual pair-specific cap and independently refitted. Failure to attain the relaxation returns an interval or lower bound only. These are sensitivity assumptions about omitted variables, not consequences of the observed controls.}
\value{A \code{robust_rv_result}, including the complete benchmark records and profile parameters when applicable. \code{CALIBRATION_UNAVAILABLE} has missing bounds. For an already nonrejecting observed test and an available calibration, loss RV is zero by convention. See \code{\link{rv_sand}} for numerical and inferential interpretation.}
\seealso{\code{\link{sandwich_rv}}, \code{\link{rv_benchmarks}}}
\examples{
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
rv_cal(fit, "wt", "increment", c("hp", "disp"), run_search = FALSE)
rv_cal(fit, "wt", "profile", c("hp", "disp"), k = 1, run_search = FALSE)
}
