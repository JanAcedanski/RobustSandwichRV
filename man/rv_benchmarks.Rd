\name{rv_benchmarks}
\alias{rv_benchmarks}
\title{Restoring-increment and profile benchmarks from observed controls}
\usage{
rv_benchmarks(model, treatment, benchmarks = NULL, vcov = "HC1", cluster = NULL)
}
\arguments{
\item{model}{An unweighted full-rank \code{lm} model with an intercept.}
\item{treatment}{Target model-matrix column name.}
\item{benchmarks}{Character vector of separate controls, or a named list of whole blocks. NULL selects eligible scalar terms.}
\item{vcov}{HC0/HC1, one-way CR0/CR1, HAC_Bartlett_L4, HAC_QS_BW5, TWCR0_raw, TWCR1_raw, or a fixed covariance rule.}
\item{cluster}{One- or two-way cluster labels aligned with the fitted sample, or a formula.}
}
\description{Deletes each prescribed control unit and refits the same target on the same rows. Records signed changes in the target's sandwich-to-conventional SE ratio and explanatory strengths in the reduced model.}
\details{For increment calibration, M equals omega in the full model plus the largest positive restoring increment, or zero added increment if all are nonpositive. Whole blocks are deleted together and their changes are not divided by rank. Scalar and block pools must be submitted separately.

For profile calibration, the reduced ratio is adjusted to the full model's chi correction, chi(v) = v*a(v). L is the positive maximum log increment; beta is the maximum positive log increment divided by sqrt(x)+sqrt(y), measured in the reduced model. The profile requires scalar units. A zero denominator and zero increment contributes zero; nonzero change at zero strength makes the profile unavailable.

Empty and incomplete pools are unavailable. Negative changes remain in the records. Variances must be positive and finite; a reduced SE ratio at or below 100 times machine epsilon is treated as numerical zero. Target and intercept deletion, overlapping units and mixed scalar/block pools are rejected. Deletions use stored model rows and matrices, avoiding changes to missing-value selection.}
\value{An \code{rv_calibration} with \code{records}, pool status, \code{M_delta}, \code{L}, \code{beta}, zero-benchmark indicators through the status and parameter values, and binding control names.}
\seealso{\code{\link{rv_cal}}}
\examples{
fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
rv_benchmarks(fit, "wt", c("hp", "disp"))
}
