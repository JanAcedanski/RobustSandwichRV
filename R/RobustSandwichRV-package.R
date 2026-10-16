#' RobustSandwichRV: certified sample robustness values
#'
#' A pure-R implementation of omitted-variable sensitivity analysis for
#' linear regression with HC0/HC1, one- or two-way CR0/CR1, and fixed-kernel
#' HAC standard errors.
#' The recommended workflow is sandwich_rv(), comparing RV_sand with increment
#' and strength-dependent profile calibration. Earlier point-conditioned and
#' symmetric-band procedures remain available under their original names.
#'
#' @importFrom utils head
#' @keywords internal
"_PACKAGE"
