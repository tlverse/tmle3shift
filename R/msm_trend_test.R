#' Test for a trend in the effect of shift interventions via working MSM
#'
#' @param tmle_fit_estimates A \code{list} corresponding to the
#'  \code{$estimates} slot of an object of class \code{tmle3_Fit}, containing
#'  estimates of a grid of posited shift interventions.
#' @param delta_grid A \code{numeric} vector giving the individual values of
#'  the shift parameter used in computing each of the TML estimates.
#' @param level The nominal coverage probability of the confidence interval.
#' @param weighting A \code{character} indicating the type of weighting used
#'  for construction of the marginal structural model. \code{"identity"}
#'  applies the same weight to all individual estimates while \code{"variance"}
#'  applies weights based on the inverse variance of the estimate. It would be
#'  expected that variance-based weighting would yield more stable estimates of
#'  the parameter of the MSM. The default is identity-based weighting.
#'
#' @importFrom stats cov qnorm pnorm
#' @importFrom methods is
#' @importFrom assertthat assert_that
#'
#' @export
trend_msm <- function(tmle_fit_estimates, delta_grid, level = 0.95,
                      weighting = c("identity", "variance")) {

  # set default weighting to identity
  weighting <- match.arg(weighting)

  # make sure more than one parameter has been estimated for trend
  assert_that(length(tmle_fit_estimates) > 1)

  # matrix of EIF(O_i) values and estimates across each parameter estimated
  eif_mat <- sapply(tmle_fit_estimates, `[[`, "IC")
  psi_vec <- sapply(tmle_fit_estimates, `[[`, "psi")

  # set weights to be the inverse of the variance of each TML estimate
  if (weighting == "identity") {
    weights <- rep(1, ncol(eif_mat))
  } else if (weighting == "variance") {
    weights <- as.numeric(1 / diag(stats::cov(eif_mat)))
  }

  # multiplier for CI construction
  ci_mult <- (c(1, -1) * stats::qnorm((1 - level) / 2))

  # compute the MSM parameters
  intercept <- rep(1, length(delta_grid))
  x_mat <- cbind(intercept, delta_grid)
  omega <- diag(weights)
  s_mat <- tcrossprod(solve(crossprod(x_mat, omega) %*% x_mat), x_mat) %*% omega
  msm_param <- as.vector(s_mat %*% psi_vec)

  # compute inference for MSM based on individual EIF(O_i) for each parameter
  msm_eif <- t(tcrossprod(s_mat, eif_mat))
  msm_var <- diag(stats::cov(msm_eif))
  msm_se <- sqrt(msm_var / nrow(msm_eif))

  # build confidence intervals and hypothesis tests for EIF(msm)
  ci_msm_param <- tcrossprod(msm_se, ci_mult) + msm_param
  pval_msm_param <- 2 * stats::pnorm(-abs(msm_param / msm_se))

  # matrix for output
  out <- cbind(
    ci_msm_param[, 1], msm_param, ci_msm_param[, 2], msm_se,
    pval_msm_param
  )
  colnames(out) <- c("ci_low", "param_est", "ci_high", "param_se", "p_value")
  rownames(out) <- names(msm_se)
  return(out)
}
