################################################################################
# MARGINAL STRUCTURAL MODELS
################################################################################

f_msm_linear <- function(psis, eifs, x) {
  # matrix of EIF(O_i) values and estimates across each parameter estimated
  eif_mat <- do.call(cbind, eifs)
  psi_vec <- do.call(cbind, psis)

  # set weights to be the inverse of the variance of each TML estimate
  if (is.null(weights)) {
    weights <- as.numeric(1 / diag(stats::cov(eif_mat)))
  }

  # compute the MSM parameters
  intercept <- rep(1, length(psi_vec))
  x_mat <- cbind(intercept, delta_grid)
  omega <- diag(weights)
  s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega
  msm_param <- as.vector(s_mat %*% psi_vec)
  return(msm_param)
}

df_msm_linear <- function(psis, eifs) {
  # matrix of EIF(O_i) values and estimates across each parameter estimated
  eif_mat <- do.call(cbind, eifs)
  psi_vec <- do.call(cbind, psis)

  # set weights to be the inverse of the variance of each TML estimate
  if (is.null(weights)) {
    weights <- as.numeric(1 / diag(stats::cov(eif_mat)))
  }

  # compute the MSM parameters
  intercept <- rep(1, length(psi_vec))
  x_mat <- cbind(intercept, delta_grid)
  omega <- diag(weights)
  s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega

  # compute inference for MSM based on individual EIF(O_i) for each parameter
  msm_eif <- t(tcrossprod(s_mat, eif_mat))
  return(msm_eif)
}

#' Linear Marginal Structural Models
#' @export
delta_param_MSM_linear <- list(type = "MSM_linear",
                               name = function(names) {
                                 sprintf("RR(%s/%s)", names[[2]], names[[1]])
                                },
                               f = f_msm_linear,
                               df = df_msm_linear
                              )

