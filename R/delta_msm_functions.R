################################################################################
# MARGINAL STRUCTURAL MODELS
################################################################################

# function factory for linear MSMs
msm_linear_factory <- function(params_list, design_matrix) {
  # bookkeeping

  psis <- lapply(estimates, `[[`, "psi")
  ICs <- lapply(estimates, `[[`, "IC")

  # ...
  f_msm_linear <- function(...) {
    ...
  }

  # ...
  df_msm_linear <- function(...) {
    ...
  }

  # create list with the f and df functions for delta method
  delta_param_MSM_linear <- list(type = "MSM_linear",
                                 name = "linear MSM via delta method",
                                 f = f_msm_linear,
                                 df = df_msm_linear
                                )

  # output the list containing the f and df functions
  return(delta_param_MSM_linear)
}

# compute parameters of working MSM via delta method
f_msm_linear <- function(psis, weights = NULL, delta_grid, ...) {
  # vector of parameter estimates
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

# compute EIFs of parameters of working MSM via delta method
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

#' Linear Working Marginal Structural Models
#' @export
delta_param_MSM_linear <- list(type = "MSM_linear",
                               name = function(names) {
                                 sprintf("MSM(%s/%s)", names[[2]], names[[1]])
                                },
                               f = f_msm_linear,
                               df = df_msm_linear
                              )

