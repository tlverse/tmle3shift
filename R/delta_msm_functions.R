#' Function factory for computing linear working MSM delta method functions
#'
#' @param design_matrix A \code{matrix} or \code{data.frame} giving the design
#'  matrix to be used in specifying terms of the linear working marginal
#'  structural model (MSM) for which parameters are to be estimated.
#' @param weighting A \code{character} indicating the type of weighting used for
#'  construction of the marginal structural model. \code{"identity"} applies the
#'  same weight to all individual estimates while \code{"variance"} applies
#'  weights based on the inverse variance of the estimate. It is expected that
#'  variance-based weighting would yield more stable estimates of the parameter
#'  of the MSM; however, the default remains the identity weighting.
#'
#' @keywords internal
#
msm_linear_factory <- function(design_matrix,
                               weighting = c("identity", "variance")) {
  # by default use identity weighting
  weighting <- match.arg(weighting)

  # NOTE: this is equivalent to a hat matrix for the linear working MSM
  get_msm_hat <- function(x, dx) {
    # vector of parameter estimates and matrix of EIF values
    psi_vec <- do.call(c, x)
    eif_mat <- do.call(cbind, dx)

    # set weights to be the inverse of the variance of each TML estimate
    if (weighting == "identity") {
      weights <- rep(1, ncol(eif_mat))
    } else if (weighting == "variance") {
      weights <- as.numeric(1 / diag(stats::cov(eif_mat)))
    }

    # compute the MSM parameters
    x_mat <- as.matrix(design_matrix)
    omega <- diag(weights)
    s_mat <- tcrossprod(solve(crossprod(x_mat, omega) %*%
      x_mat), x_mat) %*% omega
    return(s_mat)
  }

  # function for fitting parameter estimates of MSMs
  f_msm_linear <- function(x, dx) {
    # vector of parameter estimates
    psi_vec <- do.call(c, x)

    # get MSM projection matrix
    s_mat <- get_msm_hat(x = x, dx = dx)
    msm_param <- as.vector(s_mat %*% psi_vec)
    return(msm_param)
  }

  # function for computing EIF values of parameters from MSM
  df_msm_linear <- function(x, dx) {
    # matrix of EIF values
    eif_mat <- do.call(cbind, dx)

    # get MSM projection matrix
    s_mat <- get_msm_hat(x = x, dx = dx)

    # compute inference for MSM based on individual EIF(O_i) for each parameter
    msm_eif <- t(tcrossprod(s_mat, eif_mat))
    return(msm_eif)
  }

  # create list with the f and df functions for delta method
  delta_param_MSM_linear <- list(
    type = "MSM_linear",
    name = function(names) {
      sprintf("MSM(%s)", colnames(design_matrix))
    },
    f = f_msm_linear,
    df = df_msm_linear,
    hn = get_msm_hat
  )

  # output the list containing the f and df functions
  return(delta_param_MSM_linear)
}
