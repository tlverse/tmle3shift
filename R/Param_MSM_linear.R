#' Linear Working Marginal Structural Models
#'
#' Parameter definition for targeting the parameters of a linear working
#' marginal structural model (MSM): EY = beta0 + beta1 delta, to
#' summarize the variable importance results of a grid of shift interventions.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom tmle3 Param_base
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_MSM_linear, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood}}{the counterfactual likelihood for this treatment
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention
#'     }
#' }
#' @export
Param_MSM_linear <- R6Class(
  classname = "Param_MSM_linear",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3::Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list, ...,
                              outcome_node = "Y", shift_grid) {
      # initial
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.cf_likelihood <- make_CF_Likelihood(
        observed_likelihood,
        intervention_list
      )
      private$.shift_grid <- shift_grid
    },
    clever_covariates = function(tmle_task = NULL, cv_fold = -1) {
      # use training task if none provided
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- names(self$intervention_list)

      # create list of counterfactual means (parameters)
      tsm_cf_params <-
        lapply(self$intervention_list, function(x) {
          tmle3::Param_TSM$new(self$observed_likelihood, x)
        })

      browser()

      # combine clever covariates from individual parameters to target MSM
      tsm_aux_covars_Hn_list <-
        lapply(seq_along(tsm_cf_params), function(x) {
          as.numeric(tsm_cf_params[[x]]$clever_covariates(tmle_task)$Y)
        })
      tsm_aux_covars_Hn_mat <- do.call(cbind, tsm_aux_covars_Hn_list)

      # extract estimates of EIF in observed data
      tsm_eif_vals_list <-
        lapply(seq_along(tsm_cf_params), function(x) {
          as.numeric(tsm_cf_params[[x]]$estimates(tmle_task, cv_fold)$IC)
        })
      tsm_eif_vals_mat <- do.call(cbind, tsm_eif_vals_list)


      # set weights to be the inverse of the variance of each TML estimate
      wts <- as.numeric(1 / diag(stats::cov(tsm_eif_vals_mat)))
      omega <- diag(wts)

      # compute the MSM parameters
      intercept <- rep(1, length(self$shift_grid))
      x_mat <- cbind(intercept, self$shift_grid)
      s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega

      # build auxiliary covariates for each MSM parameter
      Hn_msm <- t(s_mat %*% t(tsm_aux_covars_Hn_mat))
      return(list(Y = Hn_msm))
    },
    estimates = function(tmle_task = NULL, cv_fold = -1) {
      # use training task if none provided
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }


      # build auxiliary covariates for each MSM parameter
      eif_msm <- t(s_mat %*% t(tsm_eif_mat))
      psi_msm <- as.numeric(apply(eif_msm, 2, sum))
      result <- list(psi = psi_msm, IC = eif_msm)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf(
        "E[%s_{%s}]", self$outcome_node,
        self$cf_likelihood$name
      )
      return(param_form)
    },
    intervention_list = function() {
      return(self$intervention_list)
    },
    shift_grid = function() {
      return(self$shift_grid)
    },
    update_nodes = function() {
      return(self$outcome_node)
    }
  ),
  private = list(
    .type = "MSM_linear",
    .cf_likelihood = NULL,
    .shift_grid = NULL
  )
)
