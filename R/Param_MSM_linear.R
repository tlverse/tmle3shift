#' Parameter for Linear Working Marginal Structural Model
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
  inherit = tmle3::Param_delta,
  public = list(
    clever_covariates = function(tmle_task = NULL, cv_fold = -1) {
      # use training task if none provided
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- names(self$intervention_list)

      estimates <- lapply(
        self$parent_parameters,
        function(tmle_param) {
          tmle_param$estimates(tmle_task, cv_fold)
        }
      )

      psis <- lapply(estimates, `[[`, "psi")
      eifs <- lapply(estimates, `[[`, "IC")
      hn_msm_coef <- self$delta_param$hn(x = psis, dx = eifs)

      # combine clever covariates from individual parameters to target MSM
      hn_params_list <-
        lapply(self$parent_parameters, function(tmle_param) {
          as.numeric(tmle_param$clever_covariates(tmle_task, cv_fold)$Y)
        })
      hn_params_mat <- do.call(cbind, hn_params_list)

      # build auxiliary covariates for each MSM parameter
      hn_msm <- t(hn_msm_coef %*% t(hn_params_mat))
      return(list(Y = hn_msm))
    }
  )
)
