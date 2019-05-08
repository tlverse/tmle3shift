#' Defines a TML Estimator for Variable Importance for Continuous Interventions
#'
#' Current limitations: pretty much tailored to \code{Param_TSM}
#' See TODO notes for places generalization can be added
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'  Param_TSM
#'
#' @export
#
tmle3_Spec_vimshift_delta <- R6::R6Class(
  classname = "tmle3_Spec_vimshift_delta",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec_shift,
  public = list(
    initialize = function(shift_fxn = shift_additive_bounded,
                              shift_fxn_inv = shift_additive_bounded_inv,
                              shift_grid = seq(-1, 1, by = 0.5),
                              max_shifted_ratio = 2,
                              ...) {
      options <- list(
        shift_fxn = shift_fxn,
        shift_fxn_inv = shift_fxn_inv,
        shift_grid = shift_grid,
        max_shifted_ratio = max_shifted_ratio,
        ...
      )
      do.call(super$initialize, options)
    },
    make_params = function(tmle_task, likelihood) {
      # unwrap internalized arguments
      shift_fxn <- self$options$shift_fxn
      shift_fxn_inv <- self$options$shift_fxn_inv
      shift_grid <- self$options$shift_grid
      max_shifted_ratio <- self$options$max_shifted_ratio

      # treatment likelihood bound (away from 0 for continuous a)
      A_bound <- c(1 / tmle_task$nrow, Inf)

      # define shift intervention over grid (additive only for now)
      interventions <-
        lapply(shift_grid, function(x) {
          tmle3::define_lf(LF_shift,
            name = "A",
            original_lf = likelihood$factor_list[["A"]],
            likelihood_base = likelihood,                  # initial likelihood
            shift_fxn, shift_fxn_inv,                      # shift functions
            shift_delta = x,                               # shift value in grid
            max_shifted_ratio = max_shifted_ratio,         # ratio for shifting
            bound = A_bound                                # bound shifted g
          )
        })

      # create list of parameters (counterfactual treatment-specific means)
      tsm_params_list <-
        lapply(interventions, function(x) {
          tmle3::Param_TSM$new(likelihood, x)
        })

      # MSM function factory
      design_matrix <- cbind(rep(1, length(shift_grid)), shift_grid)
      colnames(design_matrix) <- c("intercept", "slope")
      delta_param_msm <- msm_linear_factory(design_matrix)

      # create MSM via delta method
      msm <- Param_delta$new(
        likelihood, delta_param_msm,
        tsm_params_list
      )
      tmle_params <- unlist(list(tsm_params_list, msm), recursive = FALSE)

      # output should be a list
      return(tmle_params)
    },
    make_updater = function() {
      updater <- tmle3_Update$new(cvtmle = TRUE)
    }
  ),
  active = list(),
  private = list()
)

################################################################################

#' Outcome Under a Grid of Shifted Interventions via Delta Method
#'
#' O = (W, A, Y)
#' W = Covariates
#' A = Treatment (binary or categorical)
#' Y = Outcome (binary or bounded continuous)
#'
#' @param shift_fxn A \code{function} defining the type of shift to be applied
#'  to the treatment. For an example, see \code{shift_additive}.
#' @param shift_fxn_inv A \code{function} defining the inverse of the type of
#'  shift to be applied to the treatment. For an example, see
#'  \code{shift_additive_inv}.
#' @param shift_grid A \code{numeric} vector, specification of a selection of
#'  shifts (on the level of the treatment) to be applied to the intervention.
#'  This is a value passed to the \code{function}s above for computing various
#'  values of the outcome under modulated values of the treatment.
#' @param max_shifted_ratio A \code{numeric} value indicating the maximum
#'  tolerance for the ratio of the counterfactual and observed intervention
#'  densities. In particular, the shifted value of the intervention is assigned
#'  to a given observational unit when the ratio of counterfactual intervention
#'  density to the observed intervention density is below this value.
#' @param ... Additional arguments, passed to shift functions.
#'
#' @importFrom sl3 make_learner Lrnr_mean
#'
#' @export
#
tmle_vimshift_delta <- function(shift_fxn = shift_additive_bounded,
                                shift_fxn_inv = shift_additive_bounded_inv,
                                shift_grid = seq(-1, 1, by = 0.5),
                                max_shifted_ratio = 2,
                                ...) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_vimshift_delta$new(
    shift_fxn, shift_fxn_inv,
    shift_grid, max_shifted_ratio,
    ...
  )
}
