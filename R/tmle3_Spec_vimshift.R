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
tmle3_Spec_vimshift <- R6::R6Class(
  classname = "tmle3_Spec_vimshift",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(shift_fxn = shift_additive_bounded,
                          shift_fxn_inv = shift_additive_bounded_inv,
                          shift_grid = seq(-1, 1, by = 0.5),
                          ...) {
      options <- list(
        shift_fxn = shift_fxn,
        shift_fxn_inv = shift_fxn_inv,
        shift_grid = shift_grid
      )
      shift_args_extra = list(...)
      browser()
      do.call(super$initialize, options)
    },
    make_params = function(tmle_task, likelihood) {
      # TODO: export and use sl3:::get_levels
      A_vals <- tmle_task$get_tmle_node("A")
      if (is.factor(A_vals)) {
        msg <- paste(
          "This parameter is defined as a series of shifts of a continuous",
          "treatment. The treatment detected is NOT continuous."
        )
        stop(msg)
      }

      # unwrap internalized arguments
      shift_fxn <- self$options$shift_fxn
      shift_fxn_inv <- self$options$shift_fxn_inv
      delta_grid <- self$options$shift_grid

      # define shift intervention over grid (additive only for now)
      interventions <- lapply(delta_grid,
                              function(x) {
                                tmle3::define_lf(LF_shift,
                                                 name = "A",
                                                 original_lf =
                                                   likelihood$factor_list[["A"]],
                                                 shift_fxn,
                                                 shift_fxn_inv,
                                                 shift_delta = x)
                              })

      # compute mean counterfactual outcome under shifted interventions
      tmle_params <- lapply(interventions,
                            function(x) {
                              tmle3::Param_TSM$new(likelihood, x)
                            })
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

################################################################################

#' Outcome Under a Grid of Delta-Shifted Interventions
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
#' @param ... Additional arguments, passed to shift functions.
#'
#' @importFrom sl3 make_learner Lrnr_mean
#'
#' @export
#
tmle_vimshift <- function(shift_fxn = shift_additive_bounded,
                       shift_fxn_inv = shift_additive_bounded_inv,
                       shift_grid = seq(-1, 1, by = 0.5),
                       ...) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_vimshift$new(shift_fxn, shift_fxn_inv, shift_grid)
}

