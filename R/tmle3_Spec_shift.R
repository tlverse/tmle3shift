#' Defines a TML Estimator for the Outcome under a Shifted Treatment
#'
#' Current limitations: pretty much tailored to \code{Param_TSM}.
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'  Param_TSM
#'
#' @export
#
tmle3_Spec_shift <- R6::R6Class(
  classname = "tmle3_Spec_shift",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(shift_fxn = shift_additive_bounded,
                              shift_fxn_inv = shift_additive_bounded_inv,
                              shift_val = 0,
                              max_shifted_ratio = 2,
                              ...) {
      options <- list(
        shift_fxn = shift_fxn,
        shift_fxn_inv = shift_fxn_inv,
        delta_shift = shift_val,
        max_shifted_ratio = max_shifted_ratio,
        ...
      )
      do.call(super$initialize, options)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # check that the intervention is continuous-valued
      A_type <- tmle_task$npsem[["A"]]$variable_type
      if (A_type$type != "continuous") {
        msg <- paste(
          "This parameter is defined as a shift of a continuous",
          "treatment. The treatment detected is NOT continuous."
        )
        stop(msg)
      }

      # produce trained likelihood when likelihood_def provided
      # NOTE: this is used to pass in likelihood values
      likelihood_def <- self$options$likelihood_override

      # use pre-constructed (external) likelihood values
      if (!is.null(likelihood_def)) {
        likelihood <- likelihood_def$train(tmle_task)
      } else {
        # covariates
        W_factor <- define_lf(LF_emp, "W")

        # treatment, with bounded likelihood away from 0
        A_bound <- c(1 / tmle_task$nrow, Inf)
        A_factor <- define_lf(LF_fit, "A",
          learner = learner_list[["A"]],
          bound = A_bound
        )

        # outcome
        Y_factor <- define_lf(LF_fit, "Y",
          learner = learner_list[["Y"]],
          type = "mean"
        )

        # construct likelihood from factor list
        factor_list <- list(W_factor, A_factor, Y_factor)
        likelihood_def <- Likelihood$new(factor_list)

        # fit/train likelihood
        likelihood <- likelihood_def$train(tmle_task)
      }
      return(likelihood)
    },
    make_params = function(tmle_task, likelihood) {
      # unwrap internalized arguments
      shift_fxn <- self$options$shift_fxn
      shift_fxn_inv <- self$options$shift_fxn_inv
      delta_shift <- self$options$delta_shift
      max_shifted_ratio <- self$options$max_shifted_ratio

      # treatment likelihood bound (away from 0 for continuous A)
      A_bound <- c(1 / tmle_task$nrow, Inf)

      # define shift intervention (additive only for now)
      intervention <- tmle3::define_lf(LF_shift,
        name = "A",
        original_lf = likelihood$factor_list[["A"]],
        likelihood_base = likelihood,                    # initial likelihood
        shift_fxn, shift_fxn_inv,                        # shift fxns
        shift_delta = delta_shift,                       # shift magnitude
        max_shifted_ratio = max_shifted_ratio,           # max ratio difference
        bound = A_bound                                  # bound shifted g
      )

      # create parameter (counterfactual treatment-specific mean)
      shifted_mean <- tmle3::Param_TSM$new(likelihood, intervention)

      # output should be a list
      tmle_params <- list(shifted_mean)
      return(tmle_params)
    },
    make_updater = function(...) {
      updater <- tmle3_Update$new(...)
    }
  ),
  active = list(),
  private = list()
)

################################################################################

#' Outcome under Shifted Treatment
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
#' @param shift_val A \code{numeric}, specification of the magnitude of the
#'  desired shift (on the level of the treatment). This is a value passed to
#'  the \code{function}s above for modulating the treatment.
#' @param max_shifted_ratio A \code{numeric} value indicating the maximum
#'  tolerance for the ratio of the counterfactual and observed intervention
#'  densities. In particular, the shifted value of the intervention is assigned
#'  to a given observational unit when the ratio of counterfactual intervention
#'  density to the observed intervention density is below this value.
##' @param ... Additional arguments (currently unused).
#'
#' @importFrom sl3 make_learner Lrnr_mean
#'
#' @export
#
tmle_shift <- function(shift_fxn = shift_additive_bounded,
                       shift_fxn_inv = shift_additive_bounded_inv,
                       shift_val = 1, max_shifted_ratio = 2, ...) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_shift$new(
    shift_fxn, shift_fxn_inv,
    shift_val, max_shifted_ratio,
    ...
  )
}
