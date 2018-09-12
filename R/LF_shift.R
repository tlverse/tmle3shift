#' Shifted Likelihood Factor
#'
#' Shifts a likelihood factor according to a \code{shift_function} and a given
#' magnitude of the desired shift (\code{shift_delta}). In effect,
#' \code{get_likelihood(tmle_task)} from \code{tmle3} will instead be the
#' likelihood from the \code{original_lf}, but a for shifted value
#' \eqn{A'=}\code{shift_function}\eqn{(A, W)}
#'
#' @references
#' \describe{
#'   \item{"Stochastic Treatment Regimes."}{Díaz, Iván and van der Laan, Mark J
#'         (2018). In Targeted Learning in Data Science: Causal Inference for
#'         Complex Longitudinal Studies, 167–80. Springer Science & Business
#'         Media.}
#'   \item{"Population Intervention Causal Effects Based on Stochastic
#'         Interventions."}{Díaz, Iván and van der Laan, Mark J (2012).
#'         Biometrics 68 (2). Wiley Online Library: 541–49.}
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#'
#' @family Likelihood objects
#'
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_shift, name, type = "density", original_lf,
#'     shift_function, ...)}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node
#'           name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{original_lf}}{\code{\link{LF_base}} object, the likelihood
#'           factor to shift
#'     }
#'     \item{\code{shift_function}}{\code{function}, defines the shift
#'     }
#'     \item{\code{shift_inverse}}{\code{function}, the inverse of a given
#'           \code{shift_function}
#'     }
#'     \item{\code{shift_delta}}{\code{numeric}, specification of the magnitude
#'           of the desired shift (on the level of the treatment)
#'     }
#'     \item{\code{max_shifted_ratio}}{A \code{numeric} value indicating the maximum
#'           tolerance for the ratio of the counterfactual and observed
#'           intervention densities. In particular, the shifted value of the
#'           intervention is assigned to a given observational unit when the
#'           ratio of the counterfactual intervention density to the observed
#'           intervention density is below this value
##'    }
#'     \item{\code{...}}{Not currently used.
#'     }
#'   }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{original_lf}}{\code{\link{LF_base}} object, the likelihood
#'           factor to shift
#'     }
#'     \item{\code{shift_function}}{\code{function}, defines the shift
#'     }
#'     \item{\code{shift_inverse}}{\code{function}, the inverse of a given
#'           \code{shift_function}
#'     }
#'     \item{\code{shift_delta}}{\code{numeric}, specification of the magnitude
#'           of the desired shift (on the level of the treatment)
##'    }
#'     \item{\code{max_shifted_ratio}}{A \code{numeric} value indicating the maximum
#'           tolerance for the ratio of the counterfactual and observed
#'           intervention densities. In particular, the shifted value of the
#'           intervention is assigned to a given observational unit when the
#'           ratio of the counterfactual intervention density to the observed
#'           intervention density is below this value
##'    }
#'     \item{\code{...}}{Not currently used.
#'     }
#'   }
#'
#' @export
#
LF_shift <- R6::R6Class(
  classname = "LF_shift",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3::LF_base,
  public = list(
    initialize = function(name, original_lf, likelihood_base,
                              shift_function, shift_inverse, shift_delta,
                              max_shifted_ratio, ...) {
      super$initialize(name, ..., type = "density")
      private$.original_lf <- original_lf
      private$.likelihood_base <- likelihood_base
      private$.shift_function <- shift_function
      private$.shift_inverse <- shift_inverse
      private$.shift_delta <- shift_delta
      private$.max_shifted_ratio <- max_shifted_ratio
    },
    get_mean = function(tmle_task, cv_fold) {
      stop("get_mean not supported for LF_shift")
    },
    get_density = function(tmle_task, cv_fold) {
      # get shifted data
      shifted_values <- self$shift_inverse(
        tmle_task = tmle_task,
        delta = self$shift_delta,
        likelihood_base =
          self$likelihood_base,
        max_shifted_ratio =
          self$max_shifted_ratio
      )

      # generate cf_task data
      cf_data <- data.table(shifted_values)
      setnames(cf_data, self$name)

      cf_task <- tmle_task$generate_counterfactual_task(
        UUIDgenerate(),
        cf_data
      )

      # get original likelihood for shifted data
      cf_likelihood <- self$original_lf$get_likelihood(cf_task)
      return(cf_likelihood)
    },
    cf_values = function(tmle_task) {
      cf_values <- self$shift_function(
        tmle_task = tmle_task,
        delta = self$shift_delta,
        likelihood_base =
          self$likelihood_base,
        max_shifted_ratio =
          self$max_shifted_ratio
      )
      return(cf_values)
    }
  ),
  active = list(
    original_lf = function() {
      return(private$.original_lf)
    },
    likelihood_base = function() {
      return(private$.likelihood_base)
    },
    shift_function = function() {
      return(private$.shift_function)
    },
    shift_inverse = function() {
      return(private$.shift_inverse)
    },
    shift_delta = function() {
      return(private$.shift_delta)
    },
    max_shifted_ratio = function() {
      return(private$.max_shifted_ratio)
    }
  ),
  private = list(
    .name = NULL,
    .original_lf = NULL,
    .likelihood_base = NULL,
    .shift_function = NULL,
    .shift_inverse = NULL,
    .shift_delta = NULL,
    .max_shifted_ratio = NULL
  )
)
