#' Helper Functions for Additive Shifts of Continuous-Valued Interventions
#'
#' @param tmle_task A \code{tmle3_Task} object containing data and nodes, as
#'  described and implemented in the \code{tmle3} package. Please refer to the
#'  documentation and supporting materials of that package for details.
#' @param delta A \code{numeric} value giving a value of the shift to be applied
#'  the treatment. This is an additive shift of the treatment so it will merely
#'  be added to the observed value of the treatment node "A".
#'
#' @family shifting_interventions
#'
#' @rdname additive_shifting
#
shift_additive <- function(tmle_task, delta = 0.5) {
  out <- tmle_task$get_tmle_node("A") + delta
  return(out)
}

#' @param tmle_task A \code{tmle3_Task} object containing data and nodes, as
#'  described and implemented in the \code{tmle3} package. Please refer to the
#'  documentation and supporting materials of that package for details.
#' @param delta A \code{numeric} value giving a value of the shift to be applied
#'  the treatment. This is the inverse of an additive shift of the treatment so
#'  it will be subtracted from the observed value of the treatment node "A".
#'
#' @family shifting_interventions
#'
#' @rdname additive_shifting
#
shift_additive_inverse <- function(tmle_task, delta = 0.5) {
  out <- tmle_task$get_tmle_node("A") - delta
  return(out)
}

