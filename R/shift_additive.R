#' Additive Shifts of Continuous-Valued Interventions Without Bounds
#'
#' @param tmle_task A \code{tmle3_Task} object containing data and nodes, as
#'  described and implemented in the \code{tmle3} package. Please refer to the
#'  documentation and supporting materials of that package for details.
#' @param delta A \code{numeric} value giving a value of the shift to be
#'  applied to the treatment. This is an additive shift so the value is merely
#'  to be added to the observed value of the treatment node "A". In the case of
#'  the inverse additive shift, the specified value will be subtracted from the
#'  observed value of the treatment node "A".
#' @param ... Additional arguments (currently unused).
#'
#' @family shifting_interventions
#'
#' @rdname additive_shifting
#'
#' @export
shift_additive <- function(tmle_task, delta = 0, ...) {
  return(tmle_task$get_tmle_node("A") + delta)
}

#' @family shifting_interventions
#'
#' @rdname additive_shifting
#'
#' @export
shift_additive_inv <- function(tmle_task, delta = 0, ...) {
  return(tmle_task$get_tmle_node("A") - delta)
}
