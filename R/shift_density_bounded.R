#' Additive Shifts of Continuous-Valued Interventions based on Bounded Densities
#'
#' @param tmle_task A \code{tmle3_Task} object containing data and nodes, as
#'  described and implemented in the \code{tmle3} package. Please refer to the
#'  documentation and supporting materials of that package for details.
#' @param delta A \code{numeric} value giving a value of the shift to be applied
#'  to the treatment. This is an additive shift so the value is merely to be
#'  added to the observed value of the treatment node "A". In the case of the
#'  inverse additive shift, the specified value will be subtracted from the
#'  observed value of the treatment node "A".
#' @param likelihood_base The base observed data likelihood, to be used in
#'  implementing guards that ensure that the shifted treatment does not violate
#'  the bounds induced by the support of the intervention, conditional on the
#'  covariates.
#' @param max_shifted_ratio A \code{numeric} value indicating maximum tolerance
#'  for the ratio of the counterfactual and observed intervention densities. In
#'  particular, the shifted value of the intervention is assigned to a given
#'  observational unit when the ratio of the counterfactual intervention density
#'  to the observed intervention density is below this value.
#' @param ... Additional arguments (currently unused).
#'
#' @family shifting_interventions
#'
#' @rdname additive_shifting_bounded
#'
#' @export
#
shift_additive_bounded <- function(tmle_task, delta, likelihood_base,
                                   max_shifted_ratio, ...) {
  # ratio of observed and shifted intervention densities
  intervention_density_ratio <- get_density_ratio(
    tmle_task, delta,
    likelihood_base
  )

  # compute realistic value of intervention
  observed_a <- tmle_task$get_tmle_node("A")
  shift_amt <- ifelse(intervention_density_ratio < max_shifted_ratio, delta, 0)
  shifted_a <- observed_a + shift_amt
  return(shifted_a)
}

#' @family shifting_interventions
#'
#' @rdname additive_shifting_bounded
#'
#' @export
#
shift_additive_bounded_inv <- function(tmle_task, delta, likelihood_base,
                                       max_shifted_ratio, ...) {
  # ratio of observed and shifted intervention densities
  intervention_density_ratio <- get_density_ratio(
    tmle_task, delta,
    likelihood_base
  )

  # compute realistic value of intervention
  observed_a <- tmle_task$get_tmle_node("A")
  shift_amt <- ifelse(intervention_density_ratio < max_shifted_ratio, delta, 0)
  shifted_inv_a <- observed_a - shift_amt
  return(shifted_inv_a)
}


#' Compute Ratio of Observed and Counterfactual Intervention Densities
#'
#' @param tmle_task A \code{tmle3_Task} object containing data and nodes, as
#'  described and implemented in the \code{tmle3} package. Please refer to the
#'  documentation and supporting materials of that package for details.
#' @param delta A \code{numeric} value giving a value of the shift to be applied
#'  to the treatment. This is an additive shift so the value is merely to be
#'  added to the observed value of the treatment node "A". In the case of the
#'  inverse additive shift, the specified value will be subtracted from the
#'  observed value of the treatment node "A".
#' @param likelihood_base The base observed data likelihood, to be used in
#'  implementing guards that ensure that the shifted treatment does not violate
#'  the bounds induced by the support of the intervention, conditional on the
#'  covariates.
#'
#' @importFrom data.table data.table
#'
#' @family shifting_interventions
#'
#' @rdname additive_shifting_bounded
#'
#' @keywords internal
#
get_density_ratio <- function(tmle_task, delta, likelihood_base) {
  # first, extract observed natural value of treatment and find shifted values
  obs_a <- tmle_task$get_tmle_node("A")
  shifted_a <- obs_a - delta

  # generate counterfactual task from shifted values
  cf_task <- tmle_task$generate_counterfactual_task(
    UUIDgenerate(),
    data.table::data.table(A = shifted_a)
  )

  # find densities associated with tasks with observed and shifted intervention
  emp_intervention_density <- likelihood_base$get_likelihoods(tmle_task, "A")
  cf_intervention_density <- likelihood_base$get_likelihoods(cf_task, "A")

  # compute ratio of counterfactual and empirical intervention densities
  intervention_density_ratio <-
    cf_intervention_density / emp_intervention_density
  names(intervention_density_ratio) <- NULL
  return(intervention_density_ratio)
}
