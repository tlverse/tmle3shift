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
#' @param max_gn_ratio A \code{numeric} value indicating the maximum tolerance
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
                                   max_gn_ratio, ...) {
  # ratio of observed and shifted intervention densities
  gn_ratio <- get_density_ratio(tmle_task, delta, likelihood_base)

  # compute realistic value of intervention
  observed_a <- tmle_task$get_tmle_node("A")
  do_shift <- ifelse(gn_ratio < max_gn_ratio, delta, 0)
  shifted_a <- observed_a + do_shift
  return(shifted_a)
}

#' @family shifting_interventions
#'
#' @rdname additive_shifting_bounded
#'
#' @export
#
shift_additive_bounded_inv <- function(tmle_task, delta, likelihood_base,
                                       max_gn_ratio, ...) {
  # ratio of observed and shifted intervention densities
  gn_ratio <- get_density_ratio(tmle_task, delta, likelihood_base)

  # compute realistic value of intervention
  observed_a <- tmle_task$get_tmle_node("A")
  do_shift <- ifelse(gn_ratio < max_gn_ratio, delta, 0)
  shifted_a <- observed_a - do_shift
  return(shifted_a)
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
  cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(),
                                                    data.table(A = shifted_a))

  # find densities associated with tasks with observed and shifted intervention
  gn_a_obs <- likelihood_base$get_likelihoods(tmle_task, "A")
  gn_star_a_cf <- likelihood_base$get_likelihoods(cf_task, "A")
  gn_ratio <- gn_star_a_cf / gn_a_obs
  names(gn_ratio) <- NULL
  return(gn_ratio)
}

