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
#' @param fold_number Whether to use cross-validated likelihood factor estimates
#'  or not. Passed through to method \code{get_likelihoods} in \pkg{tmle3}.
#' @param ... Additional arguments (currently unused).
#'
#' @family shifting_interventions
#'
#' @rdname additive_shifting_bounded
#'
#' @export
#
shift_additive_bounded <- function(tmle_task, delta, likelihood_base,
                                   max_shifted_ratio, fold_number, ...) {
  # ratio of observed and shifted intervention densities
  intervention_density_ratio <- get_density_ratio(
    tmle_task = tmle_task,
    delta = delta,
    likelihood_base = likelihood_base,
    fold_number = fold_number
  )

  # set NA values to 0 for density comparison
  na_mask_intervention_density_ratio <- is.na(intervention_density_ratio)
  if (FALSE %in% na_mask_intervention_density_ratio) {
    intervention_density_ratio[na_mask_intervention_density_ratio] <- 0
  }

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
                                       max_shifted_ratio, fold_number, ...) {
  # ratio of observed and shifted intervention densities
  intervention_density_ratio <- get_density_ratio(
    tmle_task = tmle_task,
    delta = delta,
    likelihood_base = likelihood_base,
    fold_number = fold_number
  )

  # set NA values to 0 for density comparison
  na_mask_intervention_density_ratio <- is.na(intervention_density_ratio)
  if (FALSE %in% na_mask_intervention_density_ratio) {
    intervention_density_ratio[na_mask_intervention_density_ratio] <- 0
  }

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
#' @param fold_number Whether to use cross-validated likelihood factor estimates
#'  or not. Passed directly to method \code{get_likelihoods} in \pkg{tmle3}.
#'
#' @importFrom data.table data.table
#'
#' @family shifting_interventions
#'
#' @rdname additive_shifting_bounded
#'
#' @keywords internal
#
get_density_ratio <- function(tmle_task, delta, likelihood_base, fold_number) {
  # extract observed natural value of treatment and compute shifted values
  shifted_a <- tmle_task$get_tmle_node("A") - delta

  # generate counterfactual task from shifted values
  cf_task <- tmle_task$generate_counterfactual_task(
    UUIDgenerate(),
    data.table::data.table(A = shifted_a)
  )

  # find densities associated with tasks with observed and shifted intervention
  emp_intervention_density <-
    likelihood_base$get_likelihoods(tmle_task = tmle_task, nodes = "A",
                                    fold_number = fold_number)
  cf_intervention_density <-
    likelihood_base$get_likelihoods(tmle_task = cf_task, nodes = "A",
                                    fold_number = fold_number)

  # compute ratio of counterfactual and empirical intervention densities
  intervention_density_ratio <-
    cf_intervention_density / emp_intervention_density
  names(intervention_density_ratio) <- NULL
  return(intervention_density_ratio)
}
