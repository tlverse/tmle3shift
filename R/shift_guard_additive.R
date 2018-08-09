#' Additive Shifts of Continuous-Valued Interventions with Bound Guards
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
#' @rdname additive_shifting_guard
#'
#' @export
#
shift_additive_guard <- function(tmle_task, delta = 0, likelihood_base) {
  # compute upper and lower bounds of treatment from observed data likelihood
  levels_a <- as.list(unique(sort(tmle_task$get_tmle_node("A"))))
  gn_a <- sapply(levels_a, function(a) {
    cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(),
                                                      data.table(A = a))
    dens_this_a <- likelihood_base$get_likelihoods(cf_task, "A")
  })
  gn_a_bool <- gn_a > 1 / tmle_task$nrow
  get_uw <- apply(gn_a_bool, 1, function(x) {
    max(which(x))
  })
  names(get_uw) <- NULL
  get_lw <- apply(gn_a_bool, 1, function(x) {
    min(which(x))
  })
  names(get_lw) <- NULL
  uw <- unlist(levels_a[get_uw])
  lw <- unlist(levels_a[get_lw])
  bounds_a <- cbind(lw, uw)

  # compute realistic value of intervention
  shifted_a <- tmle_task$get_tmle_node("A") + delta
  if (delta > 0) {
    # compare against upper bound (uw) for upward shift of treatment
    shifted_a_or_uw <- cbind(shifted_a, bounds_a[, 2])
    realistic_shifted_a <- apply(shifted_a_or_uw, 1, min)
  } else if (delta < 0) {
    # compare against lower bound (lw) for downward shift of treatment
    shifted_a_or_lw <- cbind(shifted_a, bounds_a[, 1])
    realistic_shifted_a <- apply(shifted_a_or_lw, 1, max)
  }
  return(realistic_shifted_a)
}

#' @family shifting_interventions
#'
#' @rdname additive_shifting_guard
#'
#' @export
#
shift_additive_guard_inv <- function(tmle_task, delta = 0, likelihood_base) {
  # compute upper and lower bounds of treatment from observed data likelihood
  levels_a <- as.list(unique(sort(tmle_task$get_tmle_node("A"))))
  gn_a <- sapply(levels_a, function(a) {
    cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(),
                                                      data.table(A = a))
    dens_this_a <- likelihood_base$get_likelihoods(cf_task, "A")
  })
  gn_a_bool <- gn_a > 1 / tmle_task$nrow
  get_uw <- apply(gn_a_bool, 1, function(x) {
    max(which(x))
  })
  names(get_uw) <- NULL
  get_lw <- apply(gn_a_bool, 1, function(x) {
    min(which(x))
  })
  names(get_lw) <- NULL
  uw <- unlist(levels_a[get_uw])
  lw <- unlist(levels_a[get_lw])
  bounds_a <- cbind(lw, uw)

  # compute realistic value of intervention
  shifted_a <- tmle_task$get_tmle_node("A") - delta
  # CAREFUL: this is an INVERSE shift function so the logic is backwards
  if (delta > 0) {
    # compare against upper bound (uw) for upward shift of treatment
    shifted_a_or_lw <- cbind(shifted_a, bounds_a[, 1])
    realistic_shifted_a <- apply(shifted_a_or_lw, 1, max)
  } else if (delta < 0) {
    # compare against lower bound (lw) for downward shift of treatment
    shifted_a_or_uw <- cbind(shifted_a, bounds_a[, 2])
    realistic_shifted_a <- apply(shifted_a_or_uw, 1, min)
  }
  return(realistic_shifted_a)
}

