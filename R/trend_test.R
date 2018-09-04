# get EIF for each parameter after calling tmle3_fit


#' Test for a trend in effect of shifted intervention
#'
#' @param object An object of class \code{survtmle}.
#' @param level The nominal coverage probability of the confidence interval.
#'
#' @return An object of class \code{"trend_test"}.
#' \describe{
#'   \item{\code{alpha}}{The intercept from the projection onto the working
#'         model.}
#'   \item{\code{beta}}{The slope from the projection onto the working model,
#'         i.e., the "trend" parameter.}
#'   \item{\code{L_j}}{A data.frame showing the log ratio of cumulative
#'         incidences that were projected onto the working model.}
#'   \item{\code{se_beta}}{The influence function-based estimate of the standard
#'         error of \code{beta_n}.}
#'   \item{\code{ci}}{The confidence interval at the requested level for
#'         \code{beta_n}.}
#'   \item{\code{pval}}{The two-sided p-value from the test of the null
#'         hypothesis that \code{beta_n} = 0.}
#'   \item{\code{level}}{The requested level of the confidence interval.}
#' }
#'
#' @export
#
trend_test <- function(object, level = 0.95) {
  # this trend test only works for tmle3_fit type objects
  stopifnot(is(object, "tmle3_Fit"))

  # make sure more than one parameter has been estimated for trend
  stopifnot(length(object$estimates) > 1)

  # matrix of EIF(O_i) values and estimates across each parameter estimated
  D <- sapply(object$estimates, `[[`, "IC")
  F <- sapply(object$estimates, `[[`, "psi")
  n <- nrow(D)
  j_vec <- unique(as.numeric(unlist(lapply(strsplit(row.names(F), " "),"[[", 2))))

  # log ratios
  R_n <- g(F)

  # effect estimates
  est <- h(F, D, j_vec)
  alpha_n <- est[1]
  beta_n <- est[2]

  # standard error of effect estimate
  nabla_h <- grad_h(F, D, j_vec)
  se_beta_n <- sqrt(t(nabla_h) %*% cov(D) %*% nabla_h / n)

  # confidence interval
  ci <- beta_n + c(1,-1)*qnorm((1 - level)/2) * rep(se_beta_n, 2)

  # hypothesis test
  pval <- 2 * pnorm(-abs(beta_n / se_beta_n))

  # output
  out <- list(alpha = alpha_n, beta = beta_n, 
              L_j = data.frame(j = 1:length(R_n), L_jn = R_n), 
              se_beta = se_beta_n, ci = ci, pval = pval,
              level = level)
  class(out) <- "trend_test"
  return(out)
}

#' Print the output of a test for trends
#'
#' @param x An object of class \code{"trend_test"}.
#' @param digits Number of digits to round output.
#' @param ... Other options (not currently used)
#
#' @export
#
print.trend_test <- function(x, digits = 3, ...) {
  tmp <- data.frame(beta = round(x$beta, digits),
                    ci1 = round(x$ci[1], digits), 
                    ci2 = round(x$ci[2], digits), 
                    pval = round(x$pval, digits))
  colnames(tmp)[2:3] <- paste0(c("lower","upper"),"_",round(100*x$level), "%CI")
  cat("Trend in efficacy across failure type levels: \n")
  print(tmp)
}

#' Helper function to compute estimate of trend parameter
#'
#' @param F A vector of cumulative incidence estimates with alternating
#'  treatment = 0 and treatment = 1 for each type.
#' @param D A matrix of variance estimates corresponding with the vector
#'  \code{F}.
#' @param j_vec Vector of genetic distances
#
h <- function(F, D, j_vec) {
  R_n <- g(F)
  K <- length(F)/2
  nabla_g <- grad_g(F)
  Upsilon_n <- nabla_g %*% cov(D) %*% t(nabla_g)
  # design matrix
  X <- cbind(rep(1,length(j_vec)), j_vec)
  Upsilon_inv <- solve(Upsilon_n)
  S_n <- solve(t(X) %*% Upsilon_inv %*% X) %*% t(X) %*% Upsilon_inv
  S_n %*% R_n
}

#' Helper function to compute gradient of trend parameter function
#'
#' @param F A vector of cumulative incidence estimates with alternating
#'  treatment = 0 and treatment = 1 for each type.
#' @param D A matrix of variance estimates corresponding with the vector
#'  \code{F}.
#
grad_h <- function(F, D, j_vec) {
  R_n <- g(F)
  K <- length(F)/2
  nabla_g <- grad_g(F)
  Upsilon_n <- nabla_g %*% cov(D) %*% t(nabla_g) 
  # design matrix
  X <- cbind(rep(1,length(j_vec)), j_vec)
  Upsilon_inv <- solve(Upsilon_n)
  S_n <- solve(t(X) %*% Upsilon_inv %*% X) %*% t(X) %*% Upsilon_inv
  tmp <- (-1)^(2:(length(F)+1)) * (1 / F)
  tmp2 <- S_n[2, sort(rep(1:K, 2))]
  matrix(tmp*tmp2, ncol = 1)
} 

#' Helper function to got log ratio of cumulative incidences
#'
#' @param F A vector of cumulative incidence estimates with alternating
#'  treatment = 0 and treatment = 1 for each type.
#
g <- function(F) {
  matrix(log(F[seq(1, length(F), 2)] / F[seq(2, length(F), 2)]), ncol = 1)
}

#' Helper function to compute the gradient for the log ratio
#' transformation of the cumulative incidence estimates.
#'
#' @param F A vector of cumulative incidence estimates with alternating
#'  treatment = 0 and treatment = 1 for each type.
#
grad_g <- function(F) {
  K <- length(F) / 2
  tmp_vec <- rep(0, 2*K^2)
  ind <- as.numeric(paste0(sort(rep((1:K) - 1, 2)), 1:length(F)))
  ind[length(ind)] <- 2*K^2
  tmp_vec[ind] <- (-1)^(2:(length(F)+1)) * (1 / F)
  tmp <- matrix(tmp_vec, nrow = K, byrow = TRUE)
  tmp
}

