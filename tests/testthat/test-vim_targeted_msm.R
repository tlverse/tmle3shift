context("Variable importance measures directly targeting MSM parameters")

library(data.table)
library(sl3)
library(tmle3)
set.seed(429153)

## simulate simple data for tmle-shift sketch
n_obs <- 1000 # number of observations
n_w <- 1 # number of baseline covariates
tx_mult <- 2 # multiplier for effect of W = 1 on treatment
delta_grid <- seq(-1, 1, 1) # grid of shifts to consider

## baseline covariates -- simple, binary
W <- as.numeric(replicate(n_w, rbinom(n_obs, 1, 0.5)))

## create treatment based on baseline W
A <- as.numeric(rnorm(n_obs, mean = tx_mult * W, sd = 1))

## create outcome as a linear function of A, W + white noise
Y <- A + W + rnorm(n_obs, mean = 0, sd = 1)

## organize data and nodes for tmle3
data <- data.table(W, A, Y)
node_list <- list(W = "W", A = "A", Y = "Y")


# learners used for conditional expectation regression (e.g., outcome)
mean_lrnr <- Lrnr_mean$new()
glm_lrnr <- Lrnr_glm$new()
sl_lrnr <- Lrnr_sl$new(
  learners = list(mean_lrnr, glm_lrnr),
  metalearner = Lrnr_nnls$new()
)

# learners used for conditional density regression (i.e., propensity score)
haldensify_lrnr <- Lrnr_haldensify$new(
  n_bins = 5, grid_type = "equal_mass",
  lambda_seq = exp(seq(-1, -13, length = 100))
)
cv_haldensify_lrnr <- Lrnr_cv$new(haldensify_lrnr, full_fit = TRUE)

# specify outcome and treatment regressions and create learner list
Q_learner <- sl_lrnr
g_learner <- cv_haldensify_lrnr
learner_list <- list(Y = Q_learner, A = g_learner)


# initialize a tmle specification for direct targeting of MSM parameters
tmle_spec <- tmle_vimshift_msm(
  shift_grid = delta_grid,
  max_shifted_ratio = 2
)

## define data (from tmle3_Spec base class)
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

## define likelihood (from tmle3_Spec base class)
likelihood_init <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

## define update method (fluctuation submodel and loss function)
updater <- tmle_spec$make_updater()
likelihood_targeted <- Targeted_Likelihood$new(likelihood_init, updater)

## invoke params specified in spec
tmle_params <- tmle_spec$make_params(tmle_task, likelihood_targeted)
updater$tmle_params <- tmle_params

## fit TML estimator update
tmle_fit_targeted_msm <- fit_tmle3(
  tmle_task, likelihood_targeted, tmle_params,
  updater
)


# initialize a tmle specification for delta method
tmle_spec <- tmle_vimshift_delta(
  shift_grid = delta_grid,
  max_shifted_ratio = 2
)

## define data (from tmle3_Spec base class)
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

## define likelihood (from tmle3_Spec base class)
likelihood_init <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

## define update method (fluctuation submodel and loss function)
updater <- tmle_spec$make_updater()
likelihood_targeted <- Targeted_Likelihood$new(likelihood_init, updater)

## invoke params specified in spec
tmle_params <- tmle_spec$make_params(tmle_task, likelihood_targeted)
updater$tmle_params <- tmle_params

## fit TML estimator update
tmle_fit_delta_method <- fit_tmle3(
  tmle_task, likelihood_targeted, tmle_params,
  updater
)


## extract relevant tmle3 results for test and re-format appropriately
msm_delta_summary <- tmle_fit_delta_method$summary[4:5, ]
msm_targeted_summary <- tmle_fit_targeted_msm$summary

test_that("Targeted MSM approach and delta method MSM match nearly", {
  expect_equal(msm_delta_summary, msm_targeted_summary,
    tol = 0.01,
  )
})
