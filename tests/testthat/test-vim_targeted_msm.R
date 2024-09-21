context("Variable importance measures directly targeting MSM parameters")

library(data.table)
library(sl3)
library(tmle3)
set.seed(429153)

## simulate simple data for tmle-shift sketch
n_obs <- 1000 # number of observations
n_w <- 2 # number of baseline covariates
tx_mult <- 2 # multiplier for effect of W = 1 on treatment
delta_grid <- seq(-1, 1, 1) # grid of shifts to consider

## baseline covariates -- simple, binary
W <- replicate(n_w, rbinom(n_obs, 1, 0.65))

## create treatment based on baseline W
A <- as.numeric(rnorm(n_obs, mean = tx_mult * rowSums(W), sd = 0.5))

## create outcome as a linear function of A, W + white noise
Y <- A + rowSums(W) + rnorm(n_obs, mean = 0, sd = 1)

## organize data and nodes for tmle3
data <- data.table(W, A, Y)
setnames(data, c("W1", "W2", "A", "Y"))
node_list <- list(W = c("W1", "W2"), A = "A", Y = "Y")

# learners used for outcome regression (conditional expectation)
mean_learner <- Lrnr_mean$new()
fglm_learner <- Lrnr_glm_fast$new()
xgb_learner <- Lrnr_xgboost$new()
sl_or_learner <- Lrnr_sl$new(
  learners = list(mean_learner, fglm_learner, xgb_learner)
)

# learners used for generalized propensity score (conditional density)
hse_mean_learner <- make_learner(Lrnr_density_semiparametric,
  mean_learner = mean_learner
)
hse_fglm_learner <- make_learner(Lrnr_density_semiparametric,
  mean_learner = fglm_learner
)
hse_xgb_learner <- make_learner(Lrnr_density_semiparametric,
  mean_learner = xgb_learner
)
mvd_xgb_fglm_learner <- make_learner(Lrnr_density_semiparametric,
  mean_learner = xgb_learner,
  var_learner = fglm_learner
)
sl_gps_learner <- Lrnr_sl$new(
  learners = Stack$new(hse_mean_learner, hse_xgb_learner, hse_fglm_learner,
                       mvd_xgb_fglm_learner),
  metalearner = Lrnr_solnp_density$new()
)

# specify outcome and treatment regressions and create learner list
or_learner <- sl_or_learner
gps_learner <- sl_gps_learner
learner_list <- list(Y = or_learner, A = gps_learner, delta_Y = or_learner)

# initialize a tmle specification for direct targeting of MSM parameters
tmle_spec <- tmle_vimshift_msm(
  shift_grid = delta_grid,
  max_shifted_ratio = 10
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
  max_shifted_ratio = 10
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
  expect_equal(msm_delta_summary$tmle_est, msm_targeted_summary$tmle_est,
               tol = 0.05)
})
