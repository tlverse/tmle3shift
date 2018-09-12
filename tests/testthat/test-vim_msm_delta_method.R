context("Variable importance for shift interventions via delta method")

library(uuid)
library(assertthat)
library(data.table)
library(sl3)
library(tmle3)
library(ranger)
set.seed(429153)

################################################################################
# setup data and learners for tests
################################################################################

## simulate simple data for tmle-shift sketch
n_obs <- 1000 # number of observations
n_w <- 1 # number of baseline covariates
tx_mult <- 2 # multiplier for the effect of W = 1 on the treatment

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
lrn1 <- Lrnr_mean$new()
lrn2 <- Lrnr_glm$new()
lrn3 <- Lrnr_ranger$new()
sl_lrn <- Lrnr_sl$new(
  learners = list(lrn1, lrn2, lrn3),
  metalearner = Lrnr_nnls$new()
)

# learners used for conditional density regression (e.g., propensity score)
lrn1_dens <- Lrnr_condensier$new(
  nbins = 25, bin_estimator = lrn1,
  bin_method = "dhist"
)
lrn2_dens <- Lrnr_condensier$new(
  nbins = 20, bin_estimator = lrn2,
  bin_method = "dhist"
)
lrn3_dens <- Lrnr_condensier$new(
  nbins = 15, bin_estimator = lrn3,
  bin_method = "dhist"
)
sl_lrn_dens <- Lrnr_sl$new(
  learners = list(lrn1_dens, lrn2_dens, lrn3_dens),
  metalearner = Lrnr_solnp_density$new()
)

# specify outcome and treatment regressions and create learner list
Q_learner <- sl_lrn
g_learner <- sl_lrn_dens
learner_list <- list(Y = Q_learner, A = g_learner)


################################################################################
# CONSTRUCT MSM ESTIMATES FROM INDIVIDUAL TML ESTIMATES
################################################################################

# what's the grid of shifts we wish to consider?
delta_grid <- seq(-1, 1, 0.5)

# initialize a tmle specification
tmle_spec <- tmle_vimshift_delta(shift_grid = delta_grid,
                                 max_shifted_ratio = 2)

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
tmle_fit <- fit_tmle3(tmle_task, likelihood_targeted, tmle_params, updater)

## extract results from tmle3_Fit object
tmle_fit

## use MSM to summarize results
msm <- trend_msm(tmle_fit, delta_grid)
msm

