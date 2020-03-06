context("Incorporating corrections for missingness in covariates")

library(data.table)
library(assertthat)
library(uuid)
library(sl3)

# setup data for test
set.seed(34831)
data(cpp)
data <- as.data.table(cpp)
data$parity01 <- as.numeric(data$parity > 0)
data$parity01_fac <- factor(data$parity01)
data$haz01 <- as.numeric(data$haz > 0)

node_list <- list(
  W = c(
    "apgar1", "apgar5", "gagebrth", "mage",
    "meducyrs", "sexn"
  ),
  A = "waz",
  Y = "haz01"
)

# drop missing A for now, might add back to test later
missing_W <- apply(is.na(data[, c(node_list$W,node_list$A),
                         with = FALSE]), 1, any)
data <- data[!missing_W]

# learners used for conditional expectation regression (e.g., outcome)
mean_lrnr <- Lrnr_mean$new()
glm_lrnr <- Lrnr_glm$new()
logit_metalearner <- make_learner(
  Lrnr_solnp, metalearner_logistic_binomial,
  loss_loglik_binomial
)
sl_lrnr <- Lrnr_sl$new(
  learners = list(mean_lrnr, glm_lrnr),
  metalearner = logit_metalearner 
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
learner_list <- list(Y = Q_learner, A = g_learner, delta_Y = Q_learner)

# initialize a tmle specification
tmle_spec <- tmle_shift(
  shift_val = 0.5,
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)

## define data (from tmle3_Spec base class)
tmle_task <- tmle_spec$make_tmle_task(data, node_list)
Q_task <- tmle_task$get_regression_task("Y", drop_censored = TRUE)
Q_learner <- learner_list$Y
Q_fit <- Q_learner$train(Q_task)

## define likelihood (from tmle3_Spec base class)
likelihood_init <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

## define update method (submodel and loss function)
updater <- tmle_spec$make_updater()
likelihood_targeted <- Targeted_Likelihood$new(likelihood_init, updater)

## define param
tmle_params <- tmle_spec$make_params(tmle_task, likelihood_targeted)
updater$tmle_params <- tmle_params

## fit tmle update
tmle_fit <- fit_tmle3(tmle_task, likelihood_targeted, tmle_params, updater)
tmle_fit
