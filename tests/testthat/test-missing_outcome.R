context("Incorporating corrections for missingness in covariates")
library(data.table)
library(assertthat)
library(uuid)
library(sl3)
library(tmle3)
set.seed(34831)

# setup data for test
data(cpp)
data <- as.data.table(cpp)
data[, parity01 := as.numeric(data$parity > 0)]
data[, parity01_fac := factor(data$parity01)]
data[, haz01 := as.numeric(data$haz > 0)]

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
xgb_lrnr <- Lrnr_xgboost$new()
logit_metalearner <- make_learner(
  Lrnr_solnp, metalearner_logistic_binomial,
  loss_loglik_binomial
)
sl_lrnr <- Lrnr_sl$new(
  learners = list(mean_lrnr, glm_lrnr, xgb_lrnr),
  metalearner = logit_metalearner
)

# learners used for conditional density estimation (i.e., propensity score)
hse_learner <- make_learner(Lrnr_density_semiparametric,
  mean_learner = glm_lrnr
)
mvd_learner <- make_learner(Lrnr_density_semiparametric,
  mean_learner = xgb_lrnr,
  var_learner = glm_lrnr
)
sl_density_lrnr <- Lrnr_sl$new(
  learners = Stack$new(hse_learner, mvd_learner),
  metalearner = Lrnr_solnp_density$new()
)

# specify outcome and treatment regressions and create learner list
Q_learner <- sl_lrnr
g_learner <- sl_density_lrnr
learner_list <- list(Y = Q_learner, A = g_learner, delta_Y = Q_learner)

# initialize a tmle specification
tmle_spec <- tmle_shift(
  shift_val = 0.5,
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)

# define data (from tmle3_Spec base class)
tmle_task <- tmle_spec$make_tmle_task(data, node_list)
Q_task <- tmle_task$get_regression_task("Y", drop_censored = TRUE)
Q_learner <- learner_list$Y
Q_fit <- Q_learner$train(Q_task)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task,
                                                        learner_list)

# define update method: submodel + loss function
updater <- tmle3_Update$new()
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)

# define param
tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params

# fit tmle update
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
