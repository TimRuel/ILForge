
exp_id <- "exp_v5.3.0"
sim_id <- "sim_01"
iter_id <- "iter_01"

args <- c(
  
  "poisson",
  "group_rates_weighted_sum",
  "fixed_effects_regression",
  exp_id,
  sim_id,
  iter_id
)

config <- yaml::read_yaml(paste0("C:/Northwestern/ILForge/config/exps/", exp_id, ".yml"))

config <- yaml::read_yaml(paste0("C:/Northwestern/ILForge/experiments/", exp_id, "/simulations/", sim_id, "/", iter_id, "/config_snapshot.yml"))

true_params_dir <- paste0("C:/Northwestern/ILForge/experiments/", exp_id, "/true_params")

data_dir <- paste0("C:/Northwestern/ILForge/experiments/", exp_id, "/simulations/", sim_id, "/", iter_id, "/data")

Beta_0 <- readRDS(file.path(true_params_dir, "Beta_0.rds"))

theta_0 <- readRDS(file.path(true_params_dir, "theta_0.rds"))

weights <- readRDS(file.path(true_params_dir, "weights.rds"))

data <- readRDS(file.path(data_dir, "data.rds"))

X <- readRDS(file.path(data_dir, "X.rds"))

results_dir <- paste0("C:/Northwestern/ILForge/experiments/", exp_id, "/simulations/", sim_id, "/", iter_id, "/results")

branch_params_list <- readRDS(file.path(results_dir, "branch_params_list.rds"))

integrated_LL <- readRDS(file.path(results_dir, "integrated_LL.rds"))

profile_LL <- readRDS(file.path(results_dir, "profile_LL.rds"))

MLE_data <- readRDS(file.path(results_dir, "MLE_data.rds"))

conf_ints <- readRDS(file.path(results_dir, "conf_ints.rds"))

fit <- fit_model(config, data)

data |> 
  mutate(mu_hat = predict(fit, type = "response"),
         theta_hat = mu_hat / t) |> 
  group_by(group) |> 
  summarise(theta_hat = mean(theta_hat)) |> 
  deframe()

X <- model.matrix(~ 0 + group + X1 + group:X2, data = data)
rownames(X) <- data$group

eta_hat <- mat %*% coef(model)

data |> 
  mutate(eta_hat = eta_hat) |> 
  group_by(group) |> 
  summarise(theta_hat = mean(exp(eta_hat))) |> 
  add_column(theta_0)


(X_mc %*% (coef(fit2) |> matrix()))[1:1e6,] |> exp() |> mean()

test

# Design matrix R used internally
X_fit <- model.matrix(fit)

# Manual reconstruction
eta_manual <- X_fit %*% coef(fit)

# Compare
all.equal(as.numeric(eta_manual),
          predict(fit, type="link"))  # should be TRUE

num_workers <- config$optimization_specs$IL$num_workers
plan(multisession, workers = num_workers, gc = TRUE)
branch_params_list <- get_branch_params_list(config, data, X, weights)
plan(sequential)

for (branch in integrated_LL$IL_branches) {
  
  plot(branch)
}

