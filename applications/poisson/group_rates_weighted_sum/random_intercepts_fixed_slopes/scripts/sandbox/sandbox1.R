
exp_id <- "exp_v6.1.0"
sim_id <- "sim_01"
iter_id <- "iter_01"

args <- c(
  
  "poisson",
  "group_rates_weighted_sum",
  "random_intercepts_fixed_slopes",
  exp_id,
  sim_id,
  iter_id
)

config <- yaml::read_yaml(paste0("C:/Northwestern/ILForge/config/exps/", exp_id, ".yml"))

config <- yaml::read_yaml(paste0("C:/Northwestern/ILForge/experiments/", exp_id, "/simulations/", sim_id, "/", iter_id, "/config_snapshot.yml"))

true_params_dir <- paste0("C:/Northwestern/ILForge/experiments/", exp_id, "/true_params")

data_dir <- paste0("C:/Northwestern/ILForge/experiments/", exp_id, "/simulations/", sim_id, "/", iter_id, "/data")

Beta_0 <- readRDS(file.path(true_params_dir, "Beta_0.rds"))

sigma_alpha <- readRDS(file.path(true_params_dir, "sigma_alpha.rds"))

theta_0 <- readRDS(file.path(true_params_dir, "theta_0.rds"))

weights <- readRDS(file.path(true_params_dir, "weights.rds"))

data <- readRDS(file.path(data_dir, "data.rds"))

X <- readRDS(file.path(data_dir, "X.rds"))

intercepts <- readRDS(file.path(data_dir, "intercepts.rds"))

results_dir <- paste0("C:/Northwestern/ILForge/experiments/", exp_id, "/simulations/", sim_id, "/", iter_id, "/results")

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


