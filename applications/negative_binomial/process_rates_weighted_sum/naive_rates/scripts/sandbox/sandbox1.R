# Variables
experiment_name <- "negative_binomial"
estimand_name <- "process_rates_weighted_sum"
model_name <- "naive_rates"
exp_ver <- "exp_v7.0.2"
sim_num <- "sim_01"
iter_num <- "iter_01"

args <- c(
  experiment_name,
  estimand_name,
  model_name,
  exp_ver,
  sim_num,
  iter_num
)

proj_root <- "C:/Northwestern/ILForge"

exp_path <- file.path(proj_root, "experiments", exp_ver)
true_param_path <- file.path(exp_path, "true_params")
sim_iter_path <- file.path(exp_path, "simulations", sim_num, iter_num)

config <- yaml::read_yaml(file.path(sim_iter_path, "config_snapshot.yml"))

phi_0          <- readRDS(file.path(true_param_path, "phi_0.rds"))
theta_0        <- readRDS(file.path(true_param_path, "theta_0.rds"))
weights        <- readRDS(file.path(true_param_path, "weights.rds"))
n_per_process  <- readRDS(file.path(true_param_path, "n_per_process.rds"))
data           <- readRDS(file.path(sim_iter_path, "data", "data.rds"))


