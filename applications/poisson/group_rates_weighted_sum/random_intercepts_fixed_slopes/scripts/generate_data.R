# applications/poisson/group_rates_weighted_sum/random_intercepts_fixed_slopes/scripts/generate_data.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(nloptr)
  library(here)
  library(yaml)
  library(fs)
})

# -------------------------------
# ✅ Anchor project root
# -------------------------------
suppressMessages(i_am("applications/poisson/group_rates_weighted_sum/random_intercepts_fixed_slopes/scripts/generate_data.R"))

# -------------------------------
# ✅ Parse arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3) {
  exp_id <- args[[1]]
  sim_id <- args[[2]]
  iter_id <- args[[3]]
} else {
  stop("Usage: Rscript generate_data.R <exp_id> <sim_id> <iter_id>")
}

# -------------------------------
# ✅ Load helpers
# -------------------------------
model_helpers_dir <- here("applications", "poisson", "group_rates_weighted_sum", "random_intercepts_fixed_slopes", "scripts", "helpers")
common_helpers_dir  <- here("common", "scripts", "helpers")
miceadds::source.all(common_helpers_dir, print.source = FALSE)
miceadds::source.all(model_helpers_dir, print.source = FALSE)

# -------------------------------
# ✅ Setup directories
# -------------------------------
true_params_dir <- here("experiments", exp_id, "true_params")
iter_dir <- here("experiments", exp_id, "simulations", sim_id, iter_id)
data_dir <- here(iter_dir, "data")
dir_create(data_dir)
config_snapshot_path <- here(iter_dir, "config_snapshot.yml")

# -------------------------------
# ✅ Load experiment config
# -------------------------------
exp_config_path <- here("config", "exps", paste0(exp_id, ".yml"))
if (!file.exists(exp_config_path)) {
  stop("❌ Experiment config file not found at: ", exp_config_path)
}
exp_config  <- read_yaml(exp_config_path)

# -------------------------------
# ✅ Save config snapshot
# -------------------------------
config_snapshot <- exp_config
config_snapshot$experiment$sim_id <- sim_id
config_snapshot$experiment$iter_id <- iter_id
write_strict_yaml(config_snapshot, config_snapshot_path)
message("[✓] Saved config snapshot to: ", config_snapshot_path)

# -------------------------------
# ✅ Generate new data
# -------------------------------
message("[INFO] Generating new data for iteration: ", iter_id)
iter_seed <- get_seed_for_iter(config_snapshot$model$seed, config_snapshot$experiment$iter_id)
set.seed(iter_seed)
Beta_0 <- readRDS(here(true_params_dir, "Beta_0.rds"))
sigma_alpha <- readRDS(here(true_params_dir, "sigma_alpha.rds"))
data <- generate_data(config_snapshot, Beta_0, sigma_alpha)

# Save the generated data
save_list_objects(data, data_dir)
message("[✓] Saved simulated data to: ", data_dir)
message("[INFO] Saved objects: ", paste(names(data), collapse = ", "))

