# applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/generate_data.R

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
suppressMessages(i_am("applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/generate_data.R"))

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
model_helpers_dir <- here("applications", "poisson", "group_rates_weighted_sum", "naive_group_rates", "scripts", "helpers")
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
# ✅ Load true parameters (all three)
# -------------------------------
theta_0_path     <- here(true_params_dir, "theta_0.rds")
weights_path     <- here(true_params_dir, "weights.rds")
n_per_group_path <- here(true_params_dir, "n_per_group.rds")

if (!file_exists(theta_0_path))     stop("[ERROR] theta_0.rds not found at: ", theta_0_path)
if (!file_exists(weights_path))     stop("[ERROR] weights.rds not found at: ", weights_path)
if (!file_exists(n_per_group_path)) stop("[ERROR] n_per_group.rds not found at: ", n_per_group_path)

message("[INFO] Loading theta_0 from: ", theta_0_path)
message("[INFO] Loading weights from: ", weights_path)
message("[INFO] Loading n_per_group from: ", n_per_group_path)

theta_0     <- readRDS(theta_0_path)
weights     <- readRDS(weights_path)
n_per_group <- readRDS(n_per_group_path)

# -------------------------------
# ✅ Generate new data
# -------------------------------
message("[INFO] Generating new data for iteration: ", iter_id)

# Pass the whole config snapshot AND the true_params_dir so generate_data() can read what it needs
data <- generate_data(config_snapshot, true_params_dir)

# Save the generated data
saveRDS(data, here(data_dir, "data.rds"))
message("[✓] Saved simulated data to: ", here(data_dir, "data.rds"))
