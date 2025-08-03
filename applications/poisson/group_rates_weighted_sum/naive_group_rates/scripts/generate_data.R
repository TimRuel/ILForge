# applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/generate_data.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
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
data_dir <- here("experiments", exp_id, "simulations", sim_id, iter_id, "data")
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
# ✅ Load true parameters
# -------------------------------
lambda_0_path <- here(true_params_dir, "lambda_0.rds")

if (!file_exists(lambda_0_path)) {
  stop("[ERROR] lambda_0.rds not found at: ", lambda_0_path)
}

message("[INFO] Loading lambda_0 from: ", lambda_0_path)
lambda_0 <- readRDS(lambda_0_path)

# -------------------------------
# ✅ Generate new data
# -------------------------------
message("[INFO] Generating new data for iteration: ", iter_id)
data <- generate_data(config_snapshot, lambda_0)
saveRDS(data, here(data_dir, "data.rds"))
