# applications/multinomial_logistic_regression/shared_effect/scripts/generate_data.R

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
suppressMessages(i_am("applications/multinomial_logistic_regression/shared_effect/scripts/generate_data.R"))

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
# ✅ Load config
# -------------------------------
config_path <- here("config", "exps", paste0(exp_id, ".yml"))
if (!file.exists(config_path)) {
  stop("[ERROR] Config file not found at: ", config_path)
}
exp_config  <- read_yaml(config_path)

# -------------------------------
# ✅ Load helpers
# -------------------------------
estimand_helpers_dir <- here("applications", "multinomial_logistic_regression", "shared_effect", "scripts", "helpers")
common_helpers_dir  <- here("common", "scripts", "helpers")
miceadds::source.all(common_helpers_dir, print.source = FALSE)
miceadds::source.all(estimand_helpers_dir, print.source = FALSE)

# -------------------------------
# ✅ Setup directories
# -------------------------------
true_params_dir <- here("experiments", exp_id, "true_params")
iter_dir <- here("experiments", exp_id, "simulations", sim_id, iter_id)
data_dir <- here("experiments", exp_id, "simulations", sim_id, iter_id, "data")
dir_create(data_dir)
config_snapshot_path <- here(iter_dir, "config_snapshot.yml")

# -------------------------------
# ✅ Step 1: Load true parameters
# -------------------------------
Beta_0_path <- here(true_params_dir, "Beta_0.rds")
if (file_exists(Beta_0_path)) {
  message("[INFO] Loading Beta_0 from: ", Beta_0_path)
  Beta_0 <- readRDS(Beta_0_path)
} else {
  stop("[ERROR] Beta_0.rds not found at: ", Beta_0_path)
}

# -------------------------------
# ✅ Step 2: Generate new data
# -------------------------------
message("[INFO] Generating new data for iteration: ", iter_id)
seed <- get_seed_for_iter(exp_config$data_generation$seed, iter_id)
set.seed(seed)

data <- generate_data(exp_config, Beta_0)
save_list_objects(data, data_dir)

# -------------------------------
# ✅ Step 3: Save config snapshot
# -------------------------------
config_snapshot <- exp_config
config_snapshot$experiment$sim_id <- sim_id
config_snapshot$experiment$iter_id <- iter_id
write_strict_yaml(config_snapshot, config_snapshot_path)
message("[✓] Saved config snapshot to: ", config_snapshot_path)
