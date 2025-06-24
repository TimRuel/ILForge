# applications/multinomial_logistic_regression/entropy/scripts/generate_data.R

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
i_am("applications/multinomial_logistic_regression/entropy/scripts/generate_data.R")

# -------------------------------
# ✅ Parse arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  exp_id <- args[[1]]
  sim_id <- NULL
  run_id <- args[[2]]
} else if (length(args) == 3) {
  exp_id <- args[[1]]
  sim_id <- args[[2]]
  run_id <- args[[3]]
} else {
  stop("Usage: Rscript generate_data.R <exp_id> [sim_id] <run_id>")
}

# -------------------------------
# ✅ Load config
# -------------------------------
config_path <- here("config", "exps", paste0(exp_id, ".yml"))
if (!file.exists(config_path)) {
  stop("[ERROR] Config file not found at: ", config_path)
}
exp_config  <- read_yaml(config_path)
X1_levels   <- exp_config$X1_levels
model_specs <- exp_config$model_specs
opt_specs   <- exp_config$optimization_specs
app_name    <- exp_config$experiment$app
estimand    <- exp_config$experiment$estimand

# -------------------------------
# ✅ Load helpers
# -------------------------------
helper_dir <- here("applications", app_name, estimand, "scripts", "helpers")
miceadds::source.all(helper_dir, print.source = FALSE)

# -------------------------------
# ✅ Setup directories
# -------------------------------
true_params_dir <- here("experiments", exp_id, "true_params")
run_dir <- if (is.null(sim_id)) {
  here("experiments", exp_id, "individual_runs", run_id)
} else {
  here("experiments", exp_id, "simulations", sim_id, run_id)
}
data_dir <- here(run_dir, "data")
dir_create(data_dir)

config_snapshot_path <- path(run_dir, "config_snapshot.yml")

# -------------------------------
# ✅ Step 1: Load true parameters
# -------------------------------
Beta_0_path <- path(true_params_dir, "Beta_0.rds")
if (file_exists(Beta_0_path)) {
  message("[INFO] Loading Beta_0 from: ", Beta_0_path)
  Beta_0 <- readRDS(Beta_0_path)
} else {
  stop("[ERROR] Beta_0.rds not found at: ", Beta_0_path)
}

# -------------------------------
# ✅ Step 2: Generate new data
# -------------------------------
message("[INFO] Generating new data for run: ", run_id)
seed <- get_seed_for_run(opt_specs$seed, run_id)
set.seed(seed)

mm_formula <- substring(model_specs$formula, 2)
data <- get_data(X1_levels, mm_formula, Beta_0)

save_list_objects(data, data_dir)

# -------------------------------
# ✅ Step 3: Save config snapshot
# -------------------------------
config_snapshot <- exp_config
config_snapshot$experiment$sim_id <- sim_id
config_snapshot$experiment$run_id <- run_id
write_strict_yaml(config_snapshot, config_snapshot_path)

message("[✓] Saved config snapshot to: ", config_snapshot_path)
