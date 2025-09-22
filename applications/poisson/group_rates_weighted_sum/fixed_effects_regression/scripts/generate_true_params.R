# applications/poisson/group_rates_weighted_sum/fixed_effects_regression/scripts/generate_true_params.R

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
suppressMessages(i_am("applications/poisson/group_rates_weighted_sum/fixed_effects_regression/scripts/generate_true_params.R"))

# -------------------------------
# ✅ Parse args
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript generate_true_params.R <exp_id>")
exp_id <- args[1]

# -------------------------------
# ✅ Paths
# -------------------------------
model_helpers_dir <- here("applications", "poisson", "group_rates_weighted_sum", "fixed_effects_regression", "scripts", "helpers")
common_helpers_dir  <- here("common", "scripts", "helpers")
exp_config_path <- here("config", "exps", paste0(exp_id, ".yml"))
true_params_dir <- here("experiments", exp_id, "true_params")
dir_create(true_params_dir)

if (!file.exists(exp_config_path)) {
  stop("❌ Experiment config file not found at: ", exp_config_path)
}

# -------------------------------
# ✅ Load helpers
# -------------------------------
miceadds::source.all(common_helpers_dir, print.source = FALSE)
if (dir_exists(model_helpers_dir)) {
  miceadds::source.all(model_helpers_dir, print.source = FALSE)
} else {
  stop("[ERROR] Helper folder not found at ", model_helpers_dir)
}

# -------------------------------
# ✅ Load experiment config
# -------------------------------
exp_config <- read_yaml(exp_config_path)

# -------------------------------
# ✅ Set seed from config for reproducibility
# -------------------------------
if (!is.null(exp_config$model$seed)) {
  message("[INFO] Setting seed from config: ", exp_config$model$seed)
  set.seed(exp_config$model$seed)
} else {
  message("[INFO] No seed found in config; using default RNG state")
}

# -------------------------------
# ✅ Generate parameters
# -------------------------------
true_parameters <- generate_true_parameters(exp_config)

# Validate output elements
expected_names <- c("Beta_0", "theta_0", "weights", "n_per_group")
missing_names <- setdiff(expected_names, names(true_parameters))
if (length(missing_names) > 0) {
  stop("[ERROR] Missing elements in true_parameters: ", paste(missing_names, collapse = ", "))
}

# Save each element as separate .rds file
save_list_objects(true_parameters, true_params_dir)

message("[✓] Saved true parameters to: ", true_params_dir)
message("[INFO] Saved objects: ", paste(names(true_parameters), collapse = ", "))
