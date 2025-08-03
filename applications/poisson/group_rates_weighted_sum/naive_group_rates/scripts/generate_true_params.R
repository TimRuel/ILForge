# applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/generate_true_params.R

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
suppressMessages(i_am("applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/generate_true_params.R"))

# -------------------------------
# ✅ Parse args
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript generate_true_params.R <exp_id>")
exp_id <- args[1]

# -------------------------------
# ✅ Paths
# -------------------------------
model_helpers_dir <- here("applications", "poisson", "group_rates_weighted_sum", "naive_group_rates", "scripts", "helpers")
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
# ✅ Generate parameters
# -------------------------------
true_parameters <- generate_true_parameters(exp_config)

# Should return a named list like: list(lambda = ..., weights = ...)
save_list_objects(true_parameters, true_params_dir)

message("[✓] Saved true parameters to: ", true_params_dir)
