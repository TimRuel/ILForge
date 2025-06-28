# applications/multinomial_logistic_regression/entropy/scripts/generate_true_params.R

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
suppressMessages(i_am("applications/multinomial_logistic_regression/entropy/scripts/generate_true_params.R"))

# -------------------------------
# ✅ Parse args
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript generate_true_params.R <exp_id>")
exp_id <- args[1]

# -------------------------------
# ✅ Paths
# -------------------------------
estimand_helpers_dir <- here("applications", "multinomial_logistic_regression", "entropy", "scripts", "helpers")
common_helpers_dir  <- here("common", "scripts", "helpers")
config_path <- here("config", "exps", paste0(exp_id, ".yml"))
true_params_dir <- here("experiments", exp_id, "true_params")
dir_create(true_params_dir)
if (!file.exists(config_path)) {
  stop("❌ Config file not found at: ", config_path)
}

# -------------------------------
# ✅ Load helpers
# -------------------------------
miceadds::source.all(common_helpers_dir, print.source = FALSE)
if (dir_exists(estimand_helpers_dir)) {
  miceadds::source.all(estimand_helpers_dir, print.source = FALSE)
} else {
  stop("[ERROR] Helper folder not found at ", estimand_helpers_dir)
}

# -------------------------------
# ✅ Load config
# -------------------------------
experiment_config <- read_yaml(config_path)
X1_levels         <- experiment_config$X1_levels
model_specs       <- experiment_config$model_specs
opt_specs         <- experiment_config$optimization_specs
seed              <- opt_specs$seed %||% 1

# -------------------------------
# ✅ Generate parameters
# -------------------------------
set.seed(seed)
experiment_parameters <- get_experiment_parameters(X1_levels, model_specs)

save_list_objects(experiment_parameters$true_params, true_params_dir)
message("[✓] Saved true parameters to: ", true_params_dir)

# -------------------------------
# ✅ Update config with derived model specs
# -------------------------------
experiment_config$model_specs <- experiment_parameters$model_specs
write_strict_yaml(experiment_config, config_path)
message("[✓] Updated experiment config with additional model specs.")
