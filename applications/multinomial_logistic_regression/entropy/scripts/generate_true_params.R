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
i_am("applications/multinomial_logistic_regression/entropy/scripts/generate_true_params.R")

# -------------------------------
# ✅ Load helpers
# -------------------------------
helper_dir <- here("applications", "multinomial_logistic_regression", "entropy", "scripts", "helpers")
miceadds::source.all(helper_dir, print.source = FALSE)

# -------------------------------
# ✅ Parse args
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript generate_true_params.R <experiment_id>")
experiment_id <- args[1]

# -------------------------------
# ✅ Paths
# -------------------------------
config_path <- here("config", "exps", paste0(experiment_id, ".yml"))
true_params_dir <- here("experiments", experiment_id, "true_params")
dir_create(true_params_dir)

if (!file.exists(config_path)) {
  stop("❌ Config file not found at: ", config_path)
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
