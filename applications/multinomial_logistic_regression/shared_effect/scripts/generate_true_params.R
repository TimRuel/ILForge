# applications/multinomial_logistic_regression/shared_effect/scripts/generate_true_params.R

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
suppressMessages(i_am("applications/multinomial_logistic_regression/shared_effect/scripts/generate_true_params.R"))

# -------------------------------
# ✅ Parse args
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript generate_true_params.R <exp_id>")
exp_id <- args[1]

# -------------------------------
# ✅ Paths
# -------------------------------
estimand_helpers_dir <- here("applications", "multinomial_logistic_regression", "shared_effect", "scripts", "helpers")
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
exp_config <- read_yaml(config_path)

# -------------------------------
# ✅ Generate parameters
# -------------------------------
Beta_0 <- get_Beta_0_from_config(exp_config)
saveRDS(Beta_0, here(true_params_dir, "Beta_0.rds"))
message("[✓] Saved Beta_0 to: ", true_params_dir)
