# common/scripts/main.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(doFuture)
  library(here)
  library(quarto)
  library(yaml)
  library(fs)
  library(nloptr)
})

# -------------------------------
# ✅ Anchor project root
# -------------------------------
suppressMessages(i_am("common/scripts/main.R"))

# -------------------------------
# ✅ Parse arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript main.R <app_name> <estimand> <model> <exp_id> <sim_id> <iter_id> [cores] [--skip-integrated] [--skip-profile]")
}

app_name <- args[1]
estimand <- args[2]
model    <- args[3]
exp_id   <- args[4]
sim_id   <- args[5]
iter_id  <- args[6]
requested_cores <- if (length(args) >= 7 && !grepl("^--", args[7])) as.integer(args[7]) else NULL

skip_integrated <- "--skip-integrated" %in% args
skip_profile    <- "--skip-profile" %in% args
flag_args <- c(
  if (skip_integrated) "--skip-integrated",
  if (skip_profile) "--skip-profile"
)

# -------------------------------
# ✅ Path helpers
# -------------------------------
common_path <- function(...) here("common", ...)
common_script <- function(...) common_path("scripts", ...)
model_path <- function(...) here("applications", app_name, estimand, model, ...)
model_script <- function(...) model_path("scripts", ...)

# -------------------------------
# ✅ Load helpers (if present)
# -------------------------------
common_helpers_dir <- common_script("helpers")
miceadds::source.all(common_helpers_dir, print.source = FALSE)

model_helpers_dir <- model_script("helpers")
if (dir_exists(model_helpers_dir)) {
  miceadds::source.all(model_helpers_dir, print.source = FALSE)
} else {
  stop("Aborting due to missing model helpers directory.")
}

# -------------------------------
# ✅ Utility to call R scripts
# -------------------------------
run_common_script <- function(script_name, args = character()) {
  script_path <- common_script(script_name)
  system2("Rscript", c(script_path, args))
}

run_model_script <- function(script_name, args = character()) {
  script_path <- model_script(script_name)
  system2("Rscript", c(script_path, args))
}

# --------------------------------------
# ✅ Ensure experiment config exists
# --------------------------------------
exp_config_path <- here("config", "exps", paste0(exp_id, ".yml"))
if (!file.exists(exp_config_path)) {
  message("Creating experiment config...")
  run_common_script("make_experiment_config.R", c(app_name, estimand, model, exp_id))
}

# -------------------------------
# ✅ Generate true parameters
# -------------------------------
true_params_dir <- here("experiments", exp_id, "true_params")
if (length(dir_ls(true_params_dir, fail = FALSE)) == 0) {
  message("[INFO] No true parameters found — generating...")
  run_model_script("generate_true_params.R", exp_id)
} else {
  message("[INFO] True parameters already exist — skipping generation.")
}

# -------------------------------
# ✅ Set iteration directory
# -------------------------------
iter_dir <- here("experiments", exp_id, "simulations", sim_id, iter_id)
dir_create(iter_dir)

# -------------------------------
# ✅ Generate data
# -------------------------------
message("Generating data...")
run_model_script("generate_data.R", c(exp_id, sim_id, iter_id))

# -------------------------------
# ✅ Add optimization config
# -------------------------------
config_snapshot_path <- here(iter_dir, "config_snapshot.yml")
config_snapshot <- read_yaml(config_snapshot_path)

opt_config_path <- model_path("config", "opt_config.yml")
opt_config <- read_yaml(opt_config_path)
config_snapshot$optimization_specs <- c(config_snapshot$optimization_specs, opt_config)

# -------------------------------
# ✅ Set core usage
# -------------------------------
core_info <- get_core_config(requested_cores)
config_snapshot$optimization_specs$IL$max_cores <- core_info$max_cores
config_snapshot$optimization_specs$IL$num_workers <- core_info$num_workers

# -------------------------------
# ✅ Save updated config
# -------------------------------
write_strict_yaml(config_snapshot, config_snapshot_path)

# -------------------------------
# ✅ Run experiment
# -------------------------------
message("Executing iteration...")
run_model_script("execute_iteration.R", c(iter_dir, true_params_dir, flag_args))
message("✓ Iteration completed")
