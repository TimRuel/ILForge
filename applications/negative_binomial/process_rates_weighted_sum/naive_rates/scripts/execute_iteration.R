#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Execute a single ILForge iteration (Negative Binomial / Process Rates)
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(yaml)
  library(fs)
  library(future)
  library(doFuture)
  library(future.callr)
  library(foreach)
  library(nloptr)
  library(miceadds)
})

# -------------------------------
# ✅ Anchor project root
# -------------------------------
suppressMessages(
  i_am("applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/execute_iteration.R")
)

# -------------------------------
# ✅ Load helpers
# -------------------------------
source_all_helpers <- function(...) {
  dirs <- list(...)
  for (d in dirs) if (dir_exists(d)) miceadds::source.all(d, print.source = FALSE)
}

model_helpers_dir  <- here("applications", "negative_binomial",
                           "process_rates_weighted_sum", "naive_rates", "scripts", "helpers")
common_helpers_dir <- here("common", "scripts", "helpers")

source_all_helpers(model_helpers_dir, common_helpers_dir)

# -------------------------------
# ✅ Parse arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)

iter_dir        <- get_arg_or_stop(1, "iteration directory")
true_params_dir <- get_arg_or_stop(2, "true parameters directory")

check_path_exists(iter_dir, "Iteration directory")
check_path_exists(true_params_dir, "True parameters directory")

skip_integrated <- "--skip-integrated" %in% args
skip_profile    <- "--skip-profile" %in% args

# -------------------------------
# ✅ Start timer
# -------------------------------
overall_start <- Sys.time()

# -------------------------------
# ✅ Load and sanitize config snapshot
# -------------------------------
config_snapshot_path <- here(iter_dir, "config_snapshot.yml")
config <- sanitize_config(read_yaml(config_snapshot_path))

app_name <- config$experiment$app
estimand <- config$experiment$estimand
model    <- config$experiment$model

# -------------------------------
# ✅ Load input data
# -------------------------------
data_dir <- here(iter_dir, "data")
data <- readRDS(here(data_dir, "data.rds"))

# Required true parameter files
required_true <- c("weights.rds")
missing_true <- required_true[!file_exists(here(true_params_dir, required_true))]
if (length(missing_true) > 0)
  stop("Missing true parameter files: ", paste(missing_true, collapse = ", "))

weights <- readRDS(here(true_params_dir, "weights.rds"))

# -------------------------------
# ✅ Setup results directory
# -------------------------------
results_dir <- here(iter_dir, "results")
dir_create(results_dir)

# -------------------------------
# ✅ Register DoFuture backend
# -------------------------------
registerDoFuture()
options(future.globals.maxSize = 1e9)

# -------------------------------
# ✅ Integrated Likelihood
# -------------------------------
il_cfg <- config$optimization_specs$IL %||% list()
num_workers <- il_cfg$num_workers %||% 1L

integrated_path <- here(results_dir, "integrated_LL.rds")

if (!skip_integrated && !file_exists(integrated_path)) {
  message("🔍 Running integrated likelihood...")
  il_start <- Sys.time()
  
  plan(callr, workers = num_workers, gc = TRUE)
  integrated_LL <- get_integrated_LL(config, data, weights)
  plan(sequential)
  
  saveRDS(integrated_LL, file = integrated_path)
  
  il_end <- Sys.time()
  message(sprintf("✅ Integrated likelihood complete (%.2f min)",
                  as.numeric(difftime(il_end, il_start, units = "mins"))))
} else if (file_exists(integrated_path)) {
  message("⏭️  Skipping integrated likelihood (cached result found).")
} else {
  message("⏭️  Skipping integrated likelihood (flag set).")
}

# -------------------------------
# ✅ Profile Likelihood
# -------------------------------
profile_path <- here(results_dir, "profile_LL.rds")

if (!skip_profile && !file_exists(profile_path)) {
  message("📈 Running profile likelihood...")
  pl_start <- Sys.time()
  
  plan(callr, workers = 2, gc = TRUE)
  profile_LL <- get_profile_LL(config, data, weights)
  plan(sequential)
  
  saveRDS(profile_LL, file = profile_path)
  
  pl_end <- Sys.time()
  message(sprintf("✅ Profile likelihood complete (%.2f min)",
                  as.numeric(difftime(pl_end, pl_start, units = "mins"))))
} else if (file_exists(profile_path)) {
  message("⏭️  Skipping profile likelihood (cached result found).")
} else {
  message("⏭️  Skipping profile likelihood (flag set).")
}

# -------------------------------
# ✅ Save resolved config snapshot
# -------------------------------
write_yaml(config, here(results_dir, "resolved_config.yml"))

# -------------------------------
# ✅ Save report objects
# -------------------------------
report_objects <- get_report_objects(iter_dir)
save_list_objects(report_objects, results_dir)

# -------------------------------
# ✅ End timer and log
# -------------------------------
overall_end <- Sys.time()
elapsed_total <- round(as.numeric(difftime(overall_end, overall_start, units = "mins")), 2)

message("⏱️  Total experiment time: ", elapsed_total, " minutes")
message("✓ Results saved to /", sub(".*(/?experiments/.*)", "\\1", results_dir))

# -------------------------------
# ✅ Save timing + metadata
# -------------------------------
git_hash <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

metadata <- list(
  app_name         = app_name,
  estimand         = estimand,
  model            = model,
  exp_id           = config$experiment$id,
  sim_id           = config$experiment$sim_id,
  iter_id          = config$experiment$iter_id,
  slurm_array_id   = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA),
  timestamp_start  = overall_start,
  timestamp_end    = overall_end,
  elapsed_minutes  = elapsed_total,
  git_commit       = git_hash
)

if (exists("il_start") && exists("il_end"))
  metadata$integrated_likelihood_minutes <- round(as.numeric(difftime(il_end, il_start, units = "mins")), 2)
if (exists("pl_start") && exists("pl_end"))
  metadata$profile_likelihood_minutes <- round(as.numeric(difftime(pl_end, pl_start, units = "mins")), 2)

# -------------------------------
# ✅ Save metadata logs
# -------------------------------
log_dir <- here(iter_dir, "logs")
dir_create(log_dir)
save_iter_metadata(metadata, log_dir)

# -------------------------------
# ✅ Graceful teardown
# -------------------------------
plan("sequential")  # ensures all callr workers terminated
invisible(gc())
