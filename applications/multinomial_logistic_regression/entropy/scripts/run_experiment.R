# applications/multinomial_logistic_regression/entropy/scripts/run_experiment.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(yaml)
  library(fs)
  library(doFuture)
})

# -------------------------------
# ✅ Anchor project root
# -------------------------------
i_am("applications/multinomial_logistic_regression/entropy/scripts/run_experiment.R")

# -------------------------------
# ✅ Parse arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args) > 0) args[1] else stop("[ERROR] A run directory was not provided.")
if (!file.exists(run_dir)) {
  stop("[ERROR] Run directory does not exist at /", sub(".*(/?experiments/.*)", "\\1", run_dir))
}
skip_integrated <- "--skip-integrated" %in% args
skip_profile    <- "--skip-profile" %in% args

# -------------------------------
# ✅ Start timer
# -------------------------------
overall_start <- Sys.time()

# -------------------------------
# ✅ Load config snapshot
# -------------------------------
config_snapshot_path <- here("config_snapshot.yml")
config <- read_yaml(config_snapshot_path)
app_name <- config$experiment$app_name
estimand <- config$experiment$estimand

# -------------------------------
# ✅ Load helpers
# -------------------------------
common_helpers_dir <- here("common", "scripts", "helpers")
miceadds::source.all(common_helpers_dir, print.source = FALSE)
estimand_helpers_dir <- here("applications", app_name, estimand, "scripts", "helpers")
if (dir_exists(helper_dir)) {
  miceadds::source.all(estimand_helpers_dir, print.source = FALSE)
} else {
  stop("[ERROR] Helper folder not found at ", estimand_helpers_dir)
}

# -------------------------------
# ✅ Load input data
# -------------------------------
data_dir <- here(run_dir, "data")
required_files <- c("X_design.rds", "model_df.rds")
missing <- required_files[!file_exists(here(data_dir, required_files))]
if (length(missing) > 0) {
  stop("Missing required data files: ", paste(missing, collapse = ", "))
}
X_design <- readRDS(here(data_dir, "X_design.rds"))
model_df <- readRDS(here(data_dir, "model_df.rds"))

# -------------------------------
# ✅ Setup results directory
# -------------------------------
results_dir <- here(run_dir, "results")
dir_create(results_dir)

# -------------------------------
# ✅ Run integrated likelihood
# -------------------------------
if (!skip_integrated) {
  message("🔍 Running integrated likelihood...")
  il_start <- Sys.time()
  
  num_workers <- config$optimization_specs$IL$num_workers
  plan_strategy <- if (.Platform$OS.type == "unix") multicore else multisession
  plan(plan_strategy, workers = I(num_workers))
  
  integrated_LL <- get_integrated_LL(config, X_design, model_df)
  plan(sequential)
  
  saveRDS(integrated_LL, file = here(results_dir, "integrated_LL.rds"))
  
  il_end <- Sys.time()
  message(sprintf("✅ Integrated likelihood complete (%.2f min)", as.numeric(difftime(il_end, il_start, units = "mins"))))
} else {
  message("⏭️  Skipping integrated likelihood.")
}

# -------------------------------
# ✅ Run profile likelihood
# -------------------------------
if (!skip_profile) {
  message("📈 Running profile likelihood...")
  pl_start <- Sys.time()
  
  plan_strategy <- if (.Platform$OS.type == "unix") multicore else multisession
  plan(plan_strategy, workers = I(2))
  
  profile_LL <- get_profile_LL(config, X_design, model_df)
  plan(sequential)
  
  saveRDS(profile_LL, file = here(results_dir, "profile_LL.rds"))
  
  pl_end <- Sys.time()
  message(sprintf("✅ Profile likelihood complete (%.2f min)", as.numeric(difftime(pl_end, pl_start, units = "mins"))))
} else {
  message("⏭️  Skipping profile likelihood.")
}

# -------------------------------
# ✅ Save report objects
# -------------------------------
report_objects <- get_report_objects(run_dir)
save_list_objects(report_objects, results_dir)

# -------------------------------
# ✅ End timer and log
# -------------------------------
overall_end <- Sys.time()
elapsed_total <- round(as.numeric(difftime(overall_end, overall_start, units = "mins")), 2)
message("⏱️  Total experiment time: ", elapsed_total, " minutes")
message("✓ Experiment results saved to /", sub(".*(/?experiments/.*)", "\\1", results_dir))

# -------------------------------
# ✅ Save timing metadata
# -------------------------------
git_hash <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  error = function(e) NA_character_
)

metadata <- list(
  app_name         = app_name,
  estimand         = estimand,
  experiment_id    = config$experiment$id,
  sim_id           = config$experiment$sim_id,
  run_id           = config$experiment$run_id,
  slurm_array_id   = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA),
  timestamp_start  = overall_start,
  timestamp_end    = overall_end,
  elapsed_minutes  = round(as.numeric(difftime(overall_end, overall_start, units = "mins")), 2),
  git_commit       = git_hash
)

# Optional durations (if not skipped)
if (!skip_integrated) {
  metadata$integrated_likelihood_minutes <- round(as.numeric(difftime(il_end, il_start, units = "mins")), 2)
}
if (!skip_profile) {
  metadata$profile_likelihood_minutes <- round(as.numeric(difftime(pl_end, pl_start, units = "mins")), 2)
}

# -------------------------------
# ✅ Save metadata
# -------------------------------
log_dir <- here(run_dir, "logs")
dir_create(log_dir)
save_run_metadata(metadata, log_dir)

