# common/scripts/make_experiment_config.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(fs)
  library(glue)
})

# -------------------------------
# ✅ Anchor project root
# -------------------------------
i_am("common/scripts/make_experiment_config.R")

# -------------------------------
# ✅ Parse arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript make_experiment_config.R <app_name> <estimand> <experiment_id>")
}

app_name      <- args[1]
estimand      <- args[2]
experiment_id <- args[3]

# -------------------------------
# ✅ Define paths  
# -------------------------------
estimand_config_dir <- here("applications", app_name, estimand, "config")
template_path       <- path(estimand_config_dir, "exp_config.yml")

exp_config_dir      <- here("config", "exps")
resolved_config     <- path(exp_config_dir, paste0(experiment_id, ".yml"))

experiment_dir      <- here("experiments", experiment_id)
true_params_dir     <- path(experiment_dir, "true_params")
individual_runs_dir <- path(experiment_dir, "individual_runs")
simulations_dir     <- path(experiment_dir, "simulations")

# -------------------------------
# ✅ Load and validate template config
# -------------------------------
if (!file_exists(template_path)) {
  stop(glue("❌ Template config not found: {template_path}"))
}

experiment_config <- read_yaml(template_path)

if (!("experiment" %in% names(experiment_config))) {
  stop("❌ 'experiment' field is missing from the template YAML.")
}

experiment_config$experiment$id <- experiment_id

# -------------------------------
# ✅ Write resolved config
# -------------------------------
dir_create(exp_config_dir)
write_strict_yaml(experiment_config, resolved_config)
message(glue("[✓] Created resolved config at: {resolved_config}"))

# -------------------------------
# ✅ Scaffold experiment folder
# -------------------------------
dir_create(true_params_dir)
dir_create(individual_runs_dir)
dir_create(simulations_dir)

message(glue("[✓] Created experiment scaffold under: {experiment_dir}"))

# -------------------------------
# ✅ Summary
# -------------------------------
message("\n📦 Experiment Initialization Complete")
message(glue("  • Application:  {app_name}"))
message(glue("  • Estimand:     {estimand}"))
message(glue("  • Experiment:   {experiment_id}"))
