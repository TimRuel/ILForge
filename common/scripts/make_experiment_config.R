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
suppressMessages(i_am("common/scripts/make_experiment_config.R"))

# -------------------------------
# ✅ Parse arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript make_experiment_config.R <app_name> <estimand> <model> <exp_id>")
}
app_name <- args[1]
estimand <- args[2]
model    <- args[3]
exp_id   <- args[4]

# -------------------------------
# ✅ Define paths  
# -------------------------------
model_config_dir   <- here("applications", app_name, estimand, model, "config")
model_config_path  <- here(model_config_dir, "model_config.yml")
exp_config_dir     <- here("config", "exps")
exp_config_path    <- here(exp_config_dir, paste0(exp_id, ".yml"))
exp_dir            <- here("experiments", exp_id)
true_params_dir    <- here(exp_dir, "true_params")
simulations_dir    <- here(exp_dir, "simulations")
common_helpers_dir <- here("common", "scripts", "helpers")

# -------------------------------
# ✅ Load common helpers
# -------------------------------
miceadds::source.all(common_helpers_dir, print.source = FALSE)

# -------------------------------
# ✅ Load and validate template config
# -------------------------------
if (!file_exists(model_config_path)) {
  stop(glue("❌ Model template config not found: {model_config_path}"))
}

exp_config <- read_yaml(model_config_path)

if (!("experiment" %in% names(exp_config))) {
  stop("❌ 'experiment' field is missing from the template YAML.")
}

exp_config$experiment$id <- exp_id

# -------------------------------
# ✅ Write resolved config
# -------------------------------
dir_create(exp_config_dir)
write_strict_yaml(exp_config, exp_config_path)
message(glue("[✓] Created experiment config at: {exp_config_path}"))

# -------------------------------
# ✅ Scaffold experiment folder
# -------------------------------
dir_create(true_params_dir)
dir_create(simulations_dir)
message(glue("[✓] Created experiment scaffold under: {exp_dir}"))

# -------------------------------
# ✅ Summary
# -------------------------------
message("\n📦 Experiment Initialization Complete")
message(glue("  • Application:  {app_name}"))
message(glue("  • Estimand:     {estimand}"))
message(glue("  • Model:        {model}"))
message(glue("  • Experiment:   {exp_id}"))
