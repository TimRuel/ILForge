# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/data_utils.R

generate_data <- function(config, true_params_dir) {
  iter_seed <- get_seed_for_iter(config$seed, config$experiment$iter_id)
  set.seed(iter_seed)
  
  # Load precomputed truth
  theta_0       <- readRDS(file.path(true_params_dir, "theta_0.rds"))
  phi_0         <- readRDS(file.path(true_params_dir, "phi_0.rds"))
  n_per_process <- readRDS(file.path(true_params_dir, "n_per_process.rds"))
  
  J <- length(n_per_process)
  total_n <- sum(n_per_process)
  
  # Use saved labels for process IDs
  process_labels <- names(n_per_process)
  process_id <- rep(process_labels, times = n_per_process)
  
  # Exposure times
  exposure_dist <- match.fun(config$exposure$distribution)
  exposure_args <- config$exposure$args
  t <- do.call(exposure_dist, c(list(n = total_n), exposure_args))
  
  # True process-specific rates
  theta <- theta_0[process_id]
  phi <- phi_0[process_id]
  mu <- theta * t
  Y <- rnbinom(n = total_n, size = phi, mu = mu)
  
  # Aggregate to process level
  data <- tibble(process = factor(process_id, levels = process_labels),
                 t = t,
                 Y = Y)
  return(data)
}

