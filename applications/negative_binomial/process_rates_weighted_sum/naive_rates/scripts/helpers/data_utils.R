# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/data_utils.R

#' Generate simulated data for the negative binomial process model
#'
#' @param config A list-like experiment configuration (resolved YAML).
#' @param true_params_dir Directory path containing saved true parameters.
#'
#' @return A tibble with simulated data (columns: process, t, Y).
#' @details
#' This function:
#' 1. Loads the true parameters (theta_0, phi_0, n_per_process).
#' 2. Generates exposure times `t` from the specified distribution.
#' 3. Simulates counts `Y` from a Negative Binomial model with
#'    mean `mu = theta * t` and dispersion `phi`.
#' 4. Returns a tidy tibble of simulated observations.
#'
#' Only randomness depending on the current iteration (`config$experiment$iter_id`)
#' is introduced.
#'
generate_data <- function(config, true_params_dir) {
  # --- Set deterministic seed for reproducibility
  iter_seed <- get_seed_for_iter(config$seed, config$experiment$iter_id)
  set.seed(iter_seed)
  
  # --- Load precomputed true parameters
  theta_0       <- readRDS(file.path(true_params_dir, "theta_0.rds"))
  phi_0         <- readRDS(file.path(true_params_dir, "phi_0.rds"))
  n_per_process <- readRDS(file.path(true_params_dir, "n_per_process.rds"))
  
  # --- Basic quantities
  J <- length(n_per_process)
  total_n <- sum(n_per_process)
  
  # --- Process labels and IDs
  process_labels <- names(n_per_process)
  process_id <- rep(process_labels, times = n_per_process)
  
  # --- Generate exposure times
  exposure_dist <- match.fun(config$exposure$distribution)
  exposure_args <- config$exposure$args
  t <- do.call(exposure_dist, c(list(n = total_n), exposure_args))
  
  # --- Process-specific rate and dispersion
  theta <- theta_0[process_id]
  phi <- phi_0[process_id]
  mu <- theta * t
  
  # --- Simulate counts from Negative Binomial
  Y <- rnbinom(n = total_n, size = phi, mu = mu)
  
  # --- Return tidy data frame
  tibble::tibble(
    process = factor(process_id, levels = process_labels),
    t = t,
    Y = Y
  )
}
