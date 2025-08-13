# applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/helpers/data_utils.R

generate_data <- function(config, true_params_dir) {
  iter_seed <- get_seed_for_iter(config$seed, config$experiment$iter_id)
  set.seed(iter_seed)
  
  # Load precomputed truth
  theta_0      <- readRDS(file.path(true_params_dir, "theta_0.rds"))
  n_per_group  <- readRDS(file.path(true_params_dir, "n_per_group.rds"))
  
  G <- length(n_per_group)
  total_n <- sum(n_per_group)
  
  # Use saved labels for group IDs
  group_labels <- names(n_per_group)
  group_id <- rep(group_labels, times = n_per_group)
  
  # Exposure times
  exposure_dist <- match.fun(config$exposure$distribution)
  exposure_args <- config$exposure$args
  t <- do.call(exposure_dist, c(list(n = total_n), exposure_args))
  
  # True group-specific rates
  theta <- theta_0[group_id]
  mu <- theta * t
  Y <- rpois(total_n, mu)
  
  # Aggregate to group level
  data <- tibble(group = factor(group_id, levels = group_labels),
                 t = t,
                 Y = Y) |> 
    dplyr::group_by(group) |> 
    dplyr::summarise(t = sum(t),
                     Y = sum(Y),
                     .groups = "drop")
  
  return(data)
}

