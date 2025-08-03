# applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/helpers/data_utils.R

generate_data <- function(config, lambda_0) {
  iter_seed <- get_seed_for_iter(config$seed, config$experiment$iter_id)
  set.seed(iter_seed)
  
  G <- config$groups$n_groups
  n_per_group <- config$groups$n_per_group
  total_n <- sum(n_per_group)
  
  # Group ID for all observations
  group_id <- rep(LETTERS[seq_len(G)], times = n_per_group)
  
  # Exposure times
  exposure_dist <- match.fun(config$exposure$distribution)
  exposure_args <- config$exposure$args
  t <- do.call(exposure_dist, c(list(n = total_n), exposure_args))
  
  # True group-specific rates
  lambda <- lambda_0[group_id]
  mu <- lambda * t
  Y <- rpois(total_n, mu)
  
  data <- tibble(group = factor(group_id),
                 t = t,
                 Y = Y) |> 
    dplyr::group_by(group) |> 
    dplyr::summarise(t = sum(t),
                     Y = sum(Y))
  
  return(data)
}
