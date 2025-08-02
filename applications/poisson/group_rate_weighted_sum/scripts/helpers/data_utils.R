# applications/poisson_regression/weighted_sum/scripts/helpers/data_utils.R

generate_data <- function(config, Beta_0, covariates) {
  
  iter_seed <- get_seed_for_iter(config$seed, config$experiment$iter_id)
  set.seed(iter_seed)
  
  G <- config$group_structure$n_groups
  n_per_group <- config$group_structure$n_per_group
  total_n <- G * n_per_group
  
  # Group ID for all observations
  group_id <- rep(1:G, each = n_per_group)
  
  # Expand group-level covariates to observation-level
  X <- covariates[group_id, , drop = FALSE]
  
  # Exposure time (observation-level)
  exposure_dist <- match.fun(config$exposure$distribution)
  exposure_args <- config$exposure$args
  t <- do.call(exposure_dist, c(list(n = total_n), exposure_args))
  
  # Linear predictor and Poisson response
  lambda <- numeric(total_n)
  Y <- numeric(total_n)
  
  for (g in 1:G) {
    idx <- which(group_id == g)
    Xg <- as.numeric(covariates[g, , drop = FALSE])  # group-level covariates
    bg <- if (is.matrix(Beta_0)) Beta_0[, g] else Beta_0
    eta <- as.numeric(Xg %*% bg)
    lambda[idx] <- exp(eta)  # shared for group
    mu <- lambda[idx] * t[idx]
    Y[idx] <- rpois(length(mu), mu)
  }
  
  data <- tibble(group = factor(group_id),
                 t = t,
                 Y = Y) |> 
    bind_cols(as_tibble(X)) |>
    # group_by(group) |>
    # summarise(
    #   Y = sum(Y),
    #   t = sum(t),
    #   across(starts_with("X"), first)
    # ) |> 
    select(group, starts_with("X"), t, Y)
  
  return(data)
}
