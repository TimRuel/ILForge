# applications/poisson/group_rates_weighted_sum/fixed_effects_varying_slope/scripts/helpers/data_utils.R

# Generate covariates for all observations
generate_covariates <- function(config, n) {
  cols <- list()
  
  for (cov in config$model$covariates) {
    dist_fn <- match.fun(cov$distribution$name)
    col <- do.call(dist_fn, c(list(n = n), cov$distribution$args))
    cols[[cov$name]] <- col
  }
  
  covariates <- as_tibble(cols)
  return(covariates)
}

# Build design matrix compatible with stacked Beta (kappa, rho, delta_g, gamma, zeta_g)
get_X <- function(covariates, n_per_group, config) {
  group_labels <- names(n_per_group)
  total_n <- sum(n_per_group)
  G <- length(n_per_group)
  covs <- config$model$covariates
  
  shared_covs  <- Filter(function(c) is.null(c$slope$group_deviation), covs)
  varying_covs <- Filter(function(c) !is.null(c$slope$group_deviation), covs)
  
  # --- Compute column count ---
  n_cols <- 1 + length(shared_covs) + G + length(varying_covs) + G * length(varying_covs)
  X <- matrix(0, nrow = total_n, ncol = n_cols)
  
  # --- Row indices by group ---
  row_idx <- rep(seq_len(G), times = n_per_group)
  
  col_idx <- 1
  
  # 1. Common intercept (kappa)
  X[, col_idx] <- 1
  col_idx <- col_idx + 1
  
  # 2. Shared slopes (rho)
  for (cov in shared_covs) {
    X[, col_idx] <- covariates[[cov$name]]
    col_idx <- col_idx + 1
  }
  
  # 3. Group intercept deviations (delta_g)
  for (g in seq_len(G)) {
    obs_idx <- which(row_idx == g)
    X[obs_idx, col_idx] <- 1
    col_idx <- col_idx + 1
  }
  
  # 4. Baseline slopes for varying covariates (gamma)
  for (cov in varying_covs) {
    X[, col_idx] <- covariates[[cov$name]]
    col_idx <- col_idx + 1
  }
  
  # 5. Group-specific slope deviations (zeta_g)
  for (cov in varying_covs) {
    for (g in seq_len(G)) {
      obs_idx <- which(row_idx == g)
      X[obs_idx, col_idx] <- covariates[[cov$name]][obs_idx]
      col_idx <- col_idx + 1
    }
  }
  
  # --- Column names ---
  colnames_vec <- c(
    config$model$intercepts$common$symbol,
    sapply(shared_covs, function(c) c$slope$baseline$symbol),
    paste0(config$model$intercepts$group_deviation$symbol, "_", group_labels),
    sapply(varying_covs, function(c) c$slope$baseline$symbol),
    unlist(lapply(varying_covs, function(c) paste0(c$slope$group_deviation$symbol, "_", group_labels)))
  )
  
  colnames(X) <- colnames_vec
  return(X)
}

# Generate Poisson outcomes
generate_data <- function(config, true_params_dir) {
  iter_seed <- get_seed_for_iter(config$model$seed, config$experiment$iter_id)
  set.seed(iter_seed)
  
  Beta_0      <- readRDS(file.path(true_params_dir, "Beta_0.rds"))
  n_per_group <- readRDS(file.path(true_params_dir, "n_per_group.rds"))
  
  G <- length(n_per_group)
  total_n <- sum(n_per_group)
  group_labels <- names(n_per_group)
  group_id <- rep(group_labels, times = n_per_group)
  
  # Generate exposures
  exposure_dist <- match.fun(config$model$exposure$distribution$name)
  exposure_args <- config$model$exposure$distribution$args
  t <- do.call(exposure_dist, c(list(n = total_n), exposure_args))
  
  # Generate covariates
  covariates <- generate_covariates(config, total_n)
  
  # Construct design matrix compatible with stacked Beta_0
  X <- get_X(covariates, n_per_group, config)
  
  # Linear predictor and Poisson draws
  eta <- X %*% Beta_0
  mu <- exp(eta)
  Y <- rpois(total_n, t * mu)
  
  # Combine into final data frame
  data <- tibble(
    group = factor(group_id, levels = group_labels),
    t = t
  ) |> 
    bind_cols(covariates) |> 
    add_column(Y = Y)
  
  return(data)
}
