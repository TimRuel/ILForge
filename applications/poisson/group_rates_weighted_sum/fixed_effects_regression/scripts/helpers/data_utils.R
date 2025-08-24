# applications/poisson/group_rates_weighted_sum/fixed_effects_regression/scripts/helpers/data_utils.R

generate_covariates <- function(config, n) {
  
  cols <- list()
  
  for (cov in config$model$covariates) {
    if (cov$type == "categorical") {
      col <- sample(cov$levels, size = n, replace = TRUE)
    } else {
      dist_fn <- match.fun(cov$distribution)
      col <- do.call(dist_fn, c(list(n = n), cov$args))
    }
    cols[[cov$name]] <- col
  }
  
  covariates <- as_tibble(cols)
  return(covariates)
}

get_X <- function(covariates, n_per_group) {
  
  group_labels <- names(n_per_group)
  
  df <- tibble(group = factor(rep(group_labels, times = n_per_group))) |> 
    add_column(covariates)
  
  X <- model.matrix(~ . - 1 , data = df)
  colnames(X) <- c(group_labels, colnames(covariates))
  
  return(X)
}

generate_data <- function(config, true_params_dir) {
  
  iter_seed <- get_seed_for_iter(config$seed, config$experiment$iter_id)
  set.seed(iter_seed)
  
  Beta_0      <- readRDS(file.path(true_params_dir, "Beta_0.rds"))
  n_per_group <- readRDS(file.path(true_params_dir, "n_per_group.rds"))
  
  G <- length(n_per_group)
  total_n <- sum(n_per_group)
  
  group_labels <- names(n_per_group)
  group_id <- rep(group_labels, times = n_per_group)
  
  exposure_dist <- match.fun(config$exposure$distribution)
  exposure_args <- config$exposure$args
  t <- do.call(exposure_dist, c(list(n = total_n), exposure_args))
  
  covariates <- generate_covariates(config, total_n)
  X <- get_X(covariates, n_per_group)
  
  eta <- X %*% Beta_0
  mu <- exp(eta)
  Y <- rpois(total_n, t * mu)
  
  data <- tibble(group = factor(group_id, levels = group_labels),
                 t = t) |> 
    bind_cols(covariates) |> 
    add_column(Y = Y)
  
  return(data)
}

