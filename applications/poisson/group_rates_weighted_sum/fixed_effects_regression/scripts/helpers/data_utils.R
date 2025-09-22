# applications/poisson/group_rates_weighted_sum/fixed_effects_regression/scripts/helpers/data_utils.R

# Generate covariates for all observations
generate_covariate <- function(cov_cfg, n) {
  dist_fn <- match.fun(cov_cfg$distribution$name)
  covariate <- do.call(dist_fn, c(list(n = n), cov_cfg$distribution$args))
  return(covariate)
}

recover_original_covariates <- function(design_matrix, 
                                        intercept_prefix = "α", 
                                        homo_cov_prefix = "γ",
                                        hetero_cov_prefix = "ζ") {
  # design_matrix: data.frame or matrix
  # intercept_prefix: prefix for the one-hot group indicator columns
  # homo_cov_prefix: vector of names for covariates that are common across groups
  # hetero_cov_prefix: prefix for the group-specific covariate columns
  
  coefs <- colnames(design_matrix)
  
  homo_cov_coefs <- coefs[grepl(homo_cov_prefix, coefs)]
  
  num_homo_covs <- length(homo_cov_coefs)
  
  homo_covs <- design_matrix[, homo_cov_coefs, drop = FALSE]
  
  colnames(homo_covs) <- paste0("X", 1:num_homo_covs)
  
  G <- design_matrix |> 
    rownames() |> 
    unique() |> 
    length()
  
  num_hetero_covs <- sum(grepl(hetero_cov_prefix, coefs)) / G
  
  hetero_cov_coefs <- coefs[grepl(hetero_cov_prefix, coefs)]
  
  hetero_cov_dummy_mat <- design_matrix[, hetero_cov_coefs, drop = FALSE]
  
  hetero_covs <- matrix(NA, nrow = nrow(hetero_cov_dummy_mat), ncol = num_hetero_covs)
  
  for (i in 1:num_hetero_covs) {
    
    hetero_covs[, i] <- rowSums(hetero_cov_dummy_mat[, i:(i + G - 1), drop = FALSE]) 
  }
  
  colnames(hetero_covs) <- paste0("X", (num_homo_covs + 1):(num_homo_covs + num_hetero_covs))
  
  # 4. Combine into a single data.frame
  recovered <- cbind(homo_covs, hetero_covs) |> 
    as.data.frame(row.names = 1:nrow(design_matrix))
  
  return(recovered)
}

# Build design matrix compatible with stacked Beta (kappa, gamma, delta_g, zeta_g)
get_X <- function(config, Beta_0, n_per_group) {
  
  group_labels <- names(n_per_group)
  total_n <- sum(n_per_group)
  G <- length(n_per_group)
  cov_cfgs <- config$model$covariates
  row_idx <- rep(seq_len(G), times = n_per_group)
  X <- matrix(0, nrow = total_n, ncol = nrow(Beta_0))
  colnames(X) <- rownames(Beta_0)
  rownames(X) <- rep(names(n_per_group), times = n_per_group)
  
  for (symbol in colnames(X)) {
    
    if (grepl("α", symbol)) {
      
      g <- symbol |> 
        str_sub(-1, -1) |> 
        (\(g) which(LETTERS == g))()
      obs_idx <- which(row_idx == g)
      X[obs_idx, symbol] <- 1
    }
    
    else if (grepl("γ", symbol)) {
      
      cov_cfg <- Filter(function(c) c$coefficient$symbol == symbol, cov_cfgs)[[1]]
      X[, symbol] <- generate_covariate(cov_cfg$variable, total_n)
    }
    
    else if (grepl("ζ", symbol)) {
      
      g <- symbol |> 
        str_sub(-1, -1) |> 
        (\(g) which(LETTERS == g))()
      obs_idx <- which(row_idx == g)
      cov_cfg <- Filter(function(c) c$coefficient$symbol == str_sub(symbol, 1, 2), cov_cfgs)[[1]]
      X[obs_idx, symbol] <- generate_covariate(cov_cfg$variable, n_per_group[g])
    }
  }
  return(X)
}

# Generate Poisson outcomes
generate_data <- function(config, Beta_0) {
  
  n_per_group <- expand_groups(config$model$groups)
  group_labels <- names(n_per_group)
  G <- length(n_per_group)
  total_n <- sum(n_per_group)
  group_id <- rep(group_labels, times = n_per_group)
  
  exposure_dist <- match.fun(config$model$exposure$distribution$name)
  exposure_args <- config$model$exposure$distribution$args
  t <- do.call(exposure_dist, c(list(n = total_n), exposure_args))
  
  X <- get_X(config, Beta_0, n_per_group)
  eta <- X %*% Beta_0
  mu <- exp(eta)
  Y <- rpois(total_n, t * mu)
  
  covariates <- recover_original_covariates(X)
  
  # Combine into final data frame
  data <- tibble(
    group = factor(group_id, levels = group_labels),
    t = t
  ) |> 
    bind_cols(covariates) |> 
    add_column(Y = Y)
  
  return(list(data = data,
              X = X))
}
