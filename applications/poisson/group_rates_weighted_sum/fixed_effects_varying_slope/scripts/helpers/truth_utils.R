# applications/poisson/group_rates_weighted_sum/fixed_effects_varying_slope/scripts/helpers/truth_utils.R

# Generate group labels automatically
generate_group_labels <- function(n_groups, style = "auto_upper") {
  if (is.character(style) && length(style) == 1) {
    if (style == "auto_upper") {
      labels <- character(n_groups)
      for (i in seq_len(n_groups)) {
        num <- i - 1; name <- ""
        repeat {
          name <- paste0(LETTERS[(num %% 26) + 1], name)
          num <- num %/% 26 - 1
          if (num < 0) break
        }
        labels[i] <- name
      }
      return(labels)
    } else if (style == "auto_lower") {
      return(tolower(generate_group_labels(n_groups, "auto_upper")))
    } else if (style == "auto_num") {
      return(as.character(seq_len(n_groups)))
    } else {
      stop(sprintf("Unknown label style: %s", style))
    }
  } else if (length(style) == n_groups) {
    return(as.character(style))
  } else {
    stop("Labels must be a known style keyword or a vector of length n_groups.")
  }
}

# Generic sampler for config blocks like:
#   distribution: { name: rnorm, args: [0,1] }
.sample_from_cfg <- function(dist_cfg, n = 1) {
  if (is.null(dist_cfg)) stop("Distribution config is NULL")
  if (is.null(dist_cfg$name)) stop("Distribution name must be specified in config")
  
  dist_fun <- match.fun(dist_cfg$name)
  args <- dist_cfg$args
  if (is.null(args)) args <- list()
  
  # Combine n with args; if args is a vector, convert to list
  if (is.vector(args) && !is.list(args)) args <- as.list(args)
  out <- do.call(dist_fun, c(list(n = n), args))
  
  return(out)
}

.build_beta_names <- function(cfg, labels) {
  covs <- cfg$model$covariates
  G <- length(labels)
  
  names_vec <- character(0)
  
  # 1. Common intercept
  names_vec <- c(names_vec, cfg$model$intercepts$common$symbol)
  
  # 2. Shared slopes
  shared_covs <- covs[sapply(covs, function(c) is.null(c$slope$group_deviation))]
  if (length(shared_covs) > 0) {
    names_vec <- c(names_vec, sapply(shared_covs, function(c) c$slope$baseline$symbol))
  }
  
  # 3. Group intercept deviations
  names_vec <- c(names_vec, paste0(cfg$model$intercepts$group_deviation$symbol, "_", labels))
  
  # 4. Baseline slopes for varying covariates
  varying_covs <- covs[!sapply(covs, function(c) is.null(c$slope$group_deviation))]
  if (length(varying_covs) > 0) {
    names_vec <- c(names_vec, sapply(varying_covs, function(c) c$slope$baseline$symbol))
  }
  
  # 5. Group-specific slope deviations
  if (length(varying_covs) > 0) {
    for (c in varying_covs) {
      names_vec <- c(names_vec, paste0(c$slope$group_deviation$symbol, "_", labels))
    }
  }
  
  return(names_vec)
}

# Expand n_per_group spec
expand_n_per_group <- function(n_per_group, n_groups) {
  if (length(n_per_group) == 1) return(rep(n_per_group, n_groups))
  if (length(n_per_group) < n_groups) return(rep_len(n_per_group, n_groups))
  if (length(n_per_group) == n_groups) return(n_per_group)
  stop("n_per_group must be length 1, a divisor of n_groups, or exactly n_groups.")
}

# Generate group weights  (reads from cfg$model$weights)
generate_group_weights <- function(cfg) {
  weight_cfg <- cfg$model$weights
  n_groups   <- cfg$model$groups$n_groups
  
  if (is.null(weight_cfg$distribution) || is.null(weight_cfg$distribution$name))
    stop("Weight distribution must be specified as list(distribution=list(name=..., args=[...])).")
  
  raw <- .sample_from_cfg(weight_cfg$distribution, n_groups)
  
  if (!is.null(weight_cfg$normalize_mean_to)) {
    weights <- raw / mean(raw) * weight_cfg$normalize_mean_to
  } else if (!is.null(weight_cfg$normalize_sum_to)) {
    weights <- raw / sum(raw) * weight_cfg$normalize_sum_to
  } else {
    weights <- raw
  }
  
  labels <- generate_group_labels(n_groups, cfg$model$groups$labels)
  names(weights) <- labels
  weights
}

# -------------------------------------------------------------------
# compute_true_marginal_rates(cfg, Beta_0, n_mc = 100000, seed = NULL)
# Monte Carlo approximation of theta_g = E[ exp(eta_g(X)) ]
# Assumes covariates in cfg$model$covariates are independent and each
# cov has a distribution block: distribution: { name: "...", args: [...] }
# -------------------------------------------------------------------
compute_true_marginal_rates <- function(cfg, Beta_0, n_mc = 1e5) {
  G <- cfg$model$groups$n_groups
  labels <- generate_group_labels(G, cfg$model$groups$labels)
  covs <- cfg$model$covariates
  
  # ----- Extract coefficients from Beta_0 -----
  # Order: kappa, rho (shared), delta_g, gamma (baseline), zeta_g
  idx <- 1
  kappa <- Beta_0[idx, 1]; idx <- idx + 1
  
  # Shared slopes
  shared_covs <- covs[sapply(covs, function(c) is.null(c$slope$group_deviation))]
  rho <- if(length(shared_covs) > 0) Beta_0[idx:(idx + length(shared_covs) - 1), 1] else numeric(0)
  idx <- idx + length(shared_covs)
  
  # Group intercept deviations
  delta_g <- Beta_0[idx:(idx + G - 1), 1]; idx <- idx + G
  
  # Baseline slopes for varying covariates
  varying_covs <- covs[!sapply(covs, function(c) is.null(c$slope$group_deviation))]
  gamma <- if(length(varying_covs) > 0) Beta_0[idx:(idx + length(varying_covs) - 1), 1] else numeric(0)
  idx <- idx + length(varying_covs)
  
  # Group-specific slope deviations
  zeta_g <- if(length(varying_covs) > 0) Beta_0[idx:nrow(Beta_0), 1] else numeric(0)
  
  # ----- Generate Monte Carlo samples of covariates -----
  X_mc <- list()
  for (cov in covs) {
    X_mc[[cov$name]] <- .sample_from_cfg(cov$distribution, n_mc)
  }
  
  # ----- Compute marginal rates per group -----
  theta_0 <- numeric(G)
  for (g in seq_len(G)) {
    eta <- rep(kappa + delta_g[g], n_mc)  # common intercept + group dev
    
    # Add shared slopes
    if(length(shared_covs) > 0) {
      for (i in seq_along(shared_covs)) {
        cov_name <- shared_covs[[i]]$name
        eta <- eta + rho[i] * X_mc[[cov_name]]
      }
    }
    
    # Add varying slopes
    if(length(varying_covs) > 0) {
      for (i in seq_along(varying_covs)) {
        cov_name <- varying_covs[[i]]$name
        eta <- eta + (gamma[i] + zeta_g[(i-1)*G + g]) * X_mc[[cov_name]]
      }
    }
    
    # Apply Poisson link
    theta_0[g] <- mean(exp(eta))
  }
  
  names(theta_0) <- labels
  return(theta_0)
}

# -------------------------------------------------------------------
# generate_true_parameters(cfg, n_mc = 100000)
#   - builds Beta_0 stacked vector in order:
#       ("alpha", shared covariate rows..., "delta_<G>", baseline_<varying covs>..., "nu_<cov>_<G>"...)
#   - calls compute_true_marginal_rates() to get theta_0
# -------------------------------------------------------------------
generate_true_parameters <- function(cfg, n_mc = 1e5) {
  G <- cfg$model$groups$n_groups
  labels <- generate_group_labels(G, cfg$model$groups$labels)
  covs <- cfg$model$covariates
  
  # --- 1. Common intercept ---
  kappa <- .sample_from_cfg(cfg$model$intercepts$common$distribution, 1)
  
  # --- 2. Shared slopes ---
  shared_covs <- Filter(function(c) is.null(c$slope$group_deviation), covs)
  rho <- sapply(shared_covs, function(c) .sample_from_cfg(c$slope$baseline$distribution, 1))
  
  # --- 3. Group intercept deviations ---
  delta_g <- .sample_from_cfg(cfg$model$intercepts$group_deviation$distribution, G)
  
  # --- 4. Baseline slopes for varying covariates ---
  varying_covs <- Filter(function(c) !is.null(c$slope$group_deviation), covs)
  gamma <- sapply(varying_covs, function(c) .sample_from_cfg(c$slope$baseline$distribution, 1))
  
  # --- 5. Group-specific slope deviations ---
  zeta_g <- c()
  for (c in varying_covs) {
    zeta_g <- c(zeta_g, .sample_from_cfg(c$slope$group_deviation$distribution, G))
  }
  
  # --- 6. Stack Beta vector ---
  Beta_vals <- c(kappa, rho, delta_g, gamma, zeta_g)
  
  # --- 7. Assign rownames with symbols and "_group" suffixes ---
  Beta_names <- c(
    cfg$model$intercepts$common$symbol,                                    # kappa
    sapply(shared_covs, function(c) c$slope$baseline$symbol),              # rho
    paste0(cfg$model$intercepts$group_deviation$symbol, "_", labels),      # delta_g
    sapply(varying_covs, function(c) c$slope$baseline$symbol),             # gamma
    unlist(lapply(varying_covs, function(c) paste0(c$slope$group_deviation$symbol, "_", labels)))  # zeta_g
  )
  
  Beta_0 <- matrix(Beta_vals, ncol = 1)
  rownames(Beta_0) <- Beta_names
  
  # --- 8. Weights & n_per_group ---
  weights <- generate_group_weights(cfg)
  n_per_group <- expand_n_per_group(cfg$model$groups$n_per_group$args, G)
  n_per_group <- setNames(n_per_group, labels)
  
  # --- 9. True marginal group rates via Monte Carlo ---
  theta_0 <- compute_true_marginal_rates(cfg, Beta_0, n_mc = n_mc)
  
  # --- 10. Return ---
  list(
    Beta_0      = Beta_0,     # stacked coefficient vector
    theta_0     = theta_0,    # marginal rates per group
    weights     = weights,    # named by group
    n_per_group = n_per_group
  )
}
