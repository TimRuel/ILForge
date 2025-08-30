# Generate group labels automatically
generate_group_labels <- function(n_groups, style = "auto_upper") {
  if (is.character(style) && length(style) == 1) {
    if (style == "auto_upper") {
      labels <- character(n_groups)
      for (i in seq_len(n_groups)) {
        num <- i - 1
        name <- ""
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

# Generate group weights
generate_group_weights <- function(cfg) {
  weight_cfg <- cfg$weights
  n_groups <- cfg$model$groups$n_groups
  
  if (is.null(weight_cfg$distribution) || is.null(weight_cfg$args)) {
    stop("Weight distribution and args must be specified in config.")
  }
  
  dist_fun <- match.fun(weight_cfg$distribution)
  raw_weights <- do.call(dist_fun, c(list(n = n_groups), as.list(weight_cfg$args)))
  
  if (!is.null(weight_cfg$normalize_mean_to)) {
    weights <- raw_weights / mean(raw_weights) * weight_cfg$normalize_mean_to
  } else if (!is.null(weight_cfg$normalize_sum_to)) {
    weights <- raw_weights / sum(raw_weights) * weight_cfg$normalize_sum_to
  } else {
    weights <- raw_weights
  }
  
  labels <- generate_group_labels(n_groups, cfg$model$groups$labels)
  names(weights) <- labels
  return(weights)
}

expand_n_per_group <- function(n_per_group, n_groups) {
  if (length(n_per_group) == 1) return(rep(n_per_group, n_groups))
  if (length(n_per_group) < n_groups) return(rep_len(n_per_group, n_groups))
  if (length(n_per_group) == n_groups) return(n_per_group)
  stop("n_per_group must be length 1, a divisor of n_groups, or exactly n_groups.")
}

# --- Master function ---
generate_true_parameters <- function(cfg) {
  n_groups <- cfg$model$groups$n_groups
  labels <- generate_group_labels(n_groups, cfg$model$groups$labels)
  covariates <- cfg$model$covariates
  coef_cfg <- cfg$model$true_params$Beta_0$coefficients
  
  # --- Intercepts ---
  intercept_cfg <- cfg$model$true_params$Beta_0$intercepts
  dist_fun <- match.fun(intercept_cfg$distribution)
  intercepts <- do.call(dist_fun, c(list(n = n_groups), as.list(intercept_cfg$args)))
  
  if (!is.null(intercept_cfg$scale) && intercept_cfg$scale == "log") {
    theta_0 <- exp(intercepts)
  } else {
    theta_0 <- intercepts
    intercepts <- log(intercepts)
  }
  names(theta_0) <- labels
  
  # --- Stack coefficients ---
  beta_stack <- intercepts
  
  # Group-specific slopes
  for (cov in covariates) {
    cov_name <- cov$name
    if (!is.null(cov$varying_slope) && cov$varying_slope) {
      if (!is.null(coef_cfg[[cov_name]]$mode) && coef_cfg[[cov_name]]$mode == "fixed") {
        slopes <- rep(coef_cfg[[cov_name]]$args, n_groups)
      } else if (!is.null(coef_cfg[[cov_name]]$mode) && coef_cfg[[cov_name]]$mode == "sample") {
        dist_fun <- match.fun(coef_cfg[[cov_name]]$distribution)
        slopes <- do.call(dist_fun, c(list(n = n_groups), as.list(coef_cfg[[cov_name]]$args)))
      } else {
        stop(sprintf("Unknown or missing 'mode' for coefficient %s", cov_name))
      }
      beta_stack <- c(beta_stack, slopes)
    }
  }
  
  # Shared slopes
  for (cov in covariates) {
    cov_name <- cov$name
    if (!is.null(cov$varying_slope) && !cov$varying_slope) {
      if (!is.null(coef_cfg[[cov_name]]$mode) && coef_cfg[[cov_name]]$mode == "fixed") {
        slopes <- coef_cfg[[cov_name]]$args
      } else if (!is.null(coef_cfg[[cov_name]]$mode) && coef_cfg[[cov_name]]$mode == "sample") {
        dist_fun <- match.fun(coef_cfg[[cov_name]]$distribution)
        slopes <- do.call(dist_fun, c(list(n = 1), as.list(coef_cfg[[cov_name]]$args)))
      } else {
        stop(sprintf("Unknown or missing 'mode' for coefficient %s", cov_name))
      }
      beta_stack <- c(beta_stack, slopes)
    }
  }
  
  # --- Convert to matrix and add row names ---
  Beta_0 <- matrix(beta_stack, ncol = 1)
  rownames_vec <- c()
  # Intercepts
  rownames_vec <- c(rownames_vec, labels)
  # Group-specific slopes
  for (cov in covariates) {
    if (!is.null(cov$varying_slope) && cov$varying_slope) {
      rownames_vec <- c(rownames_vec, paste0(cov$name, "_", labels))
    }
  }
  # Shared slopes
  for (cov in covariates) {
    if (!is.null(cov$varying_slope) && !cov$varying_slope) {
      rownames_vec <- c(rownames_vec, cov$name)
    }
  }
  rownames(Beta_0) <- rownames_vec
  
  # --- Weights and n_per_group ---
  weights <- generate_group_weights(cfg)
  n_per_group <- expand_n_per_group(cfg$model$groups$n_per_group$args, n_groups)
  n_per_group <- setNames(n_per_group, labels)
  
  # --- Output ---
  list(
    Beta_0      = Beta_0,   # matrix with row names, 1 column
    theta_0     = theta_0,  # exp(intercepts)
    weights     = weights,
    n_per_group = n_per_group
  )
}



