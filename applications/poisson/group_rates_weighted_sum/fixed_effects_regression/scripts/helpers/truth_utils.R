# applications/poisson/group_rates_weighted_sum/fixed_effects_regression/scripts/helpers/truth_utils.R

# Generate group labels automatically
generate_group_labels <- function(n_groups, style = "auto_upper") {
  if (is.character(style) && length(style) == 1) {
    if (style == "auto_upper") {
      # e.g., A, B, ..., Z, AA, AB, ...
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
  if (length(n_per_group) == 1) {
    return(rep(n_per_group, n_groups))
  }
  if (length(n_per_group) < n_groups) {
    # Repeat pattern until length matches
    return(rep_len(n_per_group, n_groups))
  }
  if (length(n_per_group) == n_groups) {
    return(n_per_group)
  }
  stop("n_per_group must be length 1, a divisor of n_groups, or exactly n_groups.")
}

# Master function
generate_true_parameters <- function(cfg) {
  n_groups <- cfg$model$groups$n_groups
  labels <- generate_group_labels(n_groups, cfg$model$groups$labels)
  
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
  
  # --- Coefficients ---
  coef_cfg <- cfg$model$true_params$Beta_0$coefficients
  n_covariates <- length(cfg$model$covariates)
  
  if (!is.null(coef_cfg$mode) && coef_cfg$mode == "fixed") {
    # fixed vector applied to all groups
    coefficients <- coef_cfg$args
  } else if (!is.null(coef_cfg$mode) && coef_cfg$mode == "sample") {
    # draw coefficients from a distribution
    if (is.null(coef_cfg$distribution) || is.null(coef_cfg$args)) {
      stop("Coefficient sampling requires 'distribution' and 'args'")
    }
    dist_fun <- match.fun(coef_cfg$distribution)
    coefficients <- do.call(dist_fun, c(list(n = n_covariates), as.list(coef_cfg$args)))
  } else {
    stop("Unknown or missing 'mode' for coefficients")
  }
  
  names(coefficients) <- sapply(cfg$model$covariates, `[[`, "name")
  
  # --- Combine intercept + coefficients into Beta_0 matrix ---
  Beta_0 <- matrix(c(intercepts, coefficients))
  rownames(Beta_0) <- c(labels, names(coefficients))
  
  # --- Weights and n_per_group ---
  weights <- generate_group_weights(cfg)
  n_per_group <- expand_n_per_group(cfg$model$groups$n_per_group$args, n_groups)
  n_per_group <- setNames(n_per_group, labels)
  
  # --- Output ---
  list(
    Beta_0      = Beta_0,    # matrix: (# groups + # covariates) x 1
    theta_0     = theta_0,   # exp(intercepts) for Poisson rates
    weights     = weights,
    n_per_group = n_per_group
  )
}

