# applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/helpers/truth_utils.R

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

# Generate theta vector
generate_theta_values <- function(cfg) {
  n_groups <- cfg$groups$n_groups
  if (!is.null(cfg$theta$values)) {
    # Explicit theta vector
    validate_theta_vector(cfg$theta$values, n_groups)
    return(cfg$theta$values)
  }
  if (is.null(cfg$theta$distribution) || is.null(cfg$theta$args)) {
    stop("Theta must have either 'values' or ('distribution' and 'args') in config.")
  }
  dist_fun <- match.fun(cfg$theta$distribution)
  dist_args <- cfg$theta$args
  return(do.call(dist_fun, c(list(n = n_groups), as.list(dist_args))))
}

# Validate theta vector
validate_theta_vector <- function(theta, expected_length) {
  if (length(theta) != expected_length) {
    stop(sprintf("θ vector length (%d) does not match n_groups (%d)",
                 length(theta), expected_length))
  }
  if (any(theta < 0)) stop("θ values must be non-negative")
  invisible(TRUE)
}

# Generate group weights (unchanged except label handling)
generate_group_weights <- function(cfg) {
  weight_cfg <- cfg$weights
  n_groups <- cfg$groups$n_groups
  
  if (is.null(weight_cfg$distribution) || is.null(weight_cfg$args)) {
    stop("Weight distribution and args must be specified in config.")
  }
  
  dist_fun <- match.fun(weight_cfg$distribution)
  dist_args <- weight_cfg$args
  raw_weights <- do.call(dist_fun, c(list(n = n_groups), as.list(dist_args)))
  
  if (!is.null(weight_cfg$normalize_mean_to)) {
    weights <- raw_weights / mean(raw_weights) * weight_cfg$normalize_mean_to
  } else if (!is.null(weight_cfg$normalize_sum_to)) {
    weights <- raw_weights / sum(raw_weights) * weight_cfg$normalize_sum_to
  } else {
    weights <- raw_weights
  }
  
  labels <- generate_group_labels(n_groups, cfg$groups$labels)
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
  n_groups <- cfg$groups$n_groups
  labels <- generate_group_labels(n_groups, cfg$groups$labels)
  
  theta_vals <- generate_theta_values(cfg)
  validate_theta_vector(theta_vals, n_groups)
  
  weights <- generate_group_weights(cfg)
  n_per_group <- expand_n_per_group(cfg$groups$n_per_group, n_groups)
  
  names(theta_vals) <- labels
  names(weights) <- labels
  names(n_per_group) <- labels
  
  list(
    theta_0     = theta_vals,
    weights     = weights,
    n_per_group = n_per_group
  )
}


