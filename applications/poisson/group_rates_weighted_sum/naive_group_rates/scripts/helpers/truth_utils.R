# applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/helpers/truth_utils.R

validate_lambda_vector <- function(lambda, expected_length) {
  if (length(lambda) != expected_length) {
    stop(sprintf("λ vector length (%d) does not match n_groups (%d)", length(lambda), expected_length))
  }
  if (any(lambda < 0)) stop("λ values must be non-negative")
  invisible(TRUE)
}

generate_group_weights <- function(cfg) {
  weight_cfg <- cfg$weights
  n_groups <- cfg$groups$n_groups
  
  # Validate required fields
  if (is.null(weight_cfg$distribution) || is.null(weight_cfg$args)) {
    stop("Weight distribution and args must be specified in config.")
  }
  
  dist_fun <- match.fun(weight_cfg$distribution)
  dist_args <- weight_cfg$args
  
  # Sample raw weights
  raw_weights <- do.call(dist_fun, c(list(n = n_groups), as.list(dist_args)))
  
  # Normalize
  if (!is.null(weight_cfg$normalize_mean_to)) {
    target_mean <- weight_cfg$normalize_mean_to
    weights <- raw_weights / mean(raw_weights) * target_mean
  } else if (!is.null(weight_cfg$normalize_sum_to)) {
    target_sum <- weight_cfg$normalize_sum_to
    weights <- raw_weights / sum(raw_weights) * target_sum
  } else {
    weights <- raw_weights
  }
  
  names(weights) <- cfg$groups$labels
  return(weights)
}

generate_true_parameters <- function(cfg) {
  lambda_vals <- cfg$lambda$values
  n_groups <- cfg$groups$n_groups
  
  validate_lambda_vector(lambda_vals, n_groups)
  
  weights <- generate_group_weights(cfg)
  
  names(lambda_vals) <- names(weights)
  
  list(
    lambda_0 = lambda_vals,
    weights = weights
  )
}
