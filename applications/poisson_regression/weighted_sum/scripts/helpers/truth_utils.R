# applications/poisson_regression/weighted_sum/scripts/helpers/truth_utils.R

generate_true_parameters <- function(config) {
  
  set.seed(config$seed)
  
  # --- Number of covariates (p) ---
  p <- config$design_matrix$n_covariates
  if (isTRUE(config$design_matrix$intercept)) p <- p + 1
  
  # --- Generate Beta_0 ---
  Beta_0 <- do.call(match.fun(config$beta$distribution), c(list(p), config$beta$args))
  
  # --- Generate weights w ---
  n <- config$design_matrix$n_obs
  w <- do.call(match.fun(config$weights$distribution), c(list(n), config$weights$args))
  
  # --- Handle weight normalization ---
  norm_sum_to <- config$weights$normalize_sum_to %||% NULL
  norm_mean_to <- config$weights$normalize_mean_to %||% NULL
  
  if (!is.null(norm_sum_to) && !is.null(norm_mean_to)) {
    
    stop("Specify only one of 'normalize_sum_to' or 'normalize_mean_to', not both.")
  }
  
  if (!is.null(norm_sum_to)) {
    
    w <- w * (norm_sum_to / sum(w))
  } else if (!is.null(norm_mean_to)) {
    
    w <- w * (norm_mean_to / mean(w))
  }
  
  return(list(Beta_0 = Beta_0, weights = w))
}


