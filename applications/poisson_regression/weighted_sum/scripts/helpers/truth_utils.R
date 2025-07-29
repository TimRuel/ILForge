# applications/poisson_regression/weighted_sum/scripts/helpers/truth_utils.R

# ------------------------------
# Helper: Create coefficient names
# ------------------------------
get_coef_names <- function(covariates, has_intercept) {
  
  names <- c()
  if (has_intercept) names <- c(names, "intercept")
  if (length(covariates) > 0) names <- c(names, paste0("X", seq_along(covariates)))
  names
}

# ------------------------------
# Helper: Generate true values of regression coefficients
# ------------------------------
generate_Beta_0 <- function(config, p, coef_names) {
  
  G <- config$group_structure$n_groups
  pooled <- config$group_structure$pooled
  
  if (pooled) {
    beta <- do.call(match.fun(config$beta$distribution), c(list(n = p), config$beta$args))
    names(beta) <- coef_names
    return(beta)
  } else {
    beta <- replicate(
      G,
      do.call(match.fun(config$beta$distribution), c(list(n = p), config$beta$args)),
      simplify = "matrix"
    )
    rownames(beta) <- coef_names
    colnames(beta) <- paste0("group_", seq_len(G))
    return(beta)
  }
}

# ------------------------------
# Helper: Generate group-level covariates
# ------------------------------
generate_covariates <- function(covariate_specs, G, has_intercept) {
  
  X_list <- lapply(covariate_specs, function(spec) {
    dist_fn <- match.fun(spec$distribution)
    args <- spec$args
    do.call(dist_fn, c(list(G), args))
  })
  X <- do.call(cbind, X_list)
  colnames(X) <- sapply(covariate_specs, `[[`, "name")
  
  if (has_intercept) {
    X <- cbind(Intercept = 1, X)
  }
  return(X)
}

# ------------------------------
# Helper: Compute true values of Poisson group rates
# ------------------------------
get_lambda_0 <- function(X, Beta_0) {
  lambda_0 <- X %*% Beta_0 |> diag() |> exp()
  names(lambda_0) <- paste0("group_", seq_len(length(lambda_0)))
  return(lambda_0)
}

# ------------------------------
# Helper: Generate normalized weights
# ------------------------------
generate_weights <- function(config, G) {
  
  w <- do.call(match.fun(config$weights$distribution), c(list(n = G), config$weights$args))
  
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
  return(w)
}

# ------------------------------
# Helper: Compute true value of parameter of interest
# ------------------------------
get_psi_0 <- function(weights, lambda_0) sum(weights * lambda_0)

# ------------------------------
# Main wrapper function
# ------------------------------
generate_true_parameters <- function(config) {
  
  set.seed(config$seed)
  
  G <- config$group_structure$n_groups
  covariate_specs <- config$group_covariates$covariates
  has_intercept <- isTRUE(config$group_covariates$intercept)
  
  coef_names <- get_coef_names(covariate_specs, has_intercept)
  p <- length(coef_names)
  
  Beta_0 <- generate_Beta_0(config, p, coef_names)
  X <- generate_covariates(covariate_specs, G, has_intercept)
  lambda_0 <- get_lambda_0(X, Beta_0)
  w <- generate_weights(config, G)
  psi_0 <- get_psi_0(w, lambda_0)
  
  list(
    Beta_0 = Beta_0,
    covariates = X,
    weights = w,
    lambda_0 = lambda_0,
    psi_0 = psi_0
  )
}
