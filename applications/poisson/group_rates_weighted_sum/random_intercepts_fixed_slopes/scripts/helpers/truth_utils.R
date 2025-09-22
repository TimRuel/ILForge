# applications/poisson/group_rates_weighted_sum/random_intercepts_fixed_slopes/scripts/helpers/truth_utils.R

# Helper to generate group labels: A, B, ..., Z, AA, BB, ...
# Generate Excel-style column labels: A, B, ..., Z, AA, AB, ..., ZZ, AAA, etc.
make_group_labels <- function(G) {
  labels <- character(G)
  
  for (i in seq_len(G)) {
    n <- i
    s <- ""
    while (n > 0) {
      r <- (n - 1) %% 26
      s <- paste0(LETTERS[r + 1], s)
      n <- (n - 1) %/% 26
    }
    labels[i] <- s
  }
  
  labels
}

# Expand groups from config into a named vector of group sizes
expand_groups <- function(groups_config) {
  
  G <- groups_config$num_groups
  gen_type <- groups_config$n_per_group$gen_type
  
  if (gen_type == "fixed") {
    vals <- groups_config$n_per_group$fixed$values
    
    if (length(vals) == 1L) {
      n_per_group <- rep(vals, G)
    } else if (length(vals) == G) {
      n_per_group <- vals
    } else {
      stop("For fixed group sizes, 'values' must have length 1 or num_groups = ", G, ".")
    }
    
  } else if (gen_type == "random") {
    dist_config <- groups_config$n_per_group$random$distribution
    n_per_group <- sample_from_config(dist_config, G)
  } else {
    stop("Unknown group generation type: ", gen_type)
  }
  
  # Validation: positive integers
  if (any(n_per_group <= 0 | n_per_group != as.integer(n_per_group))) {
    stop("All group sizes must be positive integers. Got: ", paste(n_per_group, collapse = ", "))
  }
  
  n_per_group <- as.integer(n_per_group)
  names(n_per_group) <- make_group_labels(G)
  n_per_group
}

# Generic sampler for config blocks like:
#   distribution: { name: rnorm, args: [0,1] }
sample_from_config <- function(dist_config, n = 1) {
  dist_fun <- match.fun(dist_config$name)
  args <- dist_config$args
  if (is.null(args)) args <- list()
  
  # Combine n with args; if args is a vector, convert to list
  if (is.vector(args) && !is.list(args)) args <- as.list(args)
  out <- do.call(dist_fun, c(list(n = n), args))
  
  return(out)
}

filter_by_greek <- function(mat, greek_symbol) {
  if (is.null(rownames(mat))) {
    stop("Input matrix must have rownames.")
  }
  
  # get rownames
  rn <- rownames(mat)
  
  # keep rows where rowname starts with the greek symbol
  matches <- startsWith(rn, greek_symbol)
  
  # extract as named vector
  res <- mat[matches, , drop = TRUE]
  names(res) <- rn[matches]
  res
}

# Generate group weights  (reads from config$model$weights)
generate_group_weights <- function(config) {
  
  weight_config <- config$model$weights
  n_per_group <- expand_groups(config$model$groups)
  G <- length(n_per_group)
  
  if (is.null(weight_config$distribution) || is.null(weight_config$distribution$name))
    stop("Weight distribution must be specified as list(distribution=list(name=..., args=[...])).")
  
  raw <- sample_from_config(weight_config$distribution, G)
  
  if (!is.null(weight_config$normalize_mean_to)) {
    weights <- raw / mean(raw) * weight_config$normalize_mean_to
  } else if (!is.null(weight_config$normalize_sum_to)) {
    weights <- raw / sum(raw) * weight_config$normalize_sum_to
  } else {
    weights <- raw
  }
  
  names(weights) <- names(n_per_group)
  weights
}

get_Beta_0 <- function(config) {
  
  n_per_group <- expand_groups(config$model$groups)
  group_labels <- names(n_per_group)
  G <- length(n_per_group)
  covs <- config$model$covariates
  homo_covs <- Filter(\(c) c$type == "homogeneous", covs)
  hetero_covs <- Filter(\(c) c$type == "heterogeneous", covs)
  
  # ----- Generate coefficients for Beta_0 -----
  # Order: baseline intercept (α_0), homogeneous slopes (γ), heterogeneous slopes (ζ)
  
  # --- 1. Baseline fixed intercept ---
  alpha_0 <- config$model$intercepts$baseline
  
  # --- 2. Homogeneous slopes ---
  gamma <- sapply(homo_covs, function(c) sample_from_config(c$coefficient$distribution, 1))
  
  # --- 3. Heterogeneous slopes ---
  zeta_g <- c()
  for (cov in hetero_covs) {
    zeta_g <- c(zeta_g, sample_from_config(cov$coefficient$distribution, G))
  }  
  
  # --- 4. Concatenate Beta vector ---
  Beta_vals <- c(alpha_0, gamma, zeta_g)
  
  # --- 5. Assign rownames with symbols and "_group" suffixes ---
  Beta_names <- c(
    paste0(config$model$intercepts$symbol, "_0"), # α_0
    sapply(homo_covs, function(c) c$coefficient$symbol), # γ
    unlist(lapply(hetero_covs, function(c) paste0(c$coefficient$symbol, "_", group_labels))) # ζ_g
  )
  
  # --- 6. Filter missing values and stack into a matrix ---
  Beta_0 <- Beta_vals |> 
    setNames(Beta_names) |> 
    as.matrix(ncol = 1)
  
  return(Beta_0)
}

compute_true_marginal_rates <- function(config, Beta_0, sigma_alpha) {
  
  n_per_group <- expand_groups(config$model$groups)
  group_labels <- names(n_per_group)
  G <- length(group_labels)
  n_mc <- as.numeric(config$model$evaluation$marginal$n_mc)
  n_per_group <- rep(n_mc, G)
  names(n_per_group) <- group_labels
  
  X_mc <- get_X(config, Beta_0, n_per_group)

  exp_eta_mc <- exp(X_mc %*% Beta_0) 
  colnames(exp_eta_mc) <- "exp_eta"
    
  theta_0 <- exp_eta_mc |> 
    as_tibble(rownames = "group") |> 
    group_by(group) |> 
    summarise(mean_exp = mean(exp_eta) * exp(sigma_alpha^2 / 2)) |> 
    deframe()

  return(theta_0)
}

generate_true_parameters <- function(config) {
  
  Beta_0 <- get_Beta_0(config)
  
  sigma_alpha <- config$model$intercepts$distribution$args[2]

  weights <- generate_group_weights(config)
  
  theta_0 <- compute_true_marginal_rates(config, Beta_0, sigma_alpha)
  
  n_per_group <- expand_groups(config$model$groups)
  
  list(
    Beta_0 = Beta_0,  
    sigma_alpha = sigma_alpha, 
    theta_0 = theta_0, 
    weights = weights,  
    n_per_group = n_per_group
  )
}
