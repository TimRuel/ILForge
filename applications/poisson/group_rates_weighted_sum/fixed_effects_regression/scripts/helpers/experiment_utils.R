# applications/poisson/group_rates_weighted_sum/fixed_effects_regression/scripts/helpers/experiment_utils.R

# Miscellaneous -----------------------------------------------------------

make_formula <- function(response, covar_pattern, extras = NULL, offset = NULL, data) {
  
  covars <- grep(covar_pattern, names(data), value = TRUE)
  rhs <- paste(c(extras, covars), collapse = " + ")
  f <- as.formula(paste(response, "~", rhs))
  if (!is.null(offset)) {
    attr(f, "offset") <- offset
  }
  f
}

fit_model <- function(data) {
  
  formula <- make_formula("Y", "^X", extras = "0 + group", data = data)
  
  glm(formula, offset = log(t), family = poisson(), data = data)
}

get_Beta_MLE <- function(model) {
  
  group_labels <- model$xlevels$group
  covariate_labels <- model$terms |> 
    attr("term.labels") |> 
    setdiff("group")
  
  model |> 
    coef() |> 
    set_names(c(group_labels, covariate_labels))
}

get_theta <- function(Beta, G) exp(Beta[1:G])

get_psi <- function(Beta, weights, G) {
  
  theta <- get_theta(Beta, G)
  
  sum(theta * weights)
}

# In make_formula, incorporate offset properly
make_formula <- function(response, covar_pattern, extras = NULL, offset = NULL, data) {
  covars <- grep(covar_pattern, names(data), value = TRUE)
  rhs <- paste(c(extras, covars), collapse = " + ")
  if (!is.null(offset)) rhs <- paste(rhs, "+ offset(offset)")
  as.formula(paste(response, "~", rhs))
}

# General delta-method version of SE(psi_hat)
get_psi_MLE_SE <- function(model, weights, X_ref = NULL) {
  # X_ref: design matrix where psi_hat is evaluated (e.g., reference covariates)
  Beta_cov <- vcov(model)
  Beta_MLE <- get_Beta_MLE(model)
  
  if (is.null(X_ref)) {
    # naive: only intercepts
    G <- nlevels(model$data$group)
    grad <- c(weights * exp(Beta_MLE[1:G]), rep(0, length(Beta_MLE) - G))
  } else {
    eta <- X_ref %*% Beta_MLE
    grad <- t(X_ref) %*% (weights * exp(eta))
  }
  
  sqrt(as.numeric(t(grad) %*% Beta_cov %*% grad))
}

get_eta <- function(Beta, X) X %*% Beta

get_mu <- function(eta, t) t * exp(eta)
  
log_likelihood <- function(Beta, X, Y, t) {
  
  eta <- get_eta(Beta, X)
  
  mu <- get_mu(eta, t)
  
  sum(Y * eta - mu)
}

likelihood <- function(Beta, X, Y, t) exp(log_likelihood(Beta, X, Y, t))

# Psi Grid ----------------------------------------------------------------

safe_max <- function(x, fallback) if (length(x) == 0) fallback else max(x)
safe_min <- function(x, fallback) if (length(x) == 0) fallback else min(x)

get_psi_grid <- function(model, 
                         weights, 
                         num_std_errors, 
                         step_size,
                         fine_step_size,
                         fine_window) {
  
  Beta_MLE <- get_Beta_MLE(model)
  
  theta_MLE <- get_theta(Beta_MLE, G)
  
  psi_MLE <- get_psi(theta_MLE, weights)
  
  psi_MLE_SE <- get_psi_MLE_SE(model, weights)
  
  psi_endpoints <- psi_MLE + num_std_errors * psi_MLE_SE
  
  coarse_psi_grid <- seq(psi_endpoints[1], psi_endpoints[2], step_size)
  
  lower_bound <- safe_max(coarse_psi_grid[coarse_psi_grid <= (psi_MLE - fine_window)], min(coarse_psi_grid))
  upper_bound <- safe_min(coarse_psi_grid[coarse_psi_grid >= (psi_MLE + fine_window)], max(coarse_psi_grid))
  
  fine_psi_grid <- seq(from = lower_bound,
                       to = upper_bound,
                       by = fine_step_size) |> 
    round(8)
  
  psi_grid <- c(coarse_psi_grid, fine_psi_grid, psi_MLE) |> 
    round(8) |> 
    unique() |> 
    sort()
  
  return(psi_grid)
}

# Beta_hat ----------------------------------------------------------------

safe_auglag <- possibly(auglag)

Beta_hat_con_fn_template <- function(weights, G, psi) {
  
  function(Beta) abs(get_psi(Beta, weights, G) - psi)
}

Beta_hat_obj_fn_template <- function(omega_hat, X, Y, t) {
  
  function(Beta) {
    
    E_Y <- t * exp(get_eta(omega_hat, X))
    
    -log_likelihood(Beta, X, E_Y, t)
  }
}

get_Beta_hat <- function(Beta_hat_obj_fn, Beta_hat_con, init_guess) {
  
  safe_auglag(
    x0 = init_guess,
    fn = Beta_hat_obj_fn,
    heq = Beta_hat_con_fn,
    localsolver = "SLSQP",
    deprecatedBehavior = FALSE)$par
}
  
# omega_hat ---------------------------------------------------------------

get_omega_hat <- function(psi_MLE, weights, n_covariates) {
  
  n_groups <- length(weights)
  
  constraints <- hitandrun::simplexConstraints(n_groups)
  constraints$constr[1,] <- weights
  constraints$rhs[[1]] <- psi_MLE
  
  omega_hat <- constraints |> 
    hitandrun::hitandrun(1) |>
    log() |> 
    c(rnorm(n_covariates))
  
  return(omega_hat)
}

# Integrated Log-Likelihood ---------------------------------------------------

run_branch_side <- function(direction, 
                            omega_hat,
                            data,
                            weights,
                            psi_MLE,
                            psi_grid) {
  
  psi_working_grid <- if (direction == "left") {
    
    rev(psi_grid[psi_grid <= psi_MLE])
  } else {
    
    psi_grid[psi_grid > psi_MLE]
  }
  
  Y <- data$Y
  
  t <- data$t
  
  log_L_tilde_vals <- c()
  
  init_guess <- 0
  
  for (psi in psi_working_grid) {
    
    lambda_con_fn <- lambda_con_fn_template(omega_hat, weights, t, psi)
    
    lambda <- get_lambda(lambda_con_fn, init_guess)
    
    theta_hat <- get_theta_hat(lambda, omega_hat, weights, t)
    
    log_L_tilde <- log_likelihood(theta_hat, Y, t)
    
    log_L_tilde_vals <- c(log_L_tilde_vals, log_L_tilde)
    
    init_guess <- lambda
  }
  
  if (direction == "left") {
    
    psi_working_grid <- rev(psi_working_grid)
    log_L_tilde_vals <- rev(log_L_tilde_vals)
  }
  
  list(psi = psi_working_grid, Integrated = log_L_tilde_vals)
}

compute_IL_branch <- function(omega_hat,
                              data,
                              weights,
                              psi_MLE,
                              psi_grid) {
  
  # Evaluate branch to the left of the peak
  left_branch <- run_branch_side(
    direction = "left", 
    omega_hat = omega_hat,
    data      = data,
    weights   = weights,
    psi_MLE   = psi_MLE,
    psi_grid  = psi_grid
  )
  
  # Evaluate branch to the right of the peak
  right_branch <- run_branch_side(
    direction = "right", 
    omega_hat = omega_hat,
    data      = data,
    weights   = weights,
    psi_MLE   = psi_MLE,
    psi_grid  = psi_grid
  )
  
  # Combine and sort
  out <- rbind(
    data.frame(psi = left_branch$psi, Integrated = left_branch$Integrated),
    data.frame(psi = right_branch$psi, Integrated = right_branch$Integrated)
  )
  
  out <- out[order(out$psi), , drop = FALSE]
  rownames(out) <- NULL
  out
}

get_log_L_bar <- function(IL_branches) {
  
  merged_df <- reduce(IL_branches, full_join, by = "psi")
  
  branches_matrix <- merged_df[, -1] |>
    as.matrix() |>
    t() |>
    unname()
  
  colnames(branches_matrix) <- merged_df$psi
  
  log_R <- branches_matrix |>
    nrow() |>
    log()
  
  log_L_bar <- matrixStats::colLogSumExps(branches_matrix, na.rm = TRUE) - log_R
  
  log_L_bar_df <- data.frame(psi = merged_df$psi,
                             Integrated = log_L_bar)
  
  return(list(df = log_L_bar_df,
              branches_matrix = branches_matrix))
}

get_integrated_LL <- function(config, data, weights) {
  
  invisible(list2env(config$optimization_specs$IL, env = environment()))
  
  psi_grid <- get_psi_grid(data, 
                           weights, 
                           num_std_errors, 
                           step_size,
                           fine_step_size,
                           fine_window)
  
  theta_MLE <- get_theta_MLE(data)
  
  psi_MLE <- get_psi_MLE(theta_MLE, weights)
  
  num_branches <- chunk_size * num_workers
  
  result <- foreach(
    
    i = 1:num_branches,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = num_branches,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size,
                           packages = "nloptr")
  ) %dofuture% {
    
    omega_hat <- get_omega_hat(psi_MLE, weights)
    
    IL_branch <- compute_IL_branch(
      omega_hat = omega_hat,
      data      = data,
      weights   = weights,
      psi_MLE   = psi_MLE,
      psi_grid  = psi_grid
    )
    
    list(IL_branch = IL_branch,
         omega_hat = omega_hat)
  }
  
  IL_branches <- lapply(result, `[[`, 1)
  omega_hat <- lapply(result, `[[`, 2)
  
  log_L_bar <- get_log_L_bar(IL_branches)
  
  list(log_L_bar_df = log_L_bar$df,
       IL_branches  = IL_branches,
       omega_hat    = omega_hat)
}

# Profile Log-Likelihood --------------------------------------------------

compute_profile_branch <- function(direction,
                                   data,
                                   weights,
                                   step_size,
                                   fine_step_size,
                                   fine_window,
                                   alpha) {
  
  Y <- data$Y
  
  t <- data$t
  
  theta_MLE <- get_theta_MLE(data)
  
  log_L_p_max <- log_likelihood(theta_MLE, Y, t)
  
  crit <- qchisq(1 - alpha, df = 1) / 2
  
  stopping_val <- log_L_p_max - crit
  
  log_L_p <- log_L_p_max
  
  psi_MLE <- get_psi_MLE(theta_MLE, weights)
  
  psi_vals <- c()
  
  log_L_p_vals <- c()
  
  if (direction == "left") {
    
    step_anchor <- floor(psi_MLE / step_size) * step_size
    direction_sign <- -1
  } else {
    
    step_anchor <- ceiling(psi_MLE / step_size) * step_size
    direction_sign <- 1
  }
  
  psi <- psi_MLE + direction_sign * fine_step_size
  
  init_guess <- 0
  
  while (log_L_p >= stopping_val) {
    
    lambda_con_fn <- lambda_con_fn_template(theta_MLE, weights, t, psi)
    
    lambda <- get_lambda(lambda_con_fn, init_guess)
    
    theta_hat <- get_theta_hat(lambda, theta_MLE, weights, t)
    
    log_L_p <- log_likelihood(theta_hat, Y, t)
    
    log_L_p_vals <- c(log_L_p_vals, log_L_p)
    
    psi_vals <- c(psi_vals, psi)
    
    init_guess <- lambda
    
    dist_from_peak <- abs(psi - psi_MLE)
    
    use_fine_step <- dist_from_peak < fine_window || 
      (direction == "left" && psi > step_anchor) ||
      (direction == "right" && psi < step_anchor)
    
    current_step <- if (use_fine_step) fine_step_size else step_size
    
    psi <- psi + direction_sign * current_step
  }
  
  if (direction == "left") {
    
    psi_vals <- c(rev(psi_vals), psi_MLE)
    log_L_p_vals <- c(rev(log_L_p_vals), log_L_p_max)
  }
  
  list(psi = psi_vals, Profile = log_L_p_vals)
}

get_profile_LL <- function(config, data, weights) {
  
  invisible(list2env(config$optimization_specs$PL, environment()))
  
  alpha <- min(alpha_levels)
  
  result <- foreach(
    dir = c("left", "right"),
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = 2,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = 1,
                           packages ="nloptr")
  ) %dofuture% {
    
    compute_profile_branch(
      direction      = dir,
      data           = data,
      weights        = weights,
      step_size      = step_size,
      fine_step_size = fine_step_size,
      fine_window    = fine_window,
      alpha          = alpha
    )
  }
  
  profile_LL <- rbind |> 
    do.call(lapply(result, 
                   \(entry) {
                     data.frame(psi = entry$psi, Profile = entry$Profile)
                     }
                   )
            )
  
  return(profile_LL)
}

get_report_objects <- function(iter_dir) {
  
  results_dir <- here(iter_dir, "results")
  config_path <- here(iter_dir, "config_snapshot.yml")
  config <- read_yaml(config_path)
  
  integrated_LL <- readRDS(here(results_dir, "integrated_LL.rds"))
  profile_LL <- readRDS(here(results_dir, "profile_LL.rds"))
  
  LL_df <- integrated_LL$log_L_bar_df |>
    merge(profile_LL, all = TRUE)
  
  LL_df_long <- get_LL_df_long(LL_df)
  
  spline_models <- get_spline_models(LL_df_long)
  
  MLE_data <- get_MLE_data(spline_models, LL_df_long)
  
  pseudolikelihoods <- get_pseudolikelihoods(spline_models, MLE_data)
  
  alpha_levels <- config$optimization_specs$PL$alpha_levels
  
  conf_ints <- get_confidence_intervals(
    pseudolikelihoods = pseudolikelihoods,
    LL_df_long        = LL_df_long,
    MLE_data          = MLE_data,
    alpha_levels      = alpha_levels
  )
  
  return(list(MLE_data = MLE_data,
              conf_ints = conf_ints))
}
