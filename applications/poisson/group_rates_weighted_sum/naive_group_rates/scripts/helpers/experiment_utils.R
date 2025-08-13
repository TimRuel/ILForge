# applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/helpers/experiment_utils.R

# Miscellaneous -----------------------------------------------------------

fit_model <- function(data) {
  
  glm(Y ~ 0 + group, offset = log(t), family = poisson(), data = data)
}

get_theta_MLE_from_model <- function(model) {
  
  model |> 
    coef() |> 
    exp() |> 
    set_names(model$xlevels$group)
}

get_theta_MLE_from_data <- function(data) {
  
  data |> 
    reframe(theta_MLE = Y / t) |> 
    pull(theta_MLE) |> 
    set_names(data$group)
}

get_theta_MLE <- function(x) {
  if (inherits(x, "glm")) {
    x |>
      coef() |>
      exp() |>
      set_names(x$xlevels$group)
    
  } else if (is.data.frame(x)) {
    x |>
      reframe(theta_MLE = Y / t) |>
      pull(theta_MLE) |>
      set_names(x$group)
    
  } else {
    stop("Input must be a glm model or a data frame with columns Y, t, and group.")
  }
}

get_psi_MLE <- function(theta_MLE, weights) sum(theta_MLE * weights)

get_psi_MLE_SE <- function(theta_MLE, weights, exposure) sqrt(sum(weights^2 * theta_MLE / exposure))

log_likelihood <- function(theta, Y, t) sum(dpois(Y, t * theta, log = TRUE))

likelihood <- function(theta, Y, t) sum(dpois(Y, t * theta))

# Psi Grid ----------------------------------------------------------------

safe_max <- function(x, fallback) if (length(x) == 0) fallback else max(x)
safe_min <- function(x, fallback) if (length(x) == 0) fallback else min(x)

get_psi_grid <- function(data, 
                         weights, 
                         num_std_errors, 
                         step_size,
                         fine_step_size,
                         fine_window) {
  
  theta_MLE <- get_theta_MLE(data)
  
  psi_MLE <- get_psi_MLE(theta_MLE, weights)
  
  psi_MLE_SE <- get_psi_MLE_SE(theta_MLE, weights, data$t)
  
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

# theta_hat ----------------------------------------------------------------

get_theta_hat <- function(lambda, omega_hat, weights, t) (t * omega_hat) / (t - lambda * weights)

lambda_con_fn_template <- function(omega_hat, weights, t, psi) {
  
  function(lambda) {
    
    theta_hat <- get_theta_hat(lambda, omega_hat, weights, t)
    
    sum(weights * theta_hat) - psi
  }
}

safe_auglag <- possibly(auglag)

get_lambda <- function(lambda_con_fn, init_guess) {
  
  safe_auglag(
    x0 = init_guess,
    fn = function(lambda) 0,
    heq = lambda_con_fn,
    localsolver = "SLSQP",
    deprecatedBehavior = FALSE)$par
}

# omega_hat ---------------------------------------------------------------

get_threshold <- function(theta_MLE, Y, t, threshold_offset) ceiling(abs(log_likelihood(theta_MLE, Y, t))) + threshold_offset

get_omega_hat <- function(psi_MLE, weights) {
  
  n_groups <- length(weights)
  
  constraints <- hitandrun::simplexConstraints(n_groups)
  constraints$constr[1,] <- weights
  constraints$rhs[[1]] <- psi_MLE
  
  omega_hat <- constraints |> 
    hitandrun::hitandrun(1) |>
    c()
  
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
