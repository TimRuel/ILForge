# applications/poisson/group_rates_weighted_sum/naive_group_rates/scripts/helpers/experiment_utils.R

# Miscellaneous -----------------------------------------------------------

fit_model <- function(data) {
  
  glm(Y ~ 0 + group, offset = log(t), family = poisson(), data = data)
}

get_lambda_MLE_from_model <- function(model) {
  
  model |> 
    coef() |> 
    exp() |> 
    set_names(model$xlevels$group)
}

get_lambda_MLE_from_data <- function(data) {
  
  data |> 
    reframe(lambda_MLE = Y / t) |> 
    pull(lambda_MLE) |> 
    set_names(data$group)
}

get_lambda_MLE <- function(x) {
  if (inherits(x, "glm")) {
    x |>
      coef() |>
      exp() |>
      set_names(x$xlevels$group)
    
  } else if (is.data.frame(x)) {
    x |>
      reframe(lambda_MLE = Y / t) |>
      pull(lambda_MLE) |>
      set_names(x$group)
    
  } else {
    stop("Input must be a glm model or a data frame with columns Y, t, and group.")
  }
}

get_psi_MLE <- function(lambda_MLE, weights) sum(lambda_MLE * weights)

get_psi_MLE_SE <- function(lambda_MLE, weights, exposure) sqrt(sum(weights^2 * lambda_MLE / exposure))

log_likelihood <- function(lambda, Y, t) sum(Y * log(lambda) - t * lambda)

likelihood <- function(lambda, Y, t) exp(log_likelihood(lambda, Y, t))

# Psi Grid ----------------------------------------------------------------

get_psi_grid <- function(data, weights, num_std_errors, step_size) {
  
  lambda_MLE <- get_lambda_MLE(data)
  
  psi_MLE <- get_psi_MLE(lambda_MLE, weights)
  
  psi_MLE_SE <- get_psi_MLE_SE(lambda_MLE, weights, data$t)
  
  psi_endpoints <- psi_MLE + num_std_errors * psi_MLE_SE
  
  psi_grid <- seq(psi_endpoints[1], psi_endpoints[2], step_size)
  
  return(psi_grid)
}

safe_max <- function(x, fallback) if (length(x) == 0) fallback else max(x)
safe_min <- function(x, fallback) if (length(x) == 0) fallback else min(x)

get_fine_psi_grid <- function(psi_grid,
                              branch_mode,
                              fine_step_size,
                              fine_window) {
  
  # Ensure psi_grid is sorted and unique
  psi_grid <- sort(unique(psi_grid))
  
  # Determine left and right limits using nearest psi_grid values
  lower_bound <- safe_max(psi_grid[psi_grid <= (branch_mode - fine_window)], min(psi_grid))
  upper_bound <- safe_min(psi_grid[psi_grid >= (branch_mode + fine_window)], max(psi_grid))
  
  # Build fine psi values from lower_bound to upper_bound
  fine_grid <- seq(from = lower_bound,
                   to = upper_bound,
                   by = fine_step_size) |> 
    round(8)
  
  # Combine with original grid, drop duplicates, and sort
  full_grid <- c(psi_grid, fine_grid) |> 
    round(8) |> 
    unique() |> 
    sort()
  
  return(full_grid)
}

# Beta_hat ----------------------------------------------------------------

safe_auglag <- possibly(nloptr::auglag)

get_lambda_hat <- function(lambda_hat_obj_fn, init_guess) {
  
  auglag(
    x0 = init_guess,
    fn = lambda_hat_obj_fn,
    localsolver = "SLSQP",
    deprecatedBehavior = FALSE)$par
}

# Branch Parameters ---------------------------------------------------------------

get_threshold <- function(lambda_MLE, Y, t, threshold_offset) ceiling(abs(log_likelihood(lambda_MLE, Y, t))) + threshold_offset

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

get_branch_mode <- function(data, omega_hat, optim_interval) {
  
  t <- data$t
  Y <- data$Y
  mu <- t * omega_hat
  
  branch <- function(psi) {
    
    lambda_hat_obj_fn <- function(lambda) {
      
      return(-log_likelihood(lambda, mu, t))
    }
    
    init_guess <- omega_hat
    
    lambda_hat <- get_lambda_hat(lambda_hat_obj_fn, init_guess)
    
    return(log_likelihood(lambda_hat, Y, t))
  }
  
  optimize(f = branch,
           interval = optim_interval,
           maximum = TRUE,
           tol = 0.1)$maximum
}

get_branch_params <- function(data, weights, psi_MLE, psi_grid) {
  
  t <- data$t
  Y <- data$Y
  optim_interval <- c(head(psi_grid, 1), tail(psi_grid, 1))

  while (TRUE) {
    
    omega_hat <- get_omega_hat(psi_MLE, weights)
    
    branch_mode <- get_branch_mode(data, omega_hat, optim_interval)
    
    if (between(branch_mode, optim_interval[1], optim_interval[2])) {
      
      mu <- t * omega_hat
      
      lambda_hat_obj_fn <- function(lambda) {
        
        return(-log_likelihood(lambda, mu, t))
      } 
      
      init_guess <- omega_hat
      
      lambda_hat <- get_lambda_hat(lambda_hat_obj_fn, init_guess)
      
      return(list(omega_hat = omega_hat,
                  lambda_hat = lambda_hat,
                  branch_mode = branch_mode))
    }
  }
}

# Integrated Log-Likelihood ---------------------------------------------------

run_branch_side <- function(direction, 
                            branch_params,
                            X,
                            Y_one_hot,
                            Z,
                            psi_grid, 
                            step_size,
                            fine_step_size,
                            fine_window) {
  
  phi <- branch_params$phi
  Beta_mode <- branch_params$Beta_mode
  branch_mode <- branch_params$branch_mode
  
  fine_grid <- get_fine_psi_grid(psi_grid = psi_grid,
                                 branch_mode = branch_mode,
                                 fine_step_size = fine_step_size,
                                 fine_window = fine_window)
  
  psi_working_grid <- if (direction == "left") {
    
    rev(fine_grid[fine_grid < branch_mode])
  } else {
    
    fine_grid[fine_grid > branch_mode]
  }
  
  init_guess <- Beta_mode
  mu <- get_E_Y(X, Z, phi, psi_hat)
  
  psi_vals <- c()
  log_L_tilde_vals <- c()
  
  for (psi in psi_working_grid) {
    
    Beta_hat_obj_fn <- function(Beta) {
      
      Beta <- Beta |> 
        matrix(nrow = nrow(phi),
               ncol = ncol(phi),
               byrow = TRUE)
      
      return(-log_likelihood(X, mu, Z, psi, Beta))
    }
    
    Beta_hat <- get_Beta_hat(Beta_hat_obj_fn, init_guess) |> 
      matrix(nrow = nrow(phi),
             ncol = ncol(phi),
             byrow = TRUE)
    
    log_L_tilde <- log_likelihood(X, Y_one_hot, Z, psi, Beta_hat)
    init_guess <- Beta_hat
    
    if (any(abs(psi_grid - psi) < 1e-8)) {
      psi_vals <- c(psi_vals, psi)
      log_L_tilde_vals <- c(log_L_tilde_vals, log_L_tilde)
    }
  }
  
  if (direction == "left") {
    
    psi_vals <- rev(psi_vals)
    log_L_tilde_vals <- rev(log_L_tilde_vals)
  }
  
  list(psi = psi_vals, Integrated = log_L_tilde_vals)
}

compute_IL_branch <- function(branch_params,
                              X,
                              Y_one_hot,
                              Z,
                              psi_grid,
                              step_size,
                              fine_step_size,
                              fine_window) {
  
  # Evaluate branch to the left of the peak
  left_result <- run_branch_side(
    direction          = "left",
    branch_params      = branch_params,
    X                  = X,
    Y_one_hot          = Y_one_hot,
    Z                  = Z,
    psi_grid           = psi_grid,
    step_size          = step_size,
    fine_step_size     = fine_step_size,
    fine_window        = fine_window
  )
  
  # Evaluate branch to the right of the peak
  right_result <- run_branch_side(
    direction          = "right",
    branch_params      = branch_params,
    X                  = X,
    Y_one_hot          = Y_one_hot,
    Z                  = Z,
    psi_grid           = psi_grid,
    step_size          = step_size,
    fine_step_size     = fine_step_size,
    fine_window        = fine_window
  )
  
  # Combine and sort
  out <- rbind(
    data.frame(psi = left_result$psi, Integrated = left_result$Integrated),
    data.frame(psi = right_result$psi, Integrated = right_result$Integrated)
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

get_integrated_LL <- function(config, model_df) {
  
  invisible(list2env(config$optimization_specs$IL, env = environment()))
  
  num_classes <- config$model$response$num_classes
  
  num_effects <- num_classes - 1
  
  reference_class <- config$model$response$reference_class
  
  formula <- config$model$formula
  
  X <- formula |> 
    str_remove("Y ") |> 
    as.formula() |> 
    model.matrix(data = model_df) |> 
    as.data.frame() |> 
    select(!Z) |> 
    as.matrix()
  
  Z <- model_df |> 
    pull(Z)
  
  Y_one_hot <- model.matrix(~ Y - 1, data = model_df)[, 1:num_effects]
  
  constraints <- get_constraints(model_df, num_effects, formula)
  
  formula <- as.formula(formula)
  
  model <- fit_model(model_df, formula, constraints, reference_class)
  
  Beta_MLE <- get_Beta_MLE(model)
  
  psi_hat <- get_psi_hat(model)
  
  psi_grid <- get_psi_grid(model, num_std_errors, step_size)
  
  optim_interval <- c(head(psi_grid, 1), tail(psi_grid, 1))
  
  Sigma_hat <- get_Sigma_hat(model)
  
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
    
    branch_params <- get_branch_params(
      X               = X,
      Y_one_hot       = Y_one_hot,
      Z               = Z,
      Beta_MLE        = Beta_MLE,
      Sigma_hat       = Sigma_hat,
      psi_hat         = psi_hat,
      optim_interval  = optim_interval
    )
    
    IL_branch <- compute_IL_branch(
      branch_params  = branch_params,
      X              = X,
      Y_one_hot      = Y_one_hot,
      Z              = Z,
      psi_grid       = psi_grid,
      step_size      = step_size,
      fine_step_size = fine_step_size,
      fine_window    = fine_window
    )
    
    list(IL_branch = IL_branch,
         branch_params = branch_params)
  }
  
  IL_branches <- lapply(result, `[[`, 1)
  branch_params <- lapply(result, `[[`, 2)
  
  log_L_bar <- get_log_L_bar(IL_branches)
  
  list(log_L_bar_df = log_L_bar$df,
       IL_branches = IL_branches,
       branches_matrix = log_L_bar$branches_matrix,
       branch_params = branch_params)
}

# Profile Log-Likelihood --------------------------------------------------

compute_profile_branch <- function(direction,
                                   X,
                                   Y_one_hot,
                                   Z,
                                   step_size,
                                   fine_step_size,
                                   fine_window,
                                   psi_hat,
                                   Beta_MLE,
                                   PLL_max,
                                   stopping_val) {
  
  mu <- get_E_Y(X, Z, Beta_MLE, psi_hat)
  
  psi_vals <- list()
  
  log_L_p_vals <- list()
  
  init_guess <- Beta_MLE |> 
    gdata::unmatrix(byrow = TRUE) |> 
    unname()
  
  log_L_p <- PLL_max
  
  if (direction == "left") {
    
    step_anchor <- floor(psi_hat / step_size) * step_size
    direction_sign <- -1
  } else {
    
    step_anchor <- ceiling(psi_hat / step_size) * step_size
    direction_sign <- 1
  }
  
  psi <- psi_hat + direction_sign * fine_step_size
  
  while (log_L_p > stopping_val) {
    
    Beta_hat_obj_fn <- function(Beta) {
      
      Beta <- Beta |> 
        matrix(nrow = nrow(Beta_MLE),
               ncol = ncol(Beta_MLE),
               byrow = TRUE)
      
      return(-log_likelihood(X, mu, Z, psi, Beta))
    }
    
    Beta_hat <- get_Beta_hat(Beta_hat_obj_fn, init_guess) |> 
      matrix(nrow = nrow(Beta_MLE),
             ncol = ncol(Beta_MLE),
             byrow = TRUE)
    
    log_L_p <- log_likelihood(X, Y_one_hot, Z, psi, Beta_hat)
    
    init_guess <- Beta_hat
    
    psi_vals[[length(psi_vals) + 1]] <- psi
    
    log_L_p_vals[[length(log_L_p_vals) + 1]] <- log_L_p
    
    dist_from_peak <- abs(psi - psi_hat)
    use_fine_step <- dist_from_peak < fine_window || 
      (direction == "left" && psi > step_anchor) ||
      (direction == "right" && psi < step_anchor)
    current_step <- if (use_fine_step) fine_step_size else step_size
    psi <- psi + direction_sign * current_step
  }
  
  psi_vals <- unlist(psi_vals)
  log_L_p_vals <- unlist(log_L_p_vals)
  
  if (direction == "left") {
    psi_vals <- c(rev(psi_vals), psi_hat)
    log_L_p_vals <- c(rev(log_L_p_vals), PLL_max)
  }
  
  list(psi = psi_vals, Profile = log_L_p_vals)
}

get_profile_LL <- function(config, model_df) {
  
  invisible(list2env(config$optimization_specs$PL, environment()))
  
  num_classes <- config$model$response$num_classes
  
  num_effects <- num_classes - 1
  
  reference_class <- config$model$response$reference_class
  
  formula <- config$model$formula
  
  X <- formula |> 
    str_remove("Y ") |> 
    as.formula() |> 
    model.matrix(data = model_df) |> 
    as.data.frame() |> 
    select(!Z) |> 
    as.matrix()
  
  Z <- model_df |> 
    pull(Z)
  
  Y_one_hot <- model.matrix(~ Y - 1, data = model_df)[, 1:num_effects]
  
  constraints <- get_constraints(model_df, num_effects, formula)
  
  formula <- as.formula(formula)
  
  model <- fit_model(model_df, formula, constraints, reference_class)
  
  Beta_MLE <- get_Beta_MLE(model)
  
  psi_hat <- get_psi_hat(model)
  
  PLL_max <- logLik(model)
  alpha <- min(alpha_levels)
  crit <- qchisq(1 - alpha, df = 1) / 2
  stopping_val <- PLL_max - crit
  
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
      X              = X,
      Y_one_hot      = Y_one_hot,
      Z              = Z,
      step_size      = step_size,
      fine_step_size = fine_step_size,
      fine_window    = fine_window,
      psi_hat        = psi_hat,
      Beta_MLE       = Beta_MLE,
      PLL_max        = PLL_max,
      stopping_val   = stopping_val
    )
  }
  
  profile_LL <- do.call(rbind, lapply(result, function(entry) {
    data.frame(psi = entry$psi, Profile = entry$Profile)
  }))
  
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
  
  J <- config$model$response$num_classes
  
  conf_ints <- get_confidence_intervals(
    pseudolikelihoods = pseudolikelihoods,
    LL_df_long        = LL_df_long,
    MLE_data          = MLE_data,
    alpha_levels      = alpha_levels,
    J                 = J
  )
  
  return(list(MLE_data = MLE_data,
              conf_ints = conf_ints))
}
