# applications/poisson/group_rates_weighted_sum/fixed_effects_regression/scripts/helpers/experiment_utils.R

make_formula <- function(config) {
  # Build RHS parts
  rhs_parts <- c("0", "group")
  
  covs <- config$model$covariates
  homo_covs <- Filter(\(c) c$type == "homogeneous", covs)
  hetero_covs <- Filter(\(c) c$type == "heterogeneous", covs)
  
  # Homogeneous covariates (same slope across groups)
  if (length(homo_covs) > 0) {
    homo_covs <- sapply(homo_covs, \(x) x$variable$symbol)
    rhs_parts <- c(rhs_parts, homo_covs)
  }
  
  # Heterogeneous covariates (different slope across groups, can be interpreted as interaction with group)
  if (length(hetero_covs) > 0) {
    hetero_covs <- sapply(hetero_covs, \(x) x$variable$symbol)
    rhs_parts <- c(rhs_parts, paste0(hetero_covs, ":group"))
  }
  
  # Add group main effect (intercepts per group)
  rhs <- paste(rhs_parts, collapse = " + ")
  
  # Full formula
  f <- as.formula(paste("Y ~", rhs))
  
  return(f)
}

fit_model <- function(config, data) {
  
  formula <- make_formula(config)
  
  glm(formula, offset = log(t), family = poisson(), data = data)
}

rename_coefs <- function(coef_names) {
  # Separate intercepts, main effects, and interactions
  intercepts   <- coef_names[!grepl("X", coef_names)]              # group names
  main_effects <- coef_names[grepl("^X\\d+$", coef_names)]         # X1, X2, ...
  interactions <- coef_names[grepl(":", coef_names)]               # group:X2, ...
  
  # Figure out which X's are fixed vs varying
  # Fixed slopes: main effects only
  fixed_Xs <- gsub("^X", "", main_effects)
  # Varying slopes: from interactions
  varying_Xs <- unique(gsub(".*:X", "", interactions))
  
  # Replacements
  new_names <- character(length(coef_names))
  
  # Intercepts → α_group
  new_names[match(intercepts, coef_names)] <- paste0("α_", sub("^group", "", intercepts))
  
  # Fixed slopes → γ#
  new_names[match(main_effects, coef_names)] <- 
    paste0("γ", gsub("^X", "", main_effects))
  
  # Varying slopes → ζ#_group
  for (int in interactions) {
    parts <- strsplit(int, ":")[[1]]   # e.g., c("groupA", "X2")
    group <- sub("^group", "", parts[1])
    xnum  <- sub("^X", "", parts[2])
    new_names[match(int, coef_names)] <- paste0("ζ", xnum, "_", group)
  }
  
  return(new_names)
}

get_Beta_MLE <- function(model) {
  
  coefs <- coef(model)
  Beta_MLE <- as.matrix(coefs, ncol = 1)  # stack as single-column matrix
  
  rownames(Beta_MLE) <- rename_coefs(names(coefs))
  
  Beta_MLE
}

get_eta <- function(Beta, X) X %*% Beta

get_theta <- function(Beta, X) exp(get_eta(Beta, X))

get_theta_g <- function(Beta, X) {
  
  data.frame(theta = get_theta(Beta, X),
             group = factor(rownames(X))) |> 
    group_by(group) |> 
    summarise(theta = mean(theta)) |> 
    deframe()
}

get_mu <- function(eta, t) t * exp(eta)

get_psi <- function(Beta, X, weights) {
  
  theta <- get_theta_g(Beta, X)
  
  sum(theta * weights)
}

# General delta-method version of SE(psi_hat)
get_psi_MLE_SE <- function(model, weights, X) {
  # Extract Beta_MLE as numeric vector
  Beta_MLE <- drop(get_Beta_MLE(model))
  Beta_cov <- vcov(model)
  
  # Group membership for each observation
  groups <- model$data$group
  group_labels <- levels(groups)
  
  # Ensure weights are a named vector (one per group)
  if (is.null(names(weights))) {
    if (length(weights) != length(group_labels)) {
      stop("If 'weights' is unnamed, it must have length equal to number of groups.")
    }
    weights <- setNames(weights, group_labels)
  }
  
  # Linear predictors and exp
  eta <- drop(get_eta(Beta_MLE, X))
  exp_eta <- exp(eta)
  
  # Group sizes
  n_g <- table(groups)
  
  # Per-observation effective weights = w_g / n_g
  obs_w <- weights[as.character(groups)] / n_g[as.character(groups)]
  
  # Gradient: sum_i obs_w[i] * exp(eta_i) * x_i
  grad <- as.numeric(t(X) %*% (obs_w * exp_eta))
  names(grad) <- colnames(X)
  
  # Delta-method SE
  se <- sqrt(as.numeric(t(grad) %*% Beta_cov %*% grad))
  
  return(se)
}
  
log_likelihood <- function(Beta, X, Y, t) {
  
  eta <- get_eta(Beta, X)
  
  mu <- get_mu(eta, t)
  
  sum(Y * (log(t) + eta) - mu - lgamma(Y+1))
}

likelihood <- function(Beta, X, Y, t) exp(log_likelihood(Beta, X, Y, t))

# Psi Grid ----------------------------------------------------------------

safe_max <- function(x, fallback) if (length(x) == 0) fallback else max(x)
safe_min <- function(x, fallback) if (length(x) == 0) fallback else min(x)

get_psi_grid <- function(model,
                         X,
                         weights, 
                         num_std_errors, 
                         step_size,
                         fine_step_size,
                         fine_window) {
  
  Beta_MLE <- get_Beta_MLE(model)
  
  psi_MLE <- get_psi(Beta_MLE, X, weights)
  
  psi_MLE_SE <- get_psi_MLE_SE(model, weights, X)
  
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

get_search_interval <- function(model,
                                X,
                                weights, 
                                num_std_errors) {
  
  Beta_MLE <- get_Beta_MLE(model)
  
  psi_MLE <- get_psi(Beta_MLE, X, weights)
  
  psi_MLE_SE <- get_psi_MLE_SE(model, weights, X)
  
  psi_MLE + c(-1, 1) * num_std_errors * psi_MLE_SE
}

safe_auglag <- possibly(nloptr::auglag)

# Beta_hat ----------------------------------------------------------------

Beta_hat_con_fn_template <- function(X, weights, psi) {
  
  function(Beta) get_psi(Beta, X, weights) - psi
}

Beta_hat_con_fn_jac_template <- function(X, weights) {
  
  groups <- rownames(X)
  
  n_g <- table(groups)
  
  W <- as.numeric(weights[as.character(groups)] / n_g[as.character(groups)])
  
  function(Beta) {
    
    theta_Beta <- get_theta(Beta, X)
    
    as.numeric(t(X) %*% (W * theta_Beta))
  }
}

Beta_hat_obj_fn_template <- function(omega_hat, X, Y, t) {
  
  E_Y <- t * exp(get_eta(omega_hat, X))
  
  function(Beta) {
    
    -log_likelihood(Beta, X, E_Y, t)
  }
}

Beta_hat_obj_fn_gr_template <- function(omega_hat, X, Y, t) {
  
  theta_omega_hat <- get_theta(omega_hat, X)
  
  function(Beta) {
    
    theta_Beta <- get_theta(Beta, X)
    
    as.numeric(t(X) %*% (t * (theta_Beta - theta_omega_hat)))
  }
}

get_Beta_hat <- function(Beta_hat_obj_fn,
                         Beta_hat_obj_fn_gr,
                         Beta_hat_con_fn,
                         Beta_hat_con_fn_jac,
                         init_guess) {
  
  nloptr::auglag(
    x0 = init_guess,
    fn = Beta_hat_obj_fn,
    gr = Beta_hat_obj_fn_gr,
    heq = Beta_hat_con_fn,
    heqjac = Beta_hat_con_fn_jac,
    localsolver = "SLSQP",
    control = list(xtol_rel = 1e-8, maxeval = 1000),
    deprecatedBehavior = FALSE)$par
}
  
# omega_hat ---------------------------------------------------------------

omega_hat_con_fn_template <- function(X, weights, psi_MLE) {
  
  function(omega_hat) get_psi(omega_hat, X, weights) - psi_MLE
}

get_omega_hat <- function(omega_hat_con_fn, init_guess) {
  
  safe_auglag(
    x0 = init_guess,
    fn = function(omega_hat) 0,
    heq = omega_hat_con_fn,
    localsolver = "SLSQP",
    control = list(xtol_rel = 1e-8, maxeval = 1000),
    deprecatedBehavior = FALSE)$par
}

get_branch_params <- function(omega_hat, X, Y, t, weights, search_interval) {
  
  Beta_hat_obj_fn <- Beta_hat_obj_fn_template(omega_hat, X, Y, t)
  Beta_hat_obj_fn_gr <- Beta_hat_obj_fn_gr_template(omega_hat, X, Y, t)
  
  Beta_hat_con_fn_jac <- Beta_hat_con_fn_jac_template(X, weights)
  
  branch <- function(psi) {
    
    Beta_hat_con_fn <- Beta_hat_con_fn_template(X, weights, psi)
    
    Beta_hat <- get_Beta_hat(Beta_hat_obj_fn,
                             Beta_hat_obj_fn_gr,
                             Beta_hat_con_fn,
                             Beta_hat_con_fn_jac,
                             omega_hat)

    return(log_likelihood(Beta_hat, X, Y, t))
  }
  
  opt <- optimize(branch,
                  interval = search_interval,
                  maximum = TRUE,
                  tol = 0.1)
  
  psi_mode <- opt$maximum
  Beta_hat_mode <- get_Beta_hat(
    Beta_hat_obj_fn,
    Beta_hat_obj_fn_gr,
    Beta_hat_con_fn_template(X, weights, psi_mode),
    Beta_hat_con_fn_jac,
    omega_hat
  )
  
  list(psi = psi_mode,
       Beta_hat = Beta_hat_mode,
       omega_hat = omega_hat)
}

# Integrated Log-Likelihood ---------------------------------------------------

run_branch_side <- function(X,
                            Y,
                            t,
                            weights,
                            omega_hat,
                            init_guess, 
                            psi_grid) {
  
  Beta_hat_obj_fn <- Beta_hat_obj_fn_template(omega_hat, X, Y, t)
  Beta_hat_obj_fn_gr <- Beta_hat_obj_fn_gr_template(omega_hat, X, Y, t)
  
  Beta_hat_con_fn_jac <- Beta_hat_con_fn_jac_template(X, weights)
  
  log_L_b_vals <- c()
  
  for (psi in psi_grid) {
    
    Beta_hat_con_fn <- Beta_hat_con_fn_template(X, weights, psi)
    
    Beta_hat <- get_Beta_hat(Beta_hat_obj_fn,
                             Beta_hat_obj_fn_gr,
                             Beta_hat_con_fn,
                             Beta_hat_con_fn_jac,
                             init_guess)
    
    log_L_b <- log_likelihood(Beta_hat, X, Y, t)
    
    log_L_b_vals <- c(log_L_b_vals, log_L_b)
    
    init_guess <- Beta_hat
  }
  
  list(psi = psi_grid, log_L_b = log_L_b_vals)
}

compute_IL_branch <- function(X,         
                              Y,              
                              t,              
                              weights,       
                              branch_params,
                              psi_bar,
                              psi_grid) {
  
  psi_mode <- branch_params$psi
  Beta_hat_mode <- branch_params$Beta_hat
  omega_hat <- branch_params$omega_hat
  
  psi_grid_left <- rev(psi_grid[psi_mode <= psi_bar])
  psi_grid_right <- psi_grid[psi_mode > psi_bar]
  
  # Evaluate branch to the left of the peak
  left_branch <- run_branch_side(
    X          = X,
    Y          = Y,
    t          = t,
    weights    = weights,
    omega_hat  = omega_hat,
    init_guess = Beta_hat_mode, 
    psi_grid   = psi_grid_left
  )
  
  # Evaluate branch to the right of the peak
  right_branch <- run_branch_side(
    X          = X,
    Y          = Y,
    t          = t,
    weights    = weights,
    omega_hat  = omega_hat,
    init_guess = Beta_hat_mode, 
    psi_grid   = psi_grid_right
  )
  
  # Combine and sort
  branch <- rbind(
    data.frame(psi = left_branch$psi, log_L_b = left_branch$log_L_b),
    data.frame(psi = right_branch$psi, log_L_b = right_branch$log_L_b)
  )
  
  branch <- branch[order(branch$psi), , drop = FALSE]
  rownames(branch) <- NULL
  return(branch)
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

get_integrated_LL <- function(config, data, X, weights) {
  
  invisible(list2env(config$optimization_specs$IL, env = environment()))
  
  Y <- data$Y
  t <- data$t
  model <- fit_model(config, data)
  search_interval <- get_search_interval(model, X, weights, num_std_errors)
  Beta_MLE <- get_Beta_MLE(model)
  psi_MLE <- get_psi(Beta_MLE, X, weights)
  
  omega_hat_con_fn <- omega_hat_con_fn_template(X, weights, psi_MLE)
  num_branches <- chunk_size * num_workers
  
  # -------------------------------
  # Outer loop: compute branch parameters
  # -------------------------------
  branch_params_list <- foreach(
    i = 1:num_branches,
    .combine = "list",
    .multicombine = TRUE,
    .errorhandling = "remove",
    .options.future = list(
      seed = TRUE,
      chunk.size = chunk_size,
      packages = c("nloptr", "dplyr")
      # globals = list(
      #   get_branch_params = get_branch_params,
      #   get_omega_hat = get_omega_hat,
      #   get_Beta_hat = get_Beta_hat,
      #   Beta_hat_obj_fn_template = Beta_hat_obj_fn_template,
      #   Beta_hat_obj_fn_gr_template = Beta_hat_obj_fn_gr_template,
      #   Beta_hat_con_fn_jac_template = Beta_hat_con_fn_jac_template,
      #   Beta_hat_con_fn_template = Beta_hat_con_fn_template,
      #   log_likelihood = log_likelihood,
      #   get_eta = get_eta,
      #   get_theta = get_theta,
      #   X = X,
      #   Y = Y,
      #   t = t,
      #   weights = weights,
      #   search_interval = search_interval,
      #   omega_hat_con_fn = omega_hat_con_fn,
      #   Beta_MLE_length = length(Beta_MLE)
      # )
    )
  ) %dofuture% {
    init_guess <- rnorm(length(Beta_MLE))
    get_omega_hat(omega_hat_con_fn, init_guess)
    # get_branch_params(omega_hat, X, Y, t, weights, search_interval)
  }

  return(branch_params_list)
  
  # psi_modes <- sapply(branch_params_list, `[[`, 1)
  # Beta_hat_modes <- lapply(branch_params_list, `[[`, 2)
  # omega_hats <- lapply(branch_params_list, `[[`, 3)
  
  # # -------------------------------
  # # Define psi grid
  # # -------------------------------
  # psi_bar <- mean(psi_modes)
  # psi_bar_SE <- sd(psi_modes)
  # psi_grid_endpoints <- psi_bar + c(-1, 1) * (max(abs(psi_modes - psi_bar)) + num_std_errors * psi_bar_SE)
  # psi_grid <- seq(psi_grid_endpoints[1], psi_grid_endpoints[2], increment)
  
  # # -------------------------------
  # # Inner loop: compute integrated likelihood branches
  # # -------------------------------
  # IL_branches <- foreach(
  #   branch_params = branch_params_list,
  #   .combine = "list",
  #   .multicombine = TRUE,
  #   .errorhandling = "remove",
  #   .options.future = list(
  #     seed = TRUE,
  #     chunk.size = chunk_size,
  #     packages = c("nloptr", "dplyr"),
  #     globals = list(
  #       compute_IL_branch = compute_IL_branch,
  #       run_branch_side = run_branch_side,
  #       X = X,
  #       Y = Y,
  #       t = t,
  #       weights = weights,
  #       psi_bar = psi_bar,
  #       psi_grid = psi_grid
  #     )
  #   )
  # ) %dofuture% {
  #   compute_IL_branch(
  #     X = X,
  #     Y = Y,
  #     t = t,
  #     weights = weights,
  #     branch_params = branch_params,
  #     psi_bar = psi_bar,
  #     psi_grid = psi_grid
  #   )
  # }
  
  # log_L_bar <- get_log_L_bar(IL_branches)
  
  # list(
  #   log_L_bar_df = log_L_bar$df,
  #   branches_matrix = log_L_bar$branches_matrix,
  #   IL_branches = IL_branches,
  #   omega_hats = omega_hats,
  #   branch_modes = psi_modes
  # )
}

# Profile Log-Likelihood --------------------------------------------------

compute_profile_branch <- function(direction,
                                   X,
                                   Y,
                                   t,
                                   Beta_MLE,
                                   weights,
                                   step_size,
                                   fine_step_size,
                                   fine_window,
                                   alpha) {
  
  log_L_p_max <- log_likelihood(Beta_MLE, X, Y, t)
  crit <- qchisq(1 - alpha, df = 1) / 2
  stopping_val <- log_L_p_max - crit
  log_L_p <- log_L_p_max
  
  psi_MLE <- get_psi(Beta_MLE, X, weights)
  psi_vals <- numeric(0)
  log_L_p_vals <- numeric(0)
  
  if (direction == "left") {
    step_anchor <- floor(psi_MLE / step_size) * step_size
    direction_sign <- -1
  } else {
    step_anchor <- ceiling(psi_MLE / step_size) * step_size
    direction_sign <- 1
  }
  
  Beta_hat_obj_fn <- Beta_hat_obj_fn_template(Beta_MLE, X, Y, t)
  Beta_hat_obj_fn_gr <- Beta_hat_obj_fn_gr_template(Beta_MLE, X, Y, t)
  
  Beta_hat_con_fn_jac <- Beta_hat_con_fn_jac_template(X, weights)
  
  psi <- psi_MLE + direction_sign * fine_step_size
  init_guess <- Beta_MLE
  
  while (log_L_p >= stopping_val) {

    Beta_hat_con_fn <- Beta_hat_con_fn_template(X, weights, psi)

    Beta_hat <- get_Beta_hat(Beta_hat_obj_fn,
                             Beta_hat_obj_fn_gr,
                             Beta_hat_con_fn,
                             Beta_hat_con_fn_jac,
                             init_guess)
    
    # Evaluate profile log-likelihood using Beta_hat
    log_L_p <- log_likelihood(Beta_hat, X, Y, t)
    
    log_L_p_vals <- c(log_L_p_vals, log_L_p)
    psi_vals <- c(psi_vals, psi)
    
    # next initial guess
    init_guess <- Beta_hat
    
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

get_profile_LL <- function(config, data, X, weights) {

  invisible(list2env(config$optimization_specs$PL, environment()))
  
  Y <- data$Y
  t <- data$t
  model <- fit_model(config, data)
  Beta_MLE <- get_Beta_MLE(model)
  
  alpha <- min(alpha_levels)
  step_size_local <- step_size
  fine_step_size_local <- fine_step_size
  fine_window_local <- fine_window
  Beta_MLE_length <- length(Beta_MLE)
  
  result <- foreach(
    dir = c("left", "right"),
    .combine = "list",
    .multicombine = TRUE,
    .errorhandling = "remove",
    .options.future = list(
      seed = TRUE,
      chunk.size = 1,
      packages = c("nloptr", "dplyr"),
      globals = list(
        compute_profile_branch = compute_profile_branch,
        Beta_MLE_length = Beta_MLE_length,
        X = X,
        Y = Y,
        t = t,
        weights = weights,
        alpha = alpha,
        step_size = step_size_local,
        fine_step_size = fine_step_size_local,
        fine_window = fine_window_local
      )
    )
  ) %dofuture% {
    compute_profile_branch(
      direction      = dir,
      X              = X,
      Y              = Y,
      t              = t,
      Beta_MLE       = Beta_MLE,
      weights        = weights,
      step_size      = step_size,
      fine_step_size = fine_step_size,
      fine_window    = fine_window,
      alpha          = alpha
    )
  }
  
  # Combine results into a single data frame
  profile_LL <- do.call(
    rbind,
    lapply(result, \(entry) data.frame(psi = entry$psi, Profile = entry$Profile))
  )
  
  profile_LL
}

get_report_objects <- function(iter_dir) {
  
  results_dir <- here(iter_dir, "results")
  config_path <- here(iter_dir, "config_snapshot.yml")
  config <- read_yaml(config_path)
  
  integrated_LL <- readRDS(here(results_dir, "integrated_LL.rds"))
  profile_LL <- readRDS(here(results_dir, "profile_LL.rds"))
  
  log_L_bar_df <- integrated_LL$log_L_bar_df
  
  LL_df <- log_L_bar_df |>
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
