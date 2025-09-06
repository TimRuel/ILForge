# applications/poisson/group_rates_weighted_sum/fixed_effects_regression/scripts/helpers/experiment_utils.R

# Miscellaneous -----------------------------------------------------------

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

get_Beta_MLE <- function(model) {
  
  coefs <- coef(model)
  Beta_MLE <- as.matrix(coefs, ncol = 1)  # stack as single-column matrix
  
  # Clean the names
  names_clean <- names(coefs)
  
  # Intercepts: "groupA" → "A"
  names_clean <- sub("^group", "", names_clean)
  
  # Varying slopes: "groupA:X2" → "A:X2"
  names_clean <- sub("^group([A-Z]+):", "\\1:", names_clean)
  
  rownames(Beta_MLE) <- names_clean
  
  Beta_MLE
}

get_theta <- function(Beta, G) exp(Beta[1:G,])

get_psi <- function(Beta, weights) {
  
  theta <- get_theta(Beta, length(weights))
  
  sum(theta * weights)
}

# General delta-method version of SE(psi_hat)
get_psi_MLE_SE <- function(model, weights, X_ref = NULL, ref_covariates = NULL) {
  # Extract Beta_MLE as numeric vector
  Beta_MLE <- as.numeric(get_Beta_MLE(model))
  Beta_cov <- vcov(model)
  
  # Group labels
  group_labels <- levels(model$data$group)
  G <- length(group_labels)
  
  # If X_ref is not provided, build one for reference covariates
  if (is.null(X_ref)) {
    total_cols <- length(Beta_MLE)
    X_ref <- matrix(0, nrow = G, ncol = total_cols)
    
    # Intercept columns
    X_ref[cbind(1:G, 1:G)] <- 1
    
    # Shared covariates evaluated at ref_covariates or 0
    if (!is.null(ref_covariates)) {
      shared_idx <- (G + 1):(total_cols - (G * length(ref_covariates$varying))) 
      if (length(shared_idx) > 0) {
        X_ref[, shared_idx] <- matrix(rep(ref_covariates$shared, each = G), nrow = G)
      }
    }
    
    # Varying slopes evaluated at ref_covariates
    if (!is.null(ref_covariates) && length(ref_covariates$varying) > 0) {
      start <- total_cols - G * length(ref_covariates$varying) + 1
      for (i in seq_along(ref_covariates$varying)) {
        X_ref[, start:(start+G-1)] <- diag(ref_covariates$varying[i], G)
        start <- start + G
      }
    }
  }
  
  # Compute gradient
  eta <- X_ref %*% Beta_MLE
  grad <- t(X_ref) %*% (weights * exp(eta))
  
  # Delta-method SE
  sqrt(as.numeric(t(grad) %*% Beta_cov %*% grad))
}

get_eta <- function(Beta, X) X %*% Beta

get_mu <- function(eta, t) t * exp(eta)
  
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
                         weights, 
                         num_std_errors, 
                         step_size,
                         fine_step_size,
                         fine_window) {
  
  Beta_MLE <- get_Beta_MLE(model)
  
  psi_MLE <- get_psi(Beta_MLE, weights)
  
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

safe_auglag <- possibly(nloptr::auglag)

Beta_hat_con_fn_template <- function(weights, G, psi) {
  
  function(Beta) get_psi(Beta, weights) - psi
}

Beta_hat_obj_fn_template <- function(omega_hat, X, Y, t) {
  
  function(Beta) {
    
    E_Y <- t * exp(get_eta(omega_hat, X))
    
    -log_likelihood(Beta, X, E_Y, t)
  }
}

get_Beta_hat <- function(Beta_hat_obj_fn, Beta_hat_con_fn, init_guess) {
  
  safe_auglag(
    x0 = init_guess,
    fn = Beta_hat_obj_fn,
    heq = Beta_hat_con_fn,
    localsolver = "SLSQP",
    control = list(xtol_rel = 1e-8, maxeval = 1000),
    deprecatedBehavior = FALSE)$par
}
  
# omega_hat ---------------------------------------------------------------

get_omega_hat <- function(psi_MLE, weights, n_covariates,
                          v_sampler = function(n) rgamma(n, shape = 2, rate = 2),
                          gamma_sampler = function(n) rnorm(n, 0, 1),
                          B = 20) {
  G <- length(weights)
  # 1) draw base positives
  v <- v_sampler(G)
  denom <- sum(weights * v)
  if (denom <= 0) stop("Denominator non-positive.")
  cscale <- psi_MLE / denom
  u0 <- cscale * v
  beta0 <- log(u0)
  
  # 2) clamp / project u so that corresponding beta in [-B, B]
  lower_u <- exp(-B)
  upper_u <- exp(B)
  clamped <- (u0 < lower_u) | (u0 > upper_u)
  
  if (any(clamped)) {
    # set clamped u to bounds; solve for multiplier for the remaining indices
    u <- u0
    u[ u < lower_u ] <- lower_u
    u[ u > upper_u ] <- upper_u
    
    # Now find multiplier c2 on the unclamped indices to satisfy sum(weights * u) = psi_MLE
    unclamped_idx <- which(!clamped)
    if (length(unclamped_idx) == 0) {
      # all clamped: check feasibility
      total_clamped_sum <- sum(weights * u)
      if (abs(total_clamped_sum - psi_MLE) > 1e-12) {
        stop("Constraint infeasible after clamping.")
      }
    } else {
      # solve for c2: sum(weights_clamped * u_clamped) + c2 * sum(weights_unclamped * v_unclamped) = psi_MLE
      sum_clamped <- sum(weights[clamped] * u[clamped])
      sum_unclamped_wv <- sum(weights[unclamped_idx] * v[unclamped_idx])
      c2 <- (psi_MLE - sum_clamped) / sum_unclamped_wv
      if (c2 <= 0) {
        # If c2 <= 0, the clamp bounds are too tight; fall back to mild epsilon floor (safe)
        eps <- 1e-8
        u <- pmax(u, eps)
        c2 <- (psi_MLE - sum(weights[clamped] * u[clamped])) / sum_unclamped_wv
        if (c2 <= 0) stop("Unable to find feasible c2 — adjust B or v_sampler.")
      }
      u[unclamped_idx] <- c2 * v[unclamped_idx]
    }
  } else {
    u <- u0
  }
  
  beta <- log(u)
  gammas <- gamma_sampler(n_covariates)
  c(beta, gammas)
}

# Integrated Log-Likelihood ---------------------------------------------------

run_branch_side <- function(direction, 
                            omega_hat,
                            data,
                            weights,
                            psi_MLE,
                            psi_grid) {
  
  # Split the grid depending on direction
  psi_working_grid <- if (direction == "left") {
    rev(psi_grid[psi_grid <= psi_MLE])
  } else {
    psi_grid[psi_grid > psi_MLE]
  }
  
  covariates <- data |> select(starts_with("X"))
  n_per_group <- table(data$group)
  G <- length(n_per_group)
  
  X <- get_X(covariates, n_per_group)
  Y <- data$Y
  t <- data$t
  
  log_L_tilde_vals <- c()
  
  # --- Compute feasible Beta for first psi ---
  psi_start <- psi_working_grid[1]
  Beta_hat_con_fn <- Beta_hat_con_fn_template(weights, G, psi_start)
  Beta_hat_obj_fn <- Beta_hat_obj_fn_template(omega_hat, X, Y, t)
  
  # Initial guess: attempt feasible solution
  init_guess <- tryCatch({
    get_Beta_hat(Beta_hat_obj_fn, Beta_hat_con_fn, omega_hat)
  }, error = function(e) {
    message("[WARN] Initial guess not feasible; projecting to feasible point")
    # Simple projection: scale intercepts to satisfy psi constraint
    delta <- psi_start - sum(weights * exp(omega_hat[1:G]))
    omega_hat[1:G] <- omega_hat[1:G] + log(1 + delta / sum(exp(omega_hat[1:G])))
    omega_hat
  })
  
  # --- Loop through psi grid ---
  for (psi in psi_working_grid) {
    
    Beta_hat_con_fn <- Beta_hat_con_fn_template(weights, G, psi)
    Beta_hat_obj_fn <- Beta_hat_obj_fn_template(init_guess, X, Y, t)
    
    # Attempt to compute Beta_hat; project if infeasible
    Beta_hat <- tryCatch({
      get_Beta_hat(Beta_hat_obj_fn, Beta_hat_con_fn, init_guess)
    }, error = function(e) {
      message("[WARN] Step infeasible at psi = ", psi, "; projecting solution")
      # Simple projection: scale intercepts
      delta <- psi - sum(weights * exp(init_guess[1:G]))
      init_guess[1:G] <- init_guess[1:G] + log(1 + delta / sum(exp(init_guess[1:G])))
      init_guess
    })
    
    log_L_tilde <- log_likelihood(Beta_hat, X, Y, t)
    
    log_L_tilde_vals <- c(log_L_tilde_vals, log_L_tilde)
    
    # Update initial guess for next step
    init_guess <- Beta_hat
  }
  
  # Reverse vectors for left branch to keep increasing order of psi
  if (direction == "left") {
    psi_working_grid <- rev(psi_working_grid)
    log_L_tilde_vals <- rev(log_L_tilde_vals)
  }
  
  list(psi = psi_working_grid, Integrated = log_L_tilde_vals)
}

# run_branch_side <- function(direction, 
#                             omega_hat,
#                             data,
#                             weights,
#                             psi_MLE,
#                             psi_grid) {
#   
#   psi_working_grid <- if (direction == "left") {
#     
#     rev(psi_grid[psi_grid <= psi_MLE])
#   } else {
#     
#     psi_grid[psi_grid > psi_MLE]
#   }
#   
#   covariates <- data |> 
#     select(starts_with("X"))
#   
#   n_per_group <- table(data$group)
#   
#   X <- get_X(covariates, n_per_group)
#   
#   Y <- data$Y
#   
#   t <- data$t
#   
#   log_L_tilde_vals <- c()
#   
#   init_guess <- omega_hat
#   
#   for (psi in psi_working_grid) {
#     
#     Beta_hat_con_fn <- Beta_hat_con_fn_template(weights, G, psi)
#     
#     Beta_hat_obj_fn <- Beta_hat_obj_fn_template(omega_hat, X, Y, t)
#     
#     Beta_hat <- get_Beta_hat(Beta_hat_obj_fn, Beta_hat_con_fn, init_guess)
#     
#     log_L_tilde <- log_likelihood(Beta_hat, X, Y, t)
#     
#     log_L_tilde_vals <- c(log_L_tilde_vals, log_L_tilde)
#     
#     init_guess <- Beta_hat
#   }
#   
#   if (direction == "left") {
#     
#     psi_working_grid <- rev(psi_working_grid)
#     log_L_tilde_vals <- rev(log_L_tilde_vals)
#   }
#   
#   list(psi = psi_working_grid, Integrated = log_L_tilde_vals)
# }

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

flag_curve_outliers <- function(df, x_col = "x", y_col = "y", remove = FALSE, span = 0.3, threshold = 3) {
  stopifnot(x_col %in% names(df), y_col %in% names(df))
  x <- df[[x_col]]
  y <- df[[y_col]]
  
  # fit smooth curve
  fit <- loess(y ~ x, span = span)
  y_hat <- predict(fit, x)
  
  # standardized residuals
  res <- y - y_hat
  zres <- (res - mean(res, na.rm = TRUE)) / sd(res, na.rm = TRUE)
  
  outlier <- abs(zres) > threshold
  df$outlier <- outlier
  if (remove) df <- df[!outlier, , drop = FALSE] |> 
    dplyr::select(-outlier)
  df
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
  
  model <- fit_model(data)
  
  psi_grid <- get_psi_grid(model, 
                           weights, 
                           num_std_errors, 
                           step_size,
                           fine_step_size,
                           fine_window)
  
  Beta_MLE <- get_Beta_MLE(model)
  
  G <- nlevels(data$group)
  
  n_covariates <- data |> 
    select(starts_with("X")) |> 
    names() |> 
    length()
  
  psi_MLE <- get_psi(Beta_MLE, weights)
  
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
    
    omega_hat <- get_omega_hat(psi_MLE, weights, n_covariates)
    
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
       branches_matrix = log_L_bar$branches_matrix,
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
  
  covariates <- data |> select(starts_with("X"))
  n_per_group <- table(data$group)
  G <- length(n_per_group)
  X <- get_X(covariates, n_per_group)
  Y <- data$Y
  t <- data$t
  
  model <- fit_model(data)
  Beta_MLE <- get_Beta_MLE(model)
  
  log_L_p_max <- log_likelihood(Beta_MLE, X, Y, t)
  crit <- qchisq(1 - alpha, df = 1) / 2
  stopping_val <- log_L_p_max - crit
  log_L_p <- log_L_p_max
  
  psi_MLE <- get_psi(Beta_MLE, weights)
  psi_vals <- numeric(0)
  log_L_p_vals <- numeric(0)
  
  if (direction == "left") {
    step_anchor <- floor(psi_MLE / step_size) * step_size
    direction_sign <- -1
  } else {
    step_anchor <- ceiling(psi_MLE / step_size) * step_size
    direction_sign <- 1
  }
  
  psi <- psi_MLE + direction_sign * fine_step_size
  init_guess <- Beta_MLE
  
  while (log_L_p >= stopping_val) {
    # build objective and constraint *functions* for current psi
    Beta_hat_con_fn <- Beta_hat_con_fn_template(weights, G, psi)
    Beta_hat_obj_fn <- Beta_hat_obj_fn_template(Beta_MLE, X, Y, t)
    
    # call constrained estimator, but protect against failure
    Beta_hat_try <- tryCatch(
      {
        get_Beta_hat(Beta_hat_obj_fn, Beta_hat_con_fn, init_guess)
      },
      error = function(e) {
        message("[WARN] get_Beta_hat failed at psi = ", psi, " : ", conditionMessage(e))
        NULL
      },
      warning = function(w) {
        message("[WARN] get_Beta_hat warning at psi = ", psi, " : ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    if (is.null(Beta_hat_try)) {
      # If constraint infeasible or optimizer failed, stop extending the profile
      message("[INFO] Stopping profile branch early at psi = ", psi)
      break
    }
    
    Beta_hat <- Beta_hat_try
    
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
  
  # log_L_bar_df <- flag_curve_outliers(integrated_LL$log_L_bar_df, x_col = "psi", y_col = "Integrated", remove = TRUE)
  
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
