# applications/multinomial_logistic_regression/entropy/scripts/helpers/experiment_utils.R

# Miscellaneous -----------------------------------------------------------

fit_multinomial_logistic_model <- function(data, formula) {
  
  nnet::multinom(formula, data, trace = FALSE)
}

log_likelihood <- function(Beta, X_design, Y_design) {
  
  p <- ncol(X_design)
  
  Jm1 <- ncol(Y_design)
  
  Beta <- matrix(Beta,
                 nrow = p,
                 ncol = Jm1,
                 byrow = FALSE)
  
  Y_hat <- X_design %*% Beta
  
  sum(rowSums(Y_design * Y_hat) - matrixStats::rowLogSumExps(cbind(0, Y_hat)))
}

likelihood <- function(Beta, X_design, Y_design) exp(log_likelihood(Beta, X_design, Y_design))

get_Beta_MLE <- function(model) {
  
  Jm1 <- length(model$lev) - 1
  
  model |>
    coef() |>
    unname() |>
    matrix(ncol = Jm1,
           byrow = TRUE)
}

get_psi_hat_from_model <- function(model, X1_levels) {
  
  X1_level_names <- names(X1_levels)
  
  m <- map_int(X1_levels, \(x) x$m)
  
  h <- get_X1_level_of_interest(X1_levels)
  
  model |>
    predict(type = "probs") |>
    as.data.frame() |>
    mutate(X1_level = rep(X1_level_names, times = m)) |>
    aggregate(. ~ X1_level, data = _, FUN = mean) |>
    (\(df) setNames(data.frame(t(df[,-1])), df[,1]))() |>
    apply(2, entropy) |>
    (\(vec) vec[h])() |> 
    unname()
}

get_psi_hat_from_MLE <- function(Beta_MLE, X_h_design) {
  
  (X_h_design %*% cbind(0, Beta_MLE)) |> 
    apply(1, softmax) |> 
    t() |> 
    colMeans() |> 
    entropy()
}

safe_auglag <- purrr::possibly(nloptr::auglag)

# Psi Grid ----------------------------------------------------------------

get_psi_endpoints <- function(psi_hat, Beta_MLE, X_h_design, num_std_errors, J, n_h) {
  
  theta_MLE <- X_h_design %*% cbind(0, Beta_MLE) |>
    apply(1, softmax) |>
    t() |>
    colMeans()
  
  sigma <- theta_MLE*diag(J) - matrix(theta_MLE) %*% theta_MLE
  
  psi_hat_SE <- sqrt(sum(matrix(1 + log(theta_MLE)) %*% (1 + log(theta_MLE)) * sigma, na.rm = TRUE) / n_h)
  
  psi_endpoints <- psi_hat + c(-1, 1) * num_std_errors * psi_hat_SE
}

get_psi_grid <- function(psi_endpoints, step_size, J) {
  
  lower <- 0
  
  upper <- log(J)
  
  left <- max(psi_endpoints[1], lower)
  right <- min(psi_endpoints[2], upper)
  
  psi_grid <- seq(ceiling(left / step_size) * step_size,
                  floor(right / step_size) * step_size,
                  by = step_size)
  
  if (psi_endpoints[1] < lower) psi_grid <- c(lower, psi_grid[psi_grid > lower])
  if (psi_endpoints[2] > upper) psi_grid <- c(psi_grid[psi_grid < upper], upper)
  
  return(psi_grid)
}

safe_max <- function(x, fallback) if (length(x) == 0) fallback else max(x)
safe_min <- function(x, fallback) if (length(x) == 0) fallback else min(x)

get_fine_psi_grid <- function(psi_grid,
                              psi_mode,
                              fine_step_size,
                              fine_window) {
  
  # Ensure psi_grid is sorted and unique
  psi_grid <- sort(unique(psi_grid))
  
  # Determine left and right limits using nearest psi_grid values
  lower_bound <- safe_max(psi_grid[psi_grid <= (psi_mode - fine_window)], min(psi_grid))
  upper_bound <- safe_min(psi_grid[psi_grid >= (psi_mode + fine_window)], max(psi_grid))
  
  # Build fine psi values from lower_bound to upper_bound
  fine_grid <- seq(from = lower_bound,
                   to = upper_bound,
                   by = fine_step_size) |> 
    round(8)
  
  # Combine with original grid and deduplicate + sort
  full_grid <- c(psi_grid, fine_grid) |> 
    round(8) |> 
    unique() |> 
    sort()
  
  return(full_grid)
}

# Beta_hat ----------------------------------------------------------------

get_Beta_hat_template <- function(obj_fn, con_fn, init_guess) {
  
  safe_auglag(
    x0 = init_guess,
    fn = obj_fn,
    heq = con_fn,
    localsolver = "SLSQP",
    deprecatedBehavior = FALSE)$par
}

# Branch Parameters ---------------------------------------------------------------

get_threshold <- function(Beta_MLE, X_design, Y_design, threshold_offset) {
  
  ceiling(abs(log_likelihood(Beta_MLE, X_design, Y_design))) + threshold_offset
}

get_omega_hat_branch_mode <- function(omega_hat,
                                      X_design,
                                      Y_design,
                                      X_h_design, 
                                      Jm1,
                                      p,
                                      n) {
  
  Beta_hat_obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, omega_hat, Jm1, p, n)
  
  omega_hat_branch <- function(psi) {
    
    Beta_hat_con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
    
    Beta_hat <- get_Beta_hat_template(Beta_hat_obj_fn, Beta_hat_con_fn, c(omega_hat))
    
    return(log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n))
  }
  
  optimize(omega_hat_branch,
           interval = c(0, log(Jm1 + 1)),
           maximum = TRUE,
           tol = 0.1)$maximum
}

get_branch_params <- function(X_design,
                              Y_design,
                              X_h_design,
                              threshold,
                              psi_hat,
                              psi_grid,
                              Jm1, 
                              p, 
                              n,
                              init_guess_sd) {
  
  omega_hat_obj_fn <- function(Beta) 0
  
  omega_hat_eq_con_fn <- function(Beta) omega_hat_eq_con_fn_rcpp(Beta, X_h_design, Jm1, p, psi_hat)
  
  omega_hat_ineq_con_fn <- function(Beta) omega_hat_ineq_con_fn_rcpp(Beta, X_design, Y_design, Jm1, p, n, threshold)
  
  num_discarded <- 0
  
  while (TRUE) {
    
    init_guess <- rnorm(p * Jm1, sd = init_guess_sd)
    
    omega_hat_candidate <- safe_auglag(
      x0 = init_guess,
      fn = omega_hat_obj_fn,
      heq = omega_hat_eq_con_fn,
      hin = omega_hat_ineq_con_fn,
      localsolver = "SLSQP",
      deprecatedBehavior = FALSE,
      control = list(on.error = "ignore")
      )$par
    
    if (!is.null(omega_hat_candidate) && abs(omega_hat_eq_con_fn(omega_hat_candidate)) <= 0.1 && omega_hat_ineq_con_fn(omega_hat_candidate) <= 0) {
      
      omega_hat_candidate <- matrix(omega_hat_candidate, nrow = p, ncol = Jm1, byrow = FALSE)
      
      omega_hat_candidate_branch_mode <- get_omega_hat_branch_mode(
        omega_hat = omega_hat_candidate,
        X_design = X_design,
        Y_design = Y_design,
        X_h_design = X_h_design, 
        Jm1 = Jm1,
        p = p,
        n = n
        )
      
      if (between(omega_hat_candidate_branch_mode, head(psi_grid, 1), tail(psi_grid, 1))) {
        
        Beta_hat_obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, omega_hat_candidate, Jm1, p, n)
        
        Beta_hat_con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, omega_hat_candidate_branch_mode, Jm1, p)
        
        Beta_hat <- get_Beta_hat_template(Beta_hat_obj_fn, Beta_hat_con_fn, c(omega_hat_candidate)) |> 
          matrix(nrow = p, ncol = Jm1, byrow = FALSE)
        
        return(list(omega_hat = omega_hat_candidate,
                    branch_mode = list(Beta_MLE = Beta_hat,
                                       psi = omega_hat_candidate_branch_mode),
                    num_discarded = num_discarded))
      }
    }
    num_discarded <- num_discarded + 1
  }
}

# Integrated Log-Likelihood ---------------------------------------------------

run_branch_side <- function(direction, 
                            branch_mode, 
                            psi_grid, 
                            step_size,
                            fine_step_size,
                            fine_window,
                            log_likelihood_fn,
                            get_Beta_hat,
                            con_fn_template) {
  
  # Get finer working grid from psi_grid and psi_hat ± band
  fine_grid <- get_fine_psi_grid(
    psi_grid       = psi_grid,
    psi_mode       = branch_mode$psi,
    fine_step_size = fine_step_size,
    fine_window    = fine_window
  )
  
  # Determine which psi values to evaluate depending on direction
  psi_working_grid <- if (direction == "left") {
    rev(fine_grid[fine_grid < branch_mode$psi])
  } else {
    fine_grid[fine_grid > branch_mode$psi]
  }
  
  psi_vals <- c()
  log_L_tilde_vals <- c()
  init_guess <- branch_mode$Beta_MLE
  
  for (psi in psi_working_grid) {
    Beta_hat_con_fn <- function(Beta) con_fn_template(Beta, psi)
    Beta_hat <- get_Beta_hat(Beta_hat_con_fn, init_guess)
    log_L_tilde <- log_likelihood_fn(Beta_hat)
    init_guess <- Beta_hat
    
    # Only keep results for psi values in original psi_grid
    if (any(abs(psi_grid - psi) < 1e-8)) {
      psi_vals <- c(psi_vals, psi)
      log_L_tilde_vals <- c(log_L_tilde_vals, log_L_tilde)
    }
  }
  
  # Reverse the output if we reversed the grid
  if (direction == "left") {
    psi_vals <- rev(psi_vals)
    log_L_tilde_vals <- rev(log_L_tilde_vals)
  }
  
  list(psi = psi_vals, Integrated = log_L_tilde_vals)
}

compute_IL_branch <- function(branch_mode,
                              psi_grid,
                              step_size,
                              fine_step_size,
                              fine_window,
                              log_likelihood_fn,
                              get_Beta_hat,
                              con_fn_template) {
  
  # Evaluate branch to the left of the peak
  left_result <- run_branch_side(
    direction          = "left",
    branch_mode        = branch_mode,
    psi_grid           = psi_grid,
    step_size          = step_size,
    fine_step_size     = fine_step_size,
    fine_window        = fine_window,
    log_likelihood_fn  = log_likelihood_fn,
    get_Beta_hat       = get_Beta_hat,
    con_fn_template    = con_fn_template
  )
  
  # Evaluate branch to the right of the peak
  right_result <- run_branch_side(
    direction          = "right",
    branch_mode        = branch_mode,
    psi_grid           = psi_grid,
    step_size          = step_size,
    fine_step_size     = fine_step_size,
    fine_window        = fine_window,
    log_likelihood_fn  = log_likelihood_fn,
    get_Beta_hat       = get_Beta_hat,
    con_fn_template    = con_fn_template
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

get_integrated_LL <- function(config, X_design, model_df) {
  
  invisible(list2env(config$model_specs, environment()))
  invisible(list2env(config$optimization_specs$IL, environment()))
  
  X1_levels <- config$X1_levels
  formula <- as.formula(formula)
  ml_model <- fit_multinomial_logistic_model(model_df, formula)
  Y_design <- get_Y_design(model_df)
  Beta_MLE <- get_Beta_MLE(ml_model)
  threshold <- get_threshold(Beta_MLE, X_design, Y_design, threshold_offset)
  Jm1 <- J - 1
  h <- get_X1_level_of_interest(X1_levels)
  X_h_design <- get_X_h_design(X_design, X1_levels)
  psi_hat <- get_psi_hat_from_model(ml_model, X1_levels)
  n_h <- nrow(X_h_design)
  psi_endpoints <- get_psi_endpoints(psi_hat, Beta_MLE, X_h_design, num_std_errors, J, n_h)
  psi_grid <- get_psi_grid(psi_endpoints, step_size, J)
  num_branches <- num_workers * chunk_size
  
  log_likelihood_fn <- function(Beta) log_likelihood_rcpp(Beta, X_design, Y_design, Jm1, p, n)
  
  result <- foreach(
    
    i = 1:num_branches,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = num_branches,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size,
                           packages = c("PolytomousUtils", "nloptr"))
    
  ) %dofuture% {
    
    branch_params <- get_branch_params(
      X_design = X_design,
      Y_design = Y_design,
      X_h_design = X_h_design, 
      threshold = threshold,
      psi_hat = psi_hat,
      psi_grid = psi_grid,
      Jm1 = Jm1, 
      p = p, 
      n = n,
      init_guess_sd = init_guess_sd
      )

    Beta_hat_obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, branch_params$omega_hat, Jm1, p, n)
    get_Beta_hat <- function(con_fn, init_guess) get_Beta_hat_template(Beta_hat_obj_fn, con_fn, init_guess)
    con_fn_template <- function(Beta, psi) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
    
    IL_branch <- compute_IL_branch(
      branch_mode = branch_params$branch_mode,
      psi_grid = psi_grid,
      fine_step_size = fine_step_size,
      fine_window = fine_window,
      log_likelihood_fn = log_likelihood_fn,
      get_Beta_hat = get_Beta_hat,
      con_fn_template = con_fn_template
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
                                   step_size,
                                   fine_step_size,
                                   fine_window,
                                   psi_hat,
                                   Beta_MLE,
                                   PLL_max,
                                   stopping_val,
                                   log_likelihood_fn,
                                   get_Beta_hat,
                                   con_fn_template) {
  
  psi_vals <- list()
  log_L_p_vals <- list()
  init_guess <- Beta_MLE
  log_L_p <- PLL_max
  
  # Compute step anchor point and direction sign (determines the step direction based on branch side)
  if (direction == "left") {
    step_anchor <- floor(psi_hat / step_size) * step_size
    direction_sign <- -1
  } else {
    step_anchor <- ceiling(psi_hat / step_size) * step_size
    direction_sign <- 1
  }
  
  psi <- psi_hat + direction_sign * fine_step_size
  
  while (log_L_p > stopping_val) {
    
    # Compute constrained MLE
    Beta_hat_con_fn <- function(Beta) con_fn_template(Beta, psi)
    Beta_hat <- get_Beta_hat(Beta_hat_con_fn, init_guess)
    log_L_p <- log_likelihood_fn(Beta_hat)
    init_guess <- Beta_hat
    
    # Record value
    psi_vals[[length(psi_vals) + 1]] <- psi
    log_L_p_vals[[length(log_L_p_vals) + 1]] <- log_L_p
    
    # Adjust step size dynamically
    dist_from_peak <- abs(psi - psi_hat)
    use_fine_step <- dist_from_peak < fine_window || 
      (direction == "left" && psi > step_anchor) ||
      (direction == "right" && psi < step_anchor)
    current_step <- if (use_fine_step) fine_step_size else step_size
    psi <- psi + direction_sign * current_step
  }
  
  psi_vals <- unlist(psi_vals)
  log_L_p_vals <- unlist(log_L_p_vals)
  
  # Append peak for left branch
  if (direction == "left") {
    psi_vals <- c(rev(psi_vals), psi_hat)
    log_L_p_vals <- c(rev(log_L_p_vals), PLL_max)
  }
  
  list(psi = psi_vals, Profile = log_L_p_vals)
}

get_profile_LL <- function(config, X_design, model_df) {
  
  invisible(list2env(config$model_specs, environment()))
  invisible(list2env(config$optimization_specs$PL, environment()))
  
  X1_levels <- config$X1_levels
  formula <- as.formula(formula)
  ml_model <- fit_multinomial_logistic_model(model_df, formula)
  Beta_MLE <- get_Beta_MLE(ml_model)
  Jm1 <- J - 1
  h <- get_X1_level_of_interest(X1_levels)
  X_h_design <- get_X_h_design(X_design, X1_levels)
  psi_hat <- get_psi_hat_from_model(ml_model, X1_levels)
  Y_design <- get_Y_design(model_df)

  log_likelihood_fn <- function(Beta) log_likelihood_rcpp(Beta, X_design, Y_design, Jm1, p, n)
  Beta_hat_obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, Beta_MLE, Jm1, p, n)
  get_Beta_hat <- function(con_fn, init_guess) get_Beta_hat_template(Beta_hat_obj_fn, con_fn, init_guess)
  con_fn_template <- function(Beta, psi) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
  
  PLL_max <- log_likelihood_fn(Beta_MLE)
  alpha <- min(alpha_levels)
  crit <- qchisq(1 - alpha, df = 1) / 2
  stopping_val <- PLL_max - crit
  
  result <- foreach(
    dir = c("left", "right"),
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = 2,
    .errorhandling = "remove",
    .options.future = list(
      seed = TRUE,
      chunk.size = 1,
      packages = c("PolytomousUtils", "nloptr")
    )
  ) %dofuture% {
    
    compute_profile_branch(
      direction = dir,
      step_size = step_size,
      fine_step_size = fine_step_size,
      fine_window = fine_window,
      psi_hat = psi_hat,
      Beta_MLE = Beta_MLE,
      PLL_max = PLL_max,
      stopping_val = stopping_val,
      log_likelihood_fn = log_likelihood_fn,
      get_Beta_hat = get_Beta_hat,
      con_fn_template = con_fn_template
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
  
  J <- config$model_specs$J
  
  conf_ints <- get_confidence_intervals(
    pseudolikelihoods = pseudolikelihoods,
    LL_df_long = LL_df_long,
    MLE_data = MLE_data,
    alpha_levels = alpha_levels,
    J = J
  )
  
  return(list(MLE_data = MLE_data,
              conf_ints = conf_ints))
}




