# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/integrated_likelihood_utils.R

# ======================================================================
# Integrated Likelihood Utilities (Negative Binomial, Weighted-Sum ψ)
# ======================================================================

# ======================================================================
# One-point evaluation of the branch log-likelihood at a given ψ
# ======================================================================

.eval_branch_ll_at_psi <- function(psi, init_guess, ...) {
  closures <- .make_theta_phi_hat_closures(
    psi           = psi,
    theta_phi_0   = theta_phi_0,
    t             = t,
    n_per_process = n_per_process,
    weights       = weights,
    p_cutoff      = p_cutoff,
    max_y_cap     = max_y_cap
  )
  
  theta_phi_hat <- .get_theta_phi_hat(
    closures    = closures,
    init_guess  = init_guess,
    localsolver = localsolver,
    control     = control
  )
  
  J <- length(theta_phi_0$theta)
  
  branch_val <- log_likelihood(
    theta         = theta_phi_hat[seq_len(J)],
    phi           = theta_phi_hat[(J + 1):(2 * J)],
    Y             = Y,
    t             = t,
    n_per_process = n_per_process
  )
  
  list(
    branch_val    = branch_val,
    theta_phi_hat = theta_phi_hat
  )
}


# ======================================================================
# Branch maximization (locating the ψ-mode for a given branch)
# ======================================================================

.maximize_branch <- function(search_interval,
                             init_guess,
                             eval_psi_fun) {
  obj <- function(psi) {
    -eval_psi_fun(psi, init_guess)$branch_val
  }
  
  psi_hat <- stats::optimize(
    f        = obj,
    interval = search_interval
  )$minimum
  
  at_mode <- eval_psi_fun(psi_hat, init_guess)
  
  list(
    argmax = list(
      psi       = psi_hat,
      theta_phi = at_mode$theta_phi_hat
    ),
    max = at_mode$branch_val
  )
}

# ======================================================================
# Merge branches via log-sum-exp
# ======================================================================

.get_log_L_bar <- function(IL_branches) {
  merged <- Reduce(
    f = function(x, y) dplyr::inner_join(x, y, by = "k"),
    x = IL_branches
  )
  
  psi <- merged |>
    dplyr::select(dplyr::starts_with("psi")) |>
    dplyr::pull(1)
  
  branch_mat <- merged |>
    dplyr::arrange(.data$k) |>
    dplyr::select(dplyr::starts_with("value")) |>
    unname() |>
    as.matrix()
  
  R_eff <- ncol(branch_mat)
  lse   <- matrixStats::rowLogSumExps(branch_mat, na.rm = TRUE)
  
  log_L_bar <- lse - log(R_eff)
  
  rownames(branch_mat) <- as.character(psi)
  
  list(
    df              = tibble::tibble(psi = psi, value = as.numeric(log_L_bar)),
    branches_matrix = branch_mat
  )
}

# ======================================================================
# Required per-branch α for guaranteeing global IL cutoff coverage
# ======================================================================

.compute_required_branch_alpha <- function(R, alpha) {
  c_global      <- 0.5 * stats::qchisq(1 - alpha, df = 1)
  c_branch      <- c_global + log(R)
  alpha_branch  <- 1 - stats::pchisq(2 * c_branch, df = 1)
  alpha_branch
}

# ======================================================================
# Closure: construct eval_psi_fun for branch computations
# ======================================================================

.make_eval_psi_fun <- function(data_args,
                               likelihood_args,
                               solver_args) {
  force(data_args)
  force(likelihood_args)
  force(solver_args)
  
  function(psi, init_guess) {
    .eval_branch_ll_at_psi(
      psi           = psi,
      init_guess    = init_guess,
      Y             = data_args$Y,
      t             = data_args$t,
      n_per_process = data_args$n_per_process,
      weights       = data_args$weights,
      theta_phi_0   = likelihood_args$theta_phi_0,
      p_cutoff      = likelihood_args$p_cutoff,
      max_y_cap     = likelihood_args$max_y_cap,
      localsolver   = solver_args$localsolver,
      control       = solver_args$control
    )
  }
}

# ======================================================================
# Top-level integrated log-likelihood (public interface)
# ======================================================================

get_integrated_LL <- function(config, data, weights) {
  il_cfg <- config$optimization_specs$IL
  il_cfg <- sanitize_config(il_cfg)
  invisible(list2env(il_cfg, envir = environment()))
  
  global_alpha <- min(alpha_levels)
  num_branches <- max(1L, chunk_size * num_workers)
  
  Y             <- data$Y
  t             <- data$t
  n_per_process <- table(data$process)
  J             <- length(n_per_process)
  
  model <- fit_model(data)
  
  theta_MLE      <- get_theta_MLE(model)
  phi_MLE        <- get_phi_MLE(model)
  theta_phi_MLE  <- c(theta_MLE, phi_MLE)
  psi_MLE        <- get_psi(theta_MLE, weights)
  
  psi_MLE_SE <- .compute_psi_MLE_SE(
    theta_MLE,
    phi_MLE,
    weights,
    data
  )
  
  search_interval <- psi_MLE + c(-1, 1) * num_std_errors * psi_MLE_SE
  
  required_branch_alpha <- .compute_required_branch_alpha(
    R     = num_branches,
    alpha = global_alpha
  )
  
  crit <- 0.5 * stats::qchisq(1 - required_branch_alpha, df = 1)
  
  omega_hat_con_fn    <- .omega_hat_con_fn_closure(weights, psi_MLE)
  omega_hat_dist_fun  <- match.fun(config$omega$distribution)
  omega_hat_dist_args <- config$omega$args
  
  data_args <- list(
    psi_MLE       = psi_MLE,
    Y             = Y,
    t             = t,
    n_per_process = n_per_process,
    weights       = weights
  )
  
  likelihood_common <- list(
    p_cutoff  = p_cutoff,
    max_y_cap = max_y_cap
  )
  
  solver_args <- list(
    max_retries = max_retries,
    localsolver = localsolver,
    control     = list(
      xtol_rel = xtol_rel,
      maxeval  = maxeval
    )
  )
  
  result <- foreach::foreach(
    i               = seq_len(num_branches),
    .combine        = "list",
    .multicombine   = TRUE,
    .errorhandling  = "remove",
    .options.future = list(
      seed       = TRUE,
      chunk.size = chunk_size,
      packages   = c("nloptr", "dplyr", "matrixStats")
    )
  ) %dofuture% {
    omega_hat_init <- do.call(
      omega_hat_dist_fun,
      c(list(n = 2 * J), as.list(omega_hat_dist_args))
    )
    
    omega_hat <- .draw_omega_hat(
      omega_hat_con_fn = omega_hat_con_fn,
      init_guess       = omega_hat_init,
      localsolver      = localsolver,
      control          = list(
        xtol_rel = xtol_rel,
        maxeval  = maxeval
      )
    )
    
    likelihood_args <- c(
      list(theta_phi_0 = omega_hat),
      likelihood_common
    )
    
    eval_psi_fun <- .make_eval_psi_fun(
      data_args       = data_args,
      likelihood_args = likelihood_args,
      solver_args     = solver_args
    )
    
    branch_cutoff <- branch_params$max - crit
    
    branch_args <- list(
      increment     = increment,
      left_start    = adjacent$left,
      right_start   = adjacent$right,
      branch_cutoff = branch_cutoff,
      init_guess    = branch_params$argmax$theta_phi,
      max_retries   = solver_args$max_retries,
      psi_MLE       = data_args$psi_MLE
    )
    
    IL_branch <- .generate_branch(
      branch_args  = branch_args,
      eval_psi_fun = eval_psi_fun
    )
    
    list(
      IL_branch = IL_branch,
      omega_hat = omega_hat
    )
  }
  
  IL_branches <- lapply(result, `[[`, "IL_branch")
  omega_hats  <- lapply(result, `[[`, "omega_hat")
  
  log_L_bar <- .get_log_L_bar(IL_branches)
  
  list(
    log_L_bar   = log_L_bar,
    IL_branches = IL_branches,
    omega_hats  = omega_hats
  )
}
