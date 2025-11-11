# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/integrated_likelihood_utils.R

# -------------------------------------------------------------------------
# Integrated Likelihood Utilities (Negative Binomial, weighted-sum psi)
# -------------------------------------------------------------------------

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# -------------------------------------------------------------------------
# Small helpers (grid + one-point eval)
# -------------------------------------------------------------------------

#' @keywords internal
#' @title Compute nearest ψ grid points to a branch-specific mode
#'
#' @description
#' Internal helper that identifies the immediate left and right ψ grid points
#' surrounding a branch-specific ψ mode (`psi_hat_branch`), assuming a uniform
#' grid of the form:
#'
#' \deqn{\psi_k = \psi_{\text{MLE}} + k \cdot \text{increment}}
#'
#' where `k` is an integer grid index.
#'
#' This function does **not** construct the full grid. It computes the
#' nearest grid-aligned values algebraically, which is fast and avoids
#' unnecessary vector generation when extending branch evaluations.
#'
#' @param psi_hat_branch Numeric scalar. The ψ value at which the branch
#'   log-likelihood attains its maximum (branch-specific mode).
#' @param psi_MLE Numeric scalar. The global ψ MLE (acts as grid origin).
#' @param increment Positive numeric scalar. The grid spacing between ψ values.
#'
#' @return A named list with components:
#' \describe{
#'   \item{left}{Grid point immediately ≤ `psi_hat_branch`.}
#'   \item{right}{Grid point immediately ≥ `psi_hat_branch`.}
#'   \item{k_left}{Integer index of the left grid point relative to `psi_MLE`.}
#'   \item{k_right}{Integer index of the right grid point relative to `psi_MLE`.}
#' }
#'
#' @examples
#' .get_adjacent_grid_points(psi_hat_branch = 1.13, psi_MLE = 0, increment = 0.25)
#'
#' .get_adjacent_grid_points(psi_hat_branch = -0.68, psi_MLE = 0, increment = 0.25)
.get_adjacent_grid_points <- function(
    psi_hat_branch, 
    psi_MLE, 
    increment
    ) {
  
  if (!is.numeric(psi_hat_branch) || length(psi_hat_branch) != 1L)
    stop("`psi_hat_branch` must be a numeric scalar.")
  if (!is.numeric(psi_MLE) || length(psi_MLE) != 1L)
    stop("`psi_MLE` must be a numeric scalar.")
  if (!is.numeric(increment) || length(increment) != 1L || increment <= 0)
    stop("`increment` must be a positive numeric scalar.")
  
  # Grid index of branch_mode relative to psi_MLE
  k_float <- (psi_hat_branch - psi_MLE) / increment
  
  # Nearest integer grid indices on each side
  k_left  <- floor(k_float)
  k_right <- ceiling(k_float)
  
  list(
    left    = psi_MLE + k_left  * increment,
    right   = psi_MLE + k_right * increment,
    k_left  = k_left,
    k_right = k_right
  )
}

#' @keywords internal
#' @title Evaluate branch log-likelihood at a given ψ
#'
#' @description
#' Internal helper to evaluate the branch log-likelihood at a specific ψ value.
#' The function:
#'
#' 1. Constructs closure objects containing expected log-likelihood components
#'    conditional on `psi`.
#' 2. Solves for the maximizing (`theta`, `phi`) for that ψ.
#' 3. Computes the resulting branch log-likelihood.
#'
#' The function returns both the log-likelihood value and the updated parameter
#' vector, which is used as a warm start when stepping along the branch.
#'
#' @param psi Numeric scalar. ψ value at which to evaluate the branch.
#' @param init_guess Numeric vector. Current warm-start guess for (θ, φ).
#' @param omega_hat List with elements `theta` and `phi`, specifying the
#'   expectation-distribution parameters used in the integrated-likelihood
#'   evaluation.
#' @param Y Numeric vector of observed counts.
#' @param t Numeric vector of exposure times (one per observation).
#' @param n_per_process Integer vector giving number of observations per process.
#' @param weights Numeric vector of process weights used in the ψ definition.
#' @param ... Additional arguments passed through to the solver and closure
#'   constructors (e.g., numerical tolerances, hybrid expectation settings).
#'
#' @return A named list with:
#' \describe{
#'   \item{branch_val}{Log-likelihood value at the given ψ.}
#'   \item{theta_phi_hat}{Numeric vector of updated estimates of (θ, φ).}
#' }
.eval_branch_ll_at_psi <- function(
    psi, 
    init_guess, 
    omega_hat, 
    Y, 
    t, 
    n_per_process, 
    weights,
    p_cutoff,
    max_y_cap,
    ...
    ) {
  
  closures <- .make_theta_phi_hat_closures(
    psi           = psi, 
    omega_hat     = omega_hat, 
    t             = t, 
    n_per_process = n_per_process, 
    weights       = weights, 
    p_cutoff      = p_cutoff,
    max_y_cap     = max_y_cap
  )
  
  theta_phi_hat <- .get_theta_phi_hat(
    closures   = closures, 
    init_guess = init_guess, 
    ...
  )
  
  J <- length(omega_hat$theta)
  
  branch_val <- log_likelihood(
    theta         = theta_phi_hat[1:J],
    phi           = theta_phi_hat[(J+1):(2*J)],
    Y             = Y,
    t             = t,
    n_per_process = n_per_process
  )
  
  list(branch_val = branch_val, theta_phi_hat = theta_phi_hat)
}

# -------------------------------------------------------------------------
# One-Side Branch Evaluation
# -------------------------------------------------------------------------

#' @keywords internal
#' @title Evaluate one side of a branch outward from its mode
#'
#' @description
#' Starting from an already-identified branch mode (`branch_df$psi`), this
#' internal helper iteratively steps outward in ψ using the specified `increment`
#' (positive for right side, negative for left side). At each step it:
#'
#' 1. Solves for (`theta`, `phi`) conditional on the new ψ value.
#' 2. Computes the branch log-likelihood at that point.
#' 3. Appends the result to the branch data frame.
#'
#' The stepping continues until the branch log-likelihood falls below the
#' specified cutoff.
#'
#' Typically invoked twice from `.compute_IL_branch()`: once with positive
#' increment for the right side, once with negative increment for the left side.
#'
#' @param increment Numeric scalar. Step size and sign for ψ progression.
#'   Use a positive value to extend rightward, and a negative value to extend
#'   leftward.
#' @param branch_df Data frame containing one initial row with columns:
#'   \describe{
#'     \item{psi}{Branch mode ψ value.}
#'     \item{value}{Log-likelihood at the branch mode.}
#'   }
#' @param branch_cutoff Numeric scalar. Minimum log-likelihood threshold at
#'   which stepping stops. Usually equal to:
#'   \deqn{\ell_{\max} - \frac{1}{2}\chi^2_{1,\,\alpha}}
#'   for confidence interval construction.
#' @param init_guess Numeric vector containing (`theta`, `phi`) estimates at
#'   the branch mode. Used as the warm-start for numerical optimization as ψ
#'   steps outward.
#' @param omega_hat List with expectation parameters (`theta`, `phi`) used
#'   in integrated-likelihood computations.
#' @param Y Numeric vector of observed counts.
#' @param t Numeric vector of exposure times.
#' @param n_per_process Integer vector with number of observations per process.
#' @param weights Numeric process weights.
#' @param ... Additional arguments passed to `.eval_branch_ll_at_psi()`.
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{psi}{Evaluated ψ values (sorted ascending).}
#'   \item{value}{Corresponding branch log-likelihoods.}
#' }
#'
#' @examples
#' # Example (pseudocode):
#' # branch_df <- data.frame(psi = psi_mode, value = branch_max)
#' # out <- .run_branch_side(
#' #   increment     = 0.05,
#' #   branch_df     = branch_df,
#' #   branch_cutoff = branch_max - qchisq(0.99, df=1)/2,
#' #   init_guess    = theta_phi_at_mode,
#' #   omega_hat     = omega_hat,
#' #   Y = Y, t = t, n_per_process = n_per_process, weights = weights
#' # )
.run_branch_side <- function(
    increment,
    start,
    branch_cutoff,
    init_guess,
    omega_hat, 
    Y, 
    t, 
    n_per_process, 
    weights,
    p_cutoff,
    max_y_cap,
    ...
    ) {
  
  current_psi <- start
  current_val <- 1000

  branch_df <- data.frame(psi = numeric(0), value = numeric(0))
  
  # Step outward until cutoff reached
  while (is.finite(current_val) && current_val >= branch_cutoff) {
    
    eval <- .eval_branch_ll_at_psi(
      psi           = current_psi, 
      init_guess    = init_guess, 
      omega_hat     = omega_hat, 
      Y             = Y, 
      t             = t, 
      n_per_process = n_per_process, 
      weights       = weights, 
      p_cutoff      = p_cutoff,
      max_y_cap     = max_y_cap,
      ...
    )
    
    current_val <- eval$branch_val
    init_guess  <- eval$theta_phi_hat
    
    branch_df <- dplyr::bind_rows(branch_df, data.frame(psi = current_psi, value = current_val))
    
    current_psi <- current_psi + increment
  }
  
  # Ensure ascending ψ order (for consistent merging)
  branch_df <- branch_df |> 
    unique() |> 
    dplyr::arrange(psi)
  
  rownames(branch_df) <- NULL
  
  branch_df
}

# -------------------------------------------------------------------------
# Branch Maximization
# -------------------------------------------------------------------------

#' @keywords internal
#' @title Locate branch-specific ψ mode via bounded maximization
#'
.maximize_branch <- function(
    increment,
    search_interval,
    init_guess,
    omega_hat,
    Y,
    t,
    n_per_process,
    weights,
    p_cutoff,      
    max_y_cap,
    ...
    ) {

  # ---- Objective function for optimize() ----
  obj <- function(psi) {
    out <- .eval_branch_ll_at_psi(
      psi           = psi,
      init_guess    = init_guess,
      omega_hat     = omega_hat,
      Y             = Y,
      t             = t,
      n_per_process = n_per_process,
      weights       = weights,
      p_cutoff      = p_cutoff,
      max_y_cap     = max_y_cap,
      ...
    )
    -out$branch_val
  }
  
  psi_hat_branch <- optimize(obj, interval = search_interval)$minimum
  
  eval_at_mode <- .eval_branch_ll_at_psi(
    psi           = psi_hat_branch,
    init_guess    = init_guess,
    omega_hat     = omega_hat,
    Y             = Y,
    t             = t,
    n_per_process = n_per_process,
    weights       = weights,
    p_cutoff      = p_cutoff,
    max_y_cap     = max_y_cap,
    ...
  )
  
  list(
    psi_hat_branch = psi_hat_branch, 
    branch_max     = eval_at_mode$branch_val,
    theta_phi_hat  = eval_at_mode$theta_phi_hat
    )
}

# -------------------------------------------------------------------------
# Full Branch (Left + Right)
# -------------------------------------------------------------------------

#' Compute a full integrated likelihood branch (left + right)
#'
#' @inheritParams run_branch_side
#' @param psi_grid Numeric ψ grid (full).
#' @return data.frame with columns `psi`, `Integrated`
#' @export
.compute_IL_branch <- function(
    increment,
    alpha,
    omega_hat,
    theta_phi_MLE,
    psi_MLE,
    weights, 
    Y, 
    t, 
    n_per_process,
    search_interval,
    p_cutoff,
    max_y_cap,
    ...
    ) {
  
  branch_max_params <- .maximize_branch(
    increment       = increment,
    search_interval = search_interval,
    init_guess      = theta_phi_MLE,
    omega_hat       = omega_hat,
    Y               = Y,
    t               = t,
    n_per_process   = n_per_process,
    weights         = weights,
    p_cutoff        = p_cutoff,      
    max_y_cap       = max_y_cap,
    ...
  )
  
  psi_hat_branch <- branch_max_params$psi_hat_branch
  
  crit <- qchisq(1 - alpha, df = 1) / 2
  branch_cutoff <- branch_max_params$branch_max - crit
  
  init_guess <- branch_max_params$theta_phi_hat
  
  adjacent_grid_points <- .get_adjacent_grid_points(psi_hat_branch, psi_MLE, increment)
  
  left <- .run_branch_side(
    increment     = -increment,
    start         = adjacent_grid_points$left,
    branch_cutoff = branch_cutoff,
    init_guess    = init_guess,
    omega_hat     = omega_hat, 
    Y             = Y, 
    t             = t, 
    n_per_process = n_per_process, 
    weights       = weights,
    p_cutoff      = p_cutoff,      
    max_y_cap     = max_y_cap,
    ...
  )
  
  right <- .run_branch_side(
    increment     = increment,
    start         = adjacent_grid_points$right,
    branch_cutoff = branch_cutoff,
    init_guess    = init_guess,
    omega_hat     = omega_hat, 
    Y             = Y, 
    t             = t, 
    n_per_process = n_per_process, 
    weights       = weights,
    p_cutoff      = p_cutoff,      
    max_y_cap     = max_y_cap,
    ...
  )
  
  df <- left |> 
    dplyr::bind_rows(right) |> 
    dplyr::arrange(psi)
  
  rownames(df) <- NULL
  
  df
}

# -------------------------------------------------------------------------
# Merge branches via log-sum-exp
# -------------------------------------------------------------------------

#' Merge multiple integrated-likelihood branches via log-sum-exp
#'
#' @param IL_branches list of data.frames from `compute_IL_branch()`
#' @return list with `df` (psi, Integrated) and `branches_matrix`
#' @export
.get_log_L_bar <- function(IL_branches) {
  
  merged <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "psi"), IL_branches)
  merged <- dplyr::arrange(merged, psi)
  
  mat <- as.matrix(merged[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  
  eff_R <- ncol(mat)
  lse   <- matrixStats::rowLogSumExps(mat, na.rm = TRUE)
  log_L_bar <- lse - log(eff_R)
  
  rownames(mat) <- as.character(merged$psi)
  colnames(mat) <- NULL
  
  list(
    df = data.frame(psi = merged$psi, value = as.numeric(log_L_bar)),
    branches_matrix = mat
  )
}

# -------------------------------------------------------------------------
# Top-level: compute integrated log-likelihood across many branches
# -------------------------------------------------------------------------

#' Compute integrated log-likelihood across multiple random branches (parallel)
#'
#' @param config list with `optimization_specs$IL`
#' @param data data.frame with Y, t, process
#' @param weights numeric weights (length J)
#' @return list: `log_L_bar_df`, `IL_branches`, `omega_hat`
#' @export
get_integrated_LL <- function(config, data, weights) {
  
  il_cfg <- config$optimization_specs$IL %||% list()
  il_cfg <- sanitize_config(il_cfg)
  invisible(list2env(il_cfg, envir = environment()))
  
  alpha <- min(alpha_levels)
  num_branches <- max(1L, chunk_size * num_workers)
  
  Y <- data$Y
  t <- data$t
  n_per_process <- table(data$process)
  J <- length(n_per_process)
  
  model <- fit_model(data)
  
  theta_MLE <- get_theta_MLE(model)
  
  phi_MLE <- get_phi_MLE(model)
  
  theta_phi_MLE <- c(theta_MLE, phi_MLE)
  
  psi_MLE <- get_psi(theta_MLE, weights)
  
  psi_MLE_SE <- .compute_psi_MLE_SE(theta_MLE, phi_MLE, weights, data)
  
  search_interval <- psi_MLE + c(-1, 1) * num_std_errors * psi_MLE_SE
  
  omega_hat_con_fn <- .omega_hat_con_fn_closure(weights, psi_MLE)
  
  omega_hat_dist_fun <- match.fun(config$omega$distribution)
  omega_hat_dist_args <- config$omega$args
  
  result <- foreach::foreach(
    i               = seq_len(num_branches),
    .combine        = "list",
    .multicombine   = TRUE,
    .errorhandling  = "remove",
    .options.future = list(seed = TRUE, chunk.size = chunk_size, packages = c("nloptr","dplyr","matrixStats"))
  ) %dofuture% {
    
    omega_hat_init_guess <- do.call(omega_hat_dist_fun, c(list(n = 2*J), as.list(omega_hat_dist_args)))

    omega_hat <- .draw_omega_hat(
      omega_hat_con_fn = omega_hat_con_fn, 
      init_guess       = omega_hat_init_guess,
      localsolver      = localsolver,
      control          = list(xtol_rel = xtol_rel, maxeval = maxeval)
    )
    
    IL_branch <- .compute_IL_branch(
      increment       = increment,
      alpha           = alpha,
      omega_hat       = omega_hat,
      theta_phi_MLE   = theta_phi_MLE,
      psi_MLE         = psi_MLE,
      weights         = weights, 
      Y               = Y, 
      t               = t, 
      n_per_process   = n_per_process,
      search_interval = search_interval,
      p_cutoff        = p_cutoff,      
      max_y_cap       = max_y_cap,
      localsolver     = localsolver,
      control         = list(xtol_rel = xtol_rel, maxeval = maxeval)
    )
    
    list(IL_branch = IL_branch, omega_hat = omega_hat)
  }
  
  IL_branches <- lapply(result, `[[`, "IL_branch")
  omega_hats   <- lapply(result, `[[`, "omega_hat")
  
  log_L_bar <- .get_log_L_bar(IL_branches)
  
  list(
    log_L_bar   = log_L_bar,
    IL_branches = IL_branches,
    omega_hats  = omega_hats
  )
}

