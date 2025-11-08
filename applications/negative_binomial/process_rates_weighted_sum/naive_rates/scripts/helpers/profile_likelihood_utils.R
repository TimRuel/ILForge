# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/profile_likelihood_utils.R

# -------------------------------------------------------------------------
# Profile Likelihood Utilities 
# -------------------------------------------------------------------------

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# -------------------------------------------------------------------------
# One-side profile likelihood branch
# -------------------------------------------------------------------------

#' Compute one side of the profile likelihood branch
#'
#' Evaluates the profile log-likelihood \eqn{\ell_p(\psi)} by moving outward
#' from the MLE in one direction (left or right) until the likelihood drops
#' below the chi-square cutoff corresponding to a desired confidence level.
#'
#' The method iteratively solves the constrained optimization problem for
#' \eqn{(\theta, \phi)} given \eqn{\psi} until the profile log-likelihood
#' falls below
#' \deqn{
#'   \ell_p(\psi_{\mathrm{MLE}}) - \tfrac{1}{2}\chi^2_{1, 1 - \alpha},
#' }
#' where \eqn{\alpha} is the nominal confidence level.
#'
#' @param direction Character scalar: `"left"` or `"right"` side of
#'   \eqn{\psi_{\mathrm{MLE}}}.
#' @param data Data frame with columns `Y`, `t`, and `process`.
#' @param weights Numeric vector of process weights (length = number of processes).
#' @param step_size Coarse grid step size for ψ.
#' @param fine_step_size Step size within a neighborhood of ψ_MLE.
#' @param fine_window Radius (in ψ units) defining that neighborhood.
#' @param alpha Confidence level (e.g., 0.05 for 95% CI).
#' @param trace Logical; if TRUE, print iteration messages.
#'
#' @details
#' The search proceeds from \eqn{\psi_{\text{MLE}}} outward in the chosen
#' direction. Near the MLE (within `fine_window`), smaller steps of size
#' `fine_step_size` are used for smoothness; otherwise, `step_size` is used.
#' The branch terminates once the profile log-likelihood drops below
#' the chi-square threshold.
#'
#' @return A list containing:
#' \describe{
#'   \item{psi}{Numeric vector of ψ values.}
#'   \item{Profile}{Numeric vector of profile log-likelihoods.}
#'   \item{psi_MLE}{MLE value of ψ.}
#'   \item{log_L_p_max}{Maximum log-likelihood value at ψ_MLE.}
#'   \item{direction}{Branch direction ("left" or "right").}
#' }
#'
#' @seealso [get_profile_LL()]
#' @export
compute_profile_branch <- function(direction,
                                   data,
                                   weights,
                                   step_size,
                                   fine_step_size,
                                   fine_window,
                                   alpha,
                                   trace = FALSE) {
  
  direction <- match.arg(direction, c("left", "right"))
  
  Y <- data$Y
  t <- data$t
  n_per_process <- table(data$process)
  J <- length(n_per_process)
  
  model <- fit_model(data)
  
  # MLE anchors
  theta_MLE <- get_theta_MLE(model)
  phi_MLE   <- get_phi_MLE(model)
  psi_MLE   <- get_psi(theta_MLE, weights)
  log_L_p_max <- log_likelihood(theta_MLE, phi_MLE, Y, t, n_per_process)
  
  # Chi-square cutoff (1 df, one-sided)
  crit <- qchisq(1 - alpha, df = 1) / 2
  stopping_val <- log_L_p_max - crit
  
  # Initialize state
  psi_vals <- numeric(0)
  log_vals <- numeric(0)
  init_guess <- c(theta_MLE, phi_MLE)
  
  # Direction setup
  sign <- if (direction == "left") -1 else 1
  psi <- psi_MLE + sign * fine_step_size
  log_L_p <- log_L_p_max
  
  # Iteratively step until log-likelihood drops below threshold
  while (log_L_p >= stopping_val) {
    
    opt_out <- get_theta_phi_hat(
      omega           = list(theta = theta_MLE, phi = phi_MLE),
      t               = t,
      n_per_process   = n_per_process,
      init_guess      = init_guess,
      weights         = weights,
      psi             = psi
    )
    
    # Robust parameter extraction
    if (is.list(opt_out) && !is.null(opt_out$theta) && !is.null(opt_out$phi)) {
      theta_hat <- as.numeric(opt_out$theta)
      phi_hat   <- as.numeric(opt_out$phi)
    } else {
      raw <- if (is.list(opt_out) && !is.null(opt_out$par)) as.numeric(opt_out$par) else as.numeric(opt_out)
      if (length(raw) < 2L * J)
        stop("get_theta_phi_hat() returned a parameter vector of insufficient length.")
      theta_hat <- raw[seq_len(J)]
      phi_hat   <- raw[J + seq_len(J)]
    }
    
    init_guess <- c(theta_hat, phi_hat)
    log_L_p <- log_likelihood(theta_hat, phi_hat, Y, t, n_per_process)
    
    psi_vals <- c(psi_vals, psi)
    log_vals <- c(log_vals, log_L_p)
    
    if (isTRUE(trace)) {
      message(sprintf("[Profile %s] ψ = %.6g | logL = %.6f", direction, psi, log_L_p))
    }
    
    # Fine vs. coarse step selection
    step <- if (abs(psi - psi_MLE) <= fine_window) fine_step_size else step_size
    psi <- psi + sign * step
  }
  
  # Optional interpolation for smooth cutoff
  if (length(log_vals) >= 2) {
    last_idx <- length(log_vals)
    if (log_vals[last_idx] < stopping_val) {
      psi_prev <- psi_vals[last_idx - 1]
      log_prev <- log_vals[last_idx - 1]
      psi_last <- psi_vals[last_idx]
      log_last <- log_vals[last_idx]
      slope <- (log_last - log_prev) / (psi_last - psi_prev)
      psi_cutoff <- psi_prev + (stopping_val - log_prev) / slope
      psi_vals <- c(psi_vals, psi_cutoff)
      log_vals <- c(log_vals, stopping_val)
    }
  }
  
  # Attach ψ_MLE point to left branch if needed
  if (direction == "left") {
    psi_vals <- c(rev(psi_vals), psi_MLE)
    log_vals <- c(rev(log_vals), log_L_p_max)
  }
  
  list(
    psi = psi_vals,
    Profile = log_vals,
    psi_MLE = psi_MLE,
    log_L_p_max = log_L_p_max,
    direction = direction
  )
}

# -------------------------------------------------------------------------
# Full profile likelihood (left + right) and merging
# -------------------------------------------------------------------------

#' Compute full profile likelihood curve for ψ
#'
#' Combines left and right profile branches around \eqn{\psi_{\text{MLE}}}
#' to form a full profile likelihood curve.
#'
#' @param config Experiment configuration list, must contain
#'   `optimization_specs$PL`.
#' @param data Data frame with columns `Y`, `t`, and `process`.
#' @param weights Numeric vector of process weights.
#'
#' @details
#' Runs [compute_profile_branch()] for both `"left"` and `"right"` directions
#' and merges the resulting curves into a single, sorted ψ-grid.
#' The 95% (or other α-level) confidence interval can be extracted by finding
#' the ψ values where:
#' \deqn{
#'   \ell_p(\psi) = \ell_p(\psi_{\mathrm{MLE}}) - \tfrac{1}{2}\chi^2_{1,1-\alpha}.
#' }
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{psi}{Grid of ψ values (sorted ascending).}
#'   \item{Profile}{Corresponding profile log-likelihood values.}
#' }
#'
#' @seealso [compute_profile_branch()]
#' @export
get_profile_LL <- function(config, data, weights) {
  
  pl_cfg <- config$optimization_specs$PL %||% list()
  invisible(list2env(pl_cfg, environment()))
  
  # Default parameters
  step_size      <- step_size      %||% 0.25
  fine_step_size <- fine_step_size %||% 0.05
  fine_window    <- fine_window    %||% 1.0
  alpha_levels   <- alpha_levels   %||% c(0.05)
  alpha          <- min(alpha_levels)
  
  # Precompute MLE model once to avoid duplication across branches
  model <- fit_model(data)
  
  theta_MLE <- get_theta_MLE(model)
  phi_MLE   <- get_phi_MLE(model)
  psi_MLE   <- get_psi(theta_MLE, weights)
  log_L_p_max <- log_likelihood(theta_MLE, phi_MLE, data$Y, data$t, table(data$process))
  
  res <- foreach::foreach(
    dir = c("left", "right"),
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = 2,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE, chunk.size = 1, packages = "nloptr")
  ) %dofuture% {
    compute_profile_branch(
      direction      = dir,
      data           = data,
      weights        = weights,
      step_size      = step_size,
      fine_step_size = fine_step_size,
      fine_window    = fine_window,
      alpha          = alpha,
      trace          = trace %||% FALSE
    )
  }
  
  # Combine left/right results
  df <- do.call(
    rbind,
    lapply(res, function(x) data.frame(psi = x$psi, Profile = x$Profile))
  )
  
  df <- df[order(df$psi), , drop = FALSE]
  rownames(df) <- NULL
  df
}
