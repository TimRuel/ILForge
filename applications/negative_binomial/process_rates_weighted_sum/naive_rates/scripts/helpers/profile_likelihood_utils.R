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
#' below the chi-square cutoff for a given confidence level.
#'
#' @param direction `"left"` or `"right"` side of \eqn{\psi_{\text{MLE}}}
#' @param data Data frame with at least columns `Y`, `t`, and `process`
#' @param weights Numeric process weights (length = number of processes)
#' @param step_size Coarse grid step size for ψ
#' @param fine_step_size Step size in a neighborhood around ψ_MLE
#' @param fine_window Radius around ψ_MLE where fine steps are used
#' @param alpha Confidence level (e.g., 0.05 gives 95% CI)
#' @param trace Logical; print iteration messages
#'
#' @return List with vectors:
#' \describe{
#'   \item{psi}{ψ values}
#'   \item{Profile}{profile log-likelihood values}
#' }
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
  
  model <- fit_model(data)
  
  # MLE anchors
  theta_MLE <- get_theta_MLE(model)
  phi_MLE   <- get_phi_MLE(model)
  psi_MLE   <- get_psi(theta_MLE, weights)
  log_L_p_max <- log_likelihood(theta_MLE, phi_MLE, Y, t, n_per_process)
  
  # Chi-square cutoff for 1 df
  crit <- qchisq(1 - alpha, df = 1) / 2
  stopping_val <- log_L_p_max - crit
  
  psi_vals <- numeric(0)
  log_vals <- numeric(0)
  
  # Initial step direction and anchor
  if (direction == "left") {
    step_anchor <- floor(psi_MLE / step_size) * step_size
    sign <- -1
  } else {
    step_anchor <- ceiling(psi_MLE / step_size) * step_size
    sign <- 1
  }
  
  psi <- psi_MLE + sign * fine_step_size

  log_L_p <- log_L_p_max
  
  init_guess <- c(theta_MLE, phi_MLE)
  
  while (log_L_p >= stopping_val) {
    
    # Solve constrained problem for (theta, phi) given (omega_hat, psi)
    opt_out <- get_theta_phi_hat(
      omega           = list(theta = theta_MLE, phi = phi_MLE),
      t               = t,
      n_per_process   = n_per_process,
      init_guess      = init_guess,
      weights         = weights,
      psi             = psi
      # ratio_threshold = ratio_threshold,
      # n_mc            = n_mc,
      # p_cutoff        = p_cutoff,
      # max_y_cap       = max_y_cap,
      # localsolver     = localsolver,
      # xtol_rel        = xtol_rel,
      # maxeval         = maxeval,
      # return_full     = return_full
    )
    
    # Extract (theta_hat, phi_hat) robustly from various return shapes
    theta_hat <- NULL
    phi_hat   <- NULL
    
    if (is.list(opt_out) && !is.null(opt_out$theta) && !is.null(opt_out$phi)) {
      theta_hat <- as.numeric(opt_out$theta)
      phi_hat   <- as.numeric(opt_out$phi)
      init_guess <- c(theta_hat, phi_hat)
    } else {
      # Assume raw parameter vector (theta, phi) possibly under $par
      raw <- if (is.list(opt_out) && !is.null(opt_out$par)) as.numeric(opt_out$par) else as.numeric(opt_out)
      if (length(raw) < 2L * J) {
        stop("get_theta_phi_hat() returned a parameter vector of insufficient length.")
      }
      theta_hat <- raw[seq_len(J)]
      phi_hat   <- raw[J + seq_len(J)]
      init_guess <- raw
    }
    
    log_L_p <- log_likelihood(theta_hat, phi_hat, Y, t, n_per_process)
    
    psi_vals <- c(psi_vals, psi)
    log_vals <- c(log_vals, log_L_p)

    if (isTRUE(trace)) {
      message(sprintf("[Profile %s] psi=%.6g | logL=%.6f", direction, psi, log_L_p))
    }
    
    # Choose fine / coarse increment
    dist <- abs(psi - psi_MLE)
    fine <- dist < fine_window ||
      (direction == "left"  && psi > step_anchor) ||
      (direction == "right" && psi < step_anchor)
    
    step <- if (fine) fine_step_size else step_size
    psi <- psi + sign * step
  }
  
  # Attach MLE on left side
  if (direction == "left") {
    psi_vals <- c(rev(psi_vals), psi_MLE)
    log_vals <- c(rev(log_vals), log_L_p_max)
  }
  
  list(psi = psi_vals, Profile = log_vals)
}

# -------------------------------------------------------------------------
# Full profile likelihood (left + right) and merging
# -------------------------------------------------------------------------

#' Compute full profile likelihood curve for ψ
#'
#' @param config Experiment config, must contain `optimization_specs$PL`
#' @param data Data frame with `Y`, `t`, and `process`
#' @param weights Numeric vector of process weights
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{psi}{Grid of ψ values}
#'   \item{Profile}{Profile log-likelihood}
#' }
#' @export
get_profile_LL <- function(config, data, weights) {
  
  pl_cfg <- config$optimization_specs$PL %||% list()
  invisible(list2env(pl_cfg, environment()))
  
  # Defaults if not set in config
  step_size      <- step_size      %||% 0.25
  fine_step_size <- fine_step_size %||% 0.05
  fine_window    <- fine_window    %||% 1.0
  alpha_levels   <- alpha_levels   %||% c(0.05)
  alpha          <- min(alpha_levels)
  
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
  
  # Convert to data frame
  df <- do.call(
    rbind,
    lapply(res, function(x) data.frame(psi = x$psi, Profile = x$Profile))
  )
  
  df <- df[order(df$psi), , drop = FALSE]
  rownames(df) <- NULL
  df
}
