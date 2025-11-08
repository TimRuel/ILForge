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
.build_psi_grid <- function(psi_min, psi_max, psi_MLE, step_size, fine_step_size, fine_window) {
  if (psi_min > psi_max) {
    tmp <- psi_min; psi_min <- psi_max; psi_max <- tmp
  }
  # coarse grid outside fine window
  left_coarse  <- seq(psi_min, psi_MLE - fine_window, by = step_size)
  right_coarse <- seq(psi_MLE + fine_window, psi_max, by = step_size)
  # fine grid inside window (include endpoints to stitch cleanly)
  left_fine  <- seq(max(psi_min, psi_MLE - fine_window), psi_MLE, by = fine_step_size)
  right_fine <- seq(psi_MLE, min(psi_max, psi_MLE + fine_window), by = fine_step_size)
  grid <- sort(unique(c(left_coarse, left_fine, right_fine, right_coarse)))
  grid[grid >= psi_min & grid <= psi_max]
}

#' @keywords internal
.eval_branch_ll_at_psi <- function(omega_hat, psi, Y, t, n_per_process, weights,
                                   init_guess, ratio_threshold, n_mc, p_cutoff, max_y_cap,
                                   localsolver, xtol_rel, maxeval, return_full,
                                   theta_MLE, phi_MLE, fix_phi) {
  J <- length(theta_MLE)
  if (is.null(init_guess) || length(init_guess) != 2L * J) {
    init_guess <- c(as.numeric(theta_MLE), as.numeric(phi_MLE))
  }
  
  opt_out <- get_theta_phi_hat(
    omega           = omega_hat,
    t               = t,
    n_per_process   = n_per_process,
    init_guess      = init_guess,
    weights         = weights,
    psi             = psi,
    ratio_threshold = ratio_threshold,
    n_mc            = n_mc,
    p_cutoff        = p_cutoff,
    max_y_cap       = max_y_cap,
    localsolver     = localsolver,
    xtol_rel        = xtol_rel,
    maxeval         = maxeval,
    return_full     = return_full
  )
  
  if (is.list(opt_out) && !is.null(opt_out$theta) && !is.null(opt_out$phi)) {
    theta_hat <- as.numeric(opt_out$theta)
    phi_hat   <- as.numeric(opt_out$phi)
  } else {
    raw <- if (is.list(opt_out) && !is.null(opt_out$par)) as.numeric(opt_out$par) else as.numeric(opt_out)
    theta_hat <- raw[seq_len(J)]
    phi_hat   <- raw[J + seq_len(J)]
  }
  
  phi_for_ll <- if (isTRUE(fix_phi)) as.numeric(phi_MLE) else phi_hat
  
  ll <- log_likelihood(
    theta         = theta_hat,
    phi           = phi_for_ll,
    Y             = Y,
    t             = t,
    n_per_process = n_per_process
  )
  
  list(ll = ll, init_guess_next = c(theta_hat, phi_hat))
}

# -------------------------------------------------------------------------
# One-Side Branch Evaluation
# -------------------------------------------------------------------------

#' Evaluate one side of an integrated likelihood branch
#'
#' Computes the integrated (profile) log-likelihood \eqn{\tilde L(\psi)} along a
#' grid of \eqn{\psi} values **either to the left or right** of \eqn{\psi_{\text{MLE}}}.
#'
#' @param direction "left" or "right"
#' @param omega_hat Numeric vector; optimized \eqn{\omegâ}.
#' @param weights Numeric vector of process weights (length J).
#' @param Y Numeric vector of counts.
#' @param t Numeric vector of exposures.
#' @param n_per_process Integer vector (# obs per process, length J).
#' @param theta_MLE,phi_MLE MLE anchors (for init & optional fixed-phi plug-in).
#' @param psi_MLE Scalar ψ at MLE.
#' @param psi_grid Numeric vector of ψ values (full grid; this picks a side).
#' @param ratio_threshold,n_mc,p_cutoff,max_y_cap,localsolver,xtol_rel,maxeval,return_full Solver options.
#' @param fix_phi Logical; if TRUE, hold \eqn{\phi=\phi_{\text{MLE}}} when evaluating.
#' @param trace Logical; verbose per-point messages.
#' @return list with `psi` and `Integrated` vectors
#' @export
run_branch_side <- function(direction,
                            omega_hat,
                            weights,
                            Y,
                            t,
                            n_per_process,
                            theta_MLE,
                            phi_MLE,
                            psi_MLE,
                            psi_grid,
                            ratio_threshold = 100,
                            n_mc = 1e5,
                            p_cutoff = 1e-12,
                            max_y_cap = 1e6,
                            localsolver = "SLSQP",
                            xtol_rel = 1e-8,
                            maxeval = 1000,
                            return_full = FALSE,
                            fix_phi = TRUE,
                            trace = FALSE) {
  
  direction <- match.arg(direction, choices = c("left", "right"))
  
  if (direction == "left") {
    psi_working <- rev(psi_grid[psi_grid <= psi_MLE])
  } else {
    psi_working <- psi_grid[psi_grid >  psi_MLE]
  }
  if (length(psi_working) == 0L) {
    return(list(psi = numeric(0), Integrated = numeric(0)))
  }
  
  J <- length(theta_MLE)
  init_guess <- c(as.numeric(theta_MLE), as.numeric(phi_MLE))
  if (!is.numeric(init_guess) || length(init_guess) != 2L * J) {
    stop("init_guess must be numeric length 2*J (theta then phi).")
  }
  
  ll_vals <- numeric(length(psi_working))
  
  for (i in seq_along(psi_working)) {
    psi_i <- psi_working[i]
    
    eval <- .eval_branch_ll_at_psi(
      omega_hat       = omega_hat,
      psi             = psi_i,
      Y               = Y,
      t               = t,
      n_per_process   = n_per_process,
      weights         = weights,
      init_guess      = init_guess,
      ratio_threshold = ratio_threshold,
      n_mc            = n_mc,
      p_cutoff        = p_cutoff,
      max_y_cap       = max_y_cap,
      localsolver     = localsolver,
      xtol_rel        = xtol_rel,
      maxeval         = maxeval,
      return_full     = return_full,
      theta_MLE       = theta_MLE,
      phi_MLE         = phi_MLE,
      fix_phi         = fix_phi
    )
    
    ll_vals[i] <- eval$ll
    init_guess <- eval$init_guess_next
    
    if (isTRUE(trace)) {
      message(sprintf("[Branch %s] psi=%.6g | ll=%.6f | fix_phi=%s",
                      direction, psi_i, ll_vals[i],
                      toupper(as.character(isTRUE(fix_phi)))))
    }
  }
  
  if (direction == "left") {
    list(psi = rev(psi_working), Integrated = rev(ll_vals))
  } else {
    list(psi = psi_working, Integrated = ll_vals)
  }
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
compute_IL_branch <- function(omega_hat, weights, Y, t, n_per_process,
                              theta_MLE, phi_MLE, psi_MLE, psi_grid,
                              ...) {
  
  left  <- run_branch_side("left",  omega_hat, weights, Y, t, n_per_process,
                           theta_MLE, phi_MLE, psi_MLE, psi_grid, ...)
  right <- run_branch_side("right", omega_hat, weights, Y, t, n_per_process,
                           theta_MLE, phi_MLE, psi_MLE, psi_grid, ...)
  
  df <- dplyr::bind_rows(
    data.frame(psi = left$psi,  Integrated = left$Integrated,  row.names = NULL),
    data.frame(psi = right$psi, Integrated = right$Integrated, row.names = NULL)
  )
  df <- dplyr::arrange(df, psi)
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
get_log_L_bar <- function(IL_branches) {
  if (length(IL_branches) == 0L) {
    return(list(df = data.frame(psi = numeric(0), Integrated = numeric(0)),
                branches_matrix = matrix(numeric(0), nrow = 0, ncol = 0)))
  }
  merged <- Reduce(function(x, y) dplyr::full_join(x, y, by = "psi"), IL_branches)
  merged <- dplyr::arrange(merged, psi)
  
  mat <- as.matrix(merged[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  
  eff_R <- pmax(rowSums(!is.na(mat)), 1L)
  lse   <- matrixStats::rowLogSumExps(mat, na.rm = TRUE)
  log_L_bar <- lse - log(eff_R)
  
  rownames(mat) <- as.character(merged$psi)
  
  list(
    df = data.frame(psi = merged$psi, Integrated = as.numeric(log_L_bar)),
    branches_matrix = t(mat)
  )
}

# -------------------------------------------------------------------------
# Top-level: compute integrated log-likelihood across many branches
# -------------------------------------------------------------------------

#' Compute integrated log-likelihood across multiple random branches (parallel)
#'
#' Strategy:
#' 1) Sample several ω̂ (no ψ search).
#' 2) Estimate curvature of branch logL at ψ_MLE using finite differences (±δ).
#' 3) Use spread of SEs across a few ω̂ to size the global ψ grid.
#' 4) Evaluate full branches on that grid and average via log-sum-exp.
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
  
  # Defaults
  num_std_errors <- num_std_errors %||% 3
  step_size      <- step_size      %||% 0.25
  fine_step_size <- fine_step_size %||% 0.05
  fine_window    <- fine_window    %||% 1.0
  delta_curv     <- delta_curv     %||% fine_step_size  # curvature step
  
  chunk_size     <- as.integer(chunk_size  %||% 1L)
  num_workers    <- as.integer(num_workers %||% 1L)
  num_branches   <- max(1L, chunk_size * num_workers)
  
  ratio_threshold <- ratio_threshold %||% 100
  n_mc            <- n_mc           %||% 1e5
  p_cutoff        <- p_cutoff       %||% 1e-12
  max_y_cap       <- max_y_cap      %||% 1e6
  localsolver     <- localsolver    %||% "SLSQP"
  xtol_rel        <- xtol_rel       %||% 1e-8
  maxeval         <- maxeval        %||% 1000L
  return_full     <- isTRUE(return_full)
  fix_phi         <- isTRUE(fix_phi)
  trace           <- isTRUE(trace)
  
  Y <- as.numeric(data$Y)
  t <- as.numeric(data$t)
  if (is.null(data$process)) stop("`data$process` is required to compute n_per_process.")
  n_per_process <- as.integer(table(data$process))
  
  # MLE anchors
  model     <- fit_model(data)
  theta_MLE <- get_theta_MLE(model)
  phi_MLE   <- get_phi_MLE(model)
  psi_MLE   <- get_psi(theta_MLE, weights)
  J         <- length(theta_MLE)
  
  # ω̂ constraint closure at ψ_MLE
  omega_hat_con_fn <- omega_hat_con_fn_closure(weights, psi_MLE)
  
  # 1) Draw ω̂ for multiple branches
  branch_omegas <- foreach::foreach(
    i = seq_len(num_branches),
    .combine = "list",
    .multicombine = TRUE,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE, chunk.size = chunk_size, packages = c("nloptr"))
  ) %dofuture% {
    init_guess <- stats::rexp(2L * J, rate = 1)
    omega_hat  <- get_omega_hat(omega_hat_con_fn, init_guess,
                                localsolver = localsolver,
                                xtol_rel = xtol_rel,
                                maxeval = maxeval)
    list(omega_hat = omega_hat)
  }
  omega_hats <- lapply(branch_omegas, `[[`, "omega_hat")
  
  # 2) Estimate curvature at ψ_MLE across a few ω̂ to size ψ-range
  k_sample <- min(length(omega_hats), 5L)
  if (k_sample == 0L) stop("Failed to compute any omega_hat branches.")
  
  se_vec <- numeric(k_sample)
  for (k in seq_len(k_sample)) {
    omega_hat_k <- omega_hats[[k]]
    # finite-difference second derivative L'' at ψ_MLE
    eval0 <- .eval_branch_ll_at_psi(omega_hat_k, psi_MLE,
                                    Y, t, n_per_process, weights,
                                    init_guess = c(theta_MLE, phi_MLE),
                                    ratio_threshold, n_mc, p_cutoff, max_y_cap,
                                    localsolver, xtol_rel, maxeval, return_full,
                                    theta_MLE, phi_MLE, fix_phi)
    evalm <- .eval_branch_ll_at_psi(omega_hat_k, psi_MLE - delta_curv,
                                    Y, t, n_per_process, weights,
                                    init_guess = eval0$init_guess_next,
                                    ratio_threshold, n_mc, p_cutoff, max_y_cap,
                                    localsolver, xtol_rel, maxeval, return_full,
                                    theta_MLE, phi_MLE, fix_phi)
    evalp <- .eval_branch_ll_at_psi(omega_hat_k, psi_MLE + delta_curv,
                                    Y, t, n_per_process, weights,
                                    init_guess = evalm$init_guess_next,
                                    ratio_threshold, n_mc, p_cutoff, max_y_cap,
                                    localsolver, xtol_rel, maxeval, return_full,
                                    theta_MLE, phi_MLE, fix_phi)
    
    Lpp <- (evalp$ll - 2 * eval0$ll + evalm$ll) / (delta_curv^2)
    # guard against non-concavity / numerical artifacts
    if (!is.finite(Lpp) || Lpp >= 0) {
      se_vec[k] <- NA_real_
    } else {
      se_vec[k] <- sqrt(1 / (-Lpp))
    }
  }
  # robust SE across sampled branches
  se_vec <- se_vec[is.finite(se_vec)]
  psi_se <- if (length(se_vec)) stats::median(se_vec) else fine_window
  
  # 3) Build global ψ-grid endpoints
  range_ext <- max(num_std_errors * psi_se, fine_window)
  psi_min <- psi_MLE - range_ext
  psi_max <- psi_MLE + range_ext
  
  psi_grid <- .build_psi_grid(psi_min, psi_max, psi_MLE, step_size, fine_step_size, fine_window)
  
  # 4) Evaluate full branches on the grid
  result <- foreach::foreach(
    omega_hat = omega_hats,
    .combine        = "list",
    .multicombine   = TRUE,
    .errorhandling  = "remove",
    .options.future = list(seed = TRUE, chunk.size = chunk_size, packages = c("nloptr","dplyr","matrixStats"))
  ) %dofuture% {
    
    IL_branch <- compute_IL_branch(
      omega_hat       = omega_hat,
      weights         = weights,
      Y               = Y,
      t               = t,
      n_per_process   = n_per_process,
      theta_MLE       = theta_MLE,
      phi_MLE         = phi_MLE,
      psi_MLE         = psi_MLE,
      psi_grid        = psi_grid,
      ratio_threshold = ratio_threshold,
      n_mc            = n_mc,
      p_cutoff        = p_cutoff,
      max_y_cap       = max_y_cap,
      localsolver     = localsolver,
      xtol_rel        = xtol_rel,
      maxeval         = maxeval,
      return_full     = return_full,
      fix_phi         = fix_phi,
      trace           = trace
    )
    
    list(IL_branch = IL_branch, omega_hat = omega_hat)
  }
  
  IL_branches <- lapply(result, `[[`, "IL_branch")
  omega_hat   <- lapply(result, `[[`, "omega_hat")
  
  merged <- get_log_L_bar(IL_branches)
  
  list(
    log_L_bar_df = merged$df,
    IL_branches  = IL_branches,
    omega_hat    = omega_hat
  )
}
