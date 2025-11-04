# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/integrated_likelihood_utils.R

# -------------------------------------------------------------------------
# Integrated Likelihood Utilities (Negative Binomial, weighted-sum psi)
# -------------------------------------------------------------------------

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# -------------------------------------------------------------------------
# Branch parameters
# -------------------------------------------------------------------------

compute_branch_params <- function(omega_hat_con_fn, init_guess, Y, t, 
                                  n_per_process, weights, search_interval,
                                  localsolver, xtol_rel, max_eval) {
  
  omega_hat <- get_omega_hat(omega_hat_con_fn, init_guess, 
                             localsolver = localsolver,
                             xtol_rel = xtol_rel,
                             maxeval = maxeval)
  
  branch <- function(psi) {
    
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
      return_full     = FALSE
    )
    
    theta_hat <- opt_out$theta_hat
    phi_hat <- opt_out$phi_hat
    
    return(log_likelihood(theta_hat, phi_hat, Y, t, n_per_process))
  }
  
  opt <- optimize(branch,
                  interval = search_interval,
                  maximum = TRUE,
                  tol = 0.1)
  
  psi_mode <- opt$maximum
  theta_phi_mode <- get_theta_phi_hat(
    omega           = omega_hat,
    t               = t,
    n_per_process   = n_per_process,
    init_guess      = init_guess,
    weights         = weights,
    psi             = psi_mode,
    ratio_threshold = ratio_threshold,
    n_mc            = n_mc,
    p_cutoff        = p_cutoff,
    max_y_cap       = max_y_cap,
    localsolver     = localsolver,
    xtol_rel        = xtol_rel,
    maxeval         = maxeval,
    return_full     = FALSE
  )
  
  list(psi_mode = psi_mode,
       theta_phi_mode = theta_phi_mode,
       omega_hat = omega_hat)
}

get_branch_params_list <- function(config, data, weights) {
  
  invisible(list2env(config$optimization_specs$IL, env = environment()))
  
  Y <- data$Y
  t <- data$t
  n_per_process <- table(data$process)
  model <- fit_model(config, data)
  search_interval <- get_search_interval(model, data, weights, num_std_errors)
  theta_MLE <- get_theta_MLE(model)
  phi_MLE <- get_phi_MLE(model)
  psi_MLE <- get_psi(theta_MLE, weights)
  
  num_branches <- chunk_size * num_workers
  
  foreach(
    i = 1:num_branches,
    .combine = "c",
    .multicombine = TRUE,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size,
                           packages = c("nloptr", "dplyr"))
  ) %dofuture% {
    
    init_guess <- rnorm(2*length(theta_MLE))
    list(compute_branch_params(omega_hat_con_fn, init_guess, Y, t, 
                               n_per_process, weights, search_interval,
                               localsolver, xtol_rel, max_eval))
  }
}

# -------------------------------------------------------------------------
# One-Side Branch Evaluation
# -------------------------------------------------------------------------

#' Evaluate one side of an integrated likelihood branch
#'
#' Computes the integrated (profile) log-likelihood \eqn{\tilde L(\psi)} along a
#' grid of \eqn{\psi} values **either to the left or right** of the MLE \eqn{\psi_{\text{MLE}}}.
#' For each grid point, it solves the constrained problem for \eqn{(\theta,\phi)}
#' given \eqn{\omega} and \eqn{\psi}, then plugs the resulting parameters into
#' the NB log-likelihood.
#'
#' @param direction Character scalar: `"left"` or `"right"`. Which side of \eqn{\psi_{\text{MLE}}} to evaluate.
#' @param omega_hat Numeric vector; optimized \eqn{\omegâ} (e.g., weights-to-parameters mapping).
#' @param weights Numeric vector of process weights (length = number of processes \eqn{J}).
#' @param Y Numeric vector of observed counts (concatenated across processes).
#' @param t Numeric vector of exposure times (same length/order as `Y`).
#' @param n_per_process Integer vector of observation counts per process (length = \eqn{J}).
#' @param theta_MLE Numeric vector of MLEs for \eqn{\theta} at the global MLE (length \eqn{J}).
#' @param phi_MLE Numeric vector of MLEs for \eqn{\phi} at the global MLE (length \eqn{J} or 1, depending on your model).
#' @param psi_MLE Numeric scalar, the \eqn{\psi} value at the global MLE.
#' @param psi_grid Numeric vector of candidate \eqn{\psi} values (full grid; this function selects the side).
#' @param ratio_threshold Numeric; control parameter passed to `get_theta_phi_hat()` (MC safety / truncation).
#' @param n_mc Integer; Monte Carlo sample size for `get_theta_phi_hat()` if it uses MC.
#' @param p_cutoff Numeric; small probability cutoff for MC tail trimming inside `get_theta_phi_hat()`.
#' @param max_y_cap Integer; cap for simulated counts to avoid extreme tails.
#' @param localsolver Character; local solver name for `get_theta_phi_hat()` (e.g., `"SLSQP"` for `nloptr`).
#' @param xtol_rel Numeric; relative tolerance for `get_theta_phi_hat()`.
#' @param maxeval Integer; max evaluations for `get_theta_phi_hat()`.
#' @param return_full Logical; forwarded to `get_theta_phi_hat()` if it supports returning extra diagnostics.
#' @param fix_phi Logical; if `TRUE` (default), hold \eqn{\phi=\phi_{\text{MLE}}} when computing \eqn{\tilde L};
#'   if `FALSE`, plug in the optimized \eqn{\phî(\psi)}.
#' @param trace Logical; if `TRUE`, prints per-grid-point progress and constraint info.
#'
#' @return A list with:
#' \describe{
#'   \item{psi}{Numeric vector of \eqn{\psi} values on the requested side, sorted ascending.}
#'   \item{Integrated}{Numeric vector of integrated (profile) log-likelihood values aligned with `psi`.}
#' }
#' @export
run_branch_side <- function(direction,
                            omega_hat, 
                            weights, 
                            Y, 
                            t, 
                            n_per_process,
                            theta_MLE,      # used for dimension J and init
                            phi_MLE,        # used for init and/or fixed-phi plug-in
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
  J <- length(theta_MLE)
  
  # Select side of the grid; keep *working* order so successive init guesses are warm-started logically.
  if (direction == "left") {
    psi_working <- rev(psi_grid[psi_grid <= psi_MLE])
  } else {
    psi_working <- psi_grid[psi_grid >  psi_MLE]
  }
  
  if (length(psi_working) == 0L) {
    return(list(psi = numeric(0), Integrated = numeric(0)))
  }
  
  # Initialize optimizer guess as concatenated (theta, phi)
  init_guess <- c(as.numeric(theta_MLE), as.numeric(phi_MLE))
  if (!is.numeric(init_guess) || length(init_guess) != 2L * J) {
    stop("init_guess must be a numeric vector of length 2*J (theta then phi).")
  }
  
  ll_vals <- numeric(length(psi_working))
  
  for (i in seq_along(psi_working)) {
    
    psi_i <- psi_working[i]
    
    # Solve constrained problem for (theta, phi) given (omega_hat, psi_i)
    opt_out <- get_theta_phi_hat(
      omega           = omega_hat,
      t               = t,
      n_per_process   = n_per_process,
      init_guess      = init_guess,
      weights         = weights,
      psi             = psi_i,
      ratio_threshold = ratio_threshold,
      n_mc            = n_mc,
      p_cutoff        = p_cutoff,
      max_y_cap       = max_y_cap,
      localsolver     = localsolver,
      xtol_rel        = xtol_rel,
      maxeval         = maxeval,
      return_full     = return_full
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
    
    # Decide which phi to plug into the log-likelihood
    phi_for_ll <- if (isTRUE(fix_phi)) as.numeric(phi_MLE) else phi_hat
    
    # Compute integrated (profile) log-likelihood at psi_i
    ll_vals[i] <- log_likelihood(
      theta         = theta_hat,
      phi           = phi_for_ll,
      Y             = Y,
      t             = t,
      n_per_process = n_per_process
    )
    
    if (isTRUE(trace)) {
      message(sprintf(
        "[Branch %s] psi=%.6g | ll=%.6f | fix_phi=%s",
        direction, psi_i, ll_vals[i],
        toupper(as.character(isTRUE(fix_phi)))
      ))
    }
  }
  
  # Return psi in ascending order with aligned ll
  if (direction == "left") {
    psi_out <- rev(psi_working)
    ll_out  <- rev(ll_vals)
  } else {
    psi_out <- psi_working
    ll_out  <- ll_vals
  }
  
  list(psi = psi_out, Integrated = ll_out)
}

# -------------------------------------------------------------------------
# Full Branch (Left + Right)
# -------------------------------------------------------------------------

#' Compute a full integrated likelihood branch (left + right)
#'
#' Runs \code{\link{run_branch_side}} for both sides around \eqn{\psi_{\text{MLE}}}
#' and concatenates the results, sorted by \eqn{\psi}.
#'
#' @inheritParams run_branch_side
#' @param ... Additional arguments forwarded to \code{\link{run_branch_side}}
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{psi}{Numeric, sorted ascending.}
#'   \item{Integrated}{Numeric, integrated (profile) log-likelihood at each \eqn{\psi}.}
#' }
#' @export
compute_IL_branch <- function(branch_params, weights, Y, t, n_per_process,
                              theta_MLE, phi_MLE, psi_MLE, psi_bar, psi_grid,
                              ...) {
  
  psi_mode <- branch_params$psi_mode
  theta_phi_hat_mode <- branch_params$theta_phi_mode
  omega_hat <- branch_params$omega_hat
  
  psi_grid_left <- rev(psi_grid[psi_mode <= psi_bar])
  psi_grid_right <- psi_grid[psi_mode > psi_bar]
  
  left  = run_branch_side("left",  omega_hat, weights, Y, t, n_per_process,
                          theta_MLE, phi_MLE, psi_MLE, psi_grid_left, ...)
  right = run_branch_side("right", omega_hat, weights, Y, t, n_per_process,
                          theta_MLE, phi_MLE, psi_MLE, psi_grid_right, ...)
  
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
#' Given a list of branch data frames (each from \code{\link{compute_IL_branch}}),
#' aligns them on the union of \eqn{\psi} grids and computes the averaged
#' log-likelihood via a column-wise log-sum-exp:
#' \deqn{\bar L(\psi) = \frac{1}{R}\sum_{r=1}^R \exp\{\tilde L_r(\psi)\},}
#' where missing branch values at a given \eqn{\psi} are ignored in both the sum
#' and the divisor \eqn{R(\psi)}.
#'
#' @param IL_branches List of data frames produced by \code{\link{compute_IL_branch}}.
#'
#' @return A list with:
#' \describe{
#'   \item{df}{Data frame with columns `psi` and `Integrated` (the averaged log-likelihood).}
#'   \item{branches_matrix}{Numeric matrix of aligned branch log-likelihoods
#'         (rows = branches, columns = \eqn{\psi} in ascending order).}
#' }
#' @export
get_log_L_bar <- function(IL_branches) {
  
  if (length(IL_branches) == 0L) {
    return(list(df = data.frame(psi = numeric(0), Integrated = numeric(0)),
                branches_matrix = matrix(numeric(0), nrow = 0, ncol = 0)))
  }
  
  # Full outer join on psi to align columns across branches
  merged <- Reduce(function(x, y) dplyr::full_join(x, y, by = "psi"), IL_branches)
  merged <- dplyr::arrange(merged, psi)
  
  # Collect only the numeric columns after psi; each is a branch's Integrated values
  mat <- as.matrix(merged[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  
  # Column-wise log-sum-exp with NA-robust averaging
  # Effective R per column is number of non-NA entries
  eff_R <- pmax(colSums(!is.na(mat)), 1L)
  lse   <- matrixStats::colLogSumExps(mat, na.rm = TRUE)
  log_L_bar <- lse - log(eff_R)
  
  colnames(mat) <- as.character(merged$psi)
  
  list(
    df = data.frame(psi = merged$psi, Integrated = as.numeric(log_L_bar)),
    branches_matrix = t(mat) # rows=branches → transpose so rows=branches (as before)
  )
}

# -------------------------------------------------------------------------
# Top-level: compute integrated log-likelihood across many branches
# -------------------------------------------------------------------------

#' Compute integrated log-likelihood across multiple random branches (parallel)
#'
#' Builds a \eqn{\psi} grid from data and weights, samples multiple random
#' starting points for the \eqn{\omegâ} constraint solve, evaluates the full
#' left+right branch for each, and merges via \code{\link{get_log_L_bar}}.
#'
#' @param config List-like configuration. Expected to contain
#'   \code{config$optimization_specs$IL} with fields (all optional; sensible defaults used):
#'   \describe{
#'     \item{num_std_errors}{Numeric, width (in SE units) for coarse grid (default 3).}
#'     \item{step_size}{Numeric, coarse grid spacing (default 0.25).}
#'     \item{fine_step_size}{Numeric, fine grid spacing near \eqn{\psi_{\text{MLE}}} (default 0.05).}
#'     \item{fine_window}{Numeric, half-width of fine window (default 1.0).}
#'     \item{chunk_size}{Integer, number of branches per worker chunk (default 1).}
#'     \item{num_workers}{Integer, number of chunks (multiplies with chunk_size) (default 1).}
#'     \item{ratio_threshold, n_mc, p_cutoff, max_y_cap, localsolver, xtol_rel, maxeval, return_full, fix_phi}{Passed through to optimizers.}
#'     \item{trace}{Logical; print progress per grid point (default FALSE).}
#'   }
#' @param data Data frame with at least columns:
#'   \describe{
#'     \item{Y}{Observed counts.}
#'     \item{t}{Exposure times.}
#'     \item{process}{Factor or integer process ID used to build \code{n_per_process}.}
#'   }
#' @param weights Numeric vector of process weights (length = number of processes).
#'
#' @return A list with:
#' \describe{
#'   \item{log_L_bar_df}{Data frame with `psi` and averaged integrated log-likelihood.}
#'   \item{IL_branches}{List of branch data frames.}
#'   \item{omega_hat}{List of \eqn{\omegâ} solutions used per branch.}
#' }
#' @examples
#' \dontrun{
#' out <- get_integrated_LL(config, data = df, weights = w)
#' plot(out$log_L_bar_df$psi, out$log_L_bar_df$Integrated, type = "l")
#' }
#' @export
get_integrated_LL <- function(config, branch_params_list, data, weights) {
  
  # Pull IL options into env (if present), then provide defaults
  il_cfg <- config$optimization_specs$IL %||% list()
  invisible(list2env(il_cfg, envir = environment()))
  
  num_std_errors <- num_std_errors %||% 3
  step_size      <- step_size      %||% 0.25
  fine_step_size <- fine_step_size %||% 0.05
  fine_window    <- fine_window    %||% 1.0
  
  chunk_size     <- as.integer(chunk_size  %||% 1L)
  num_workers    <- as.integer(num_workers %||% 1L)
  
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
  
  model <- fit_model(data)
  
  # MLE anchors
  theta_MLE <- get_theta_MLE(model)
  phi_MLE   <- get_phi_MLE(model)
  psi_MLE   <- get_psi(theta_MLE, weights)
  
  psi_modes <- sapply(branch_params_list, `[[`, 1)
  theta_phi_modes <- lapply(branch_params_list, `[[`, 2)
  omega_hats <- lapply(branch_params_list, `[[`, 3)
  
  psi_bar <- mean(psi_modes)
  psi_bar_SE <- sd(psi_modes)
  psi_grid_endpoints <- psi_bar + c(-1, 1) * (max(abs(psi_modes - psi_bar)) + num_std_errors * psi_bar_SE)
  psi_grid <- seq(psi_grid_endpoints[1], psi_grid_endpoints[2], increment)
  
  Y <- as.numeric(data$Y)
  t <- as.numeric(data$t)
  if (is.null(data$process)) {
    stop("`data$process` is required to compute n_per_process.")
  }
  n_per_process <- as.integer(table(data$process))
  
  # Parallel branch evaluation
  result <- foreach::foreach(
    branch_params = branch_params_list,,
    .combine        = "c",
    .multicombine   = TRUE,
    .errorhandling  = "remove",
    .options.future = list(seed = TRUE, chunk.size = chunk_size, packages = c("nloptr", "dplyr", "matrixStats"))
  ) %dofuture% {
    
    # Simple reproducible random init for (theta, phi) given J
    J <- length(theta_MLE)
    init_guess <- stats::rexp(2L * J, rate = 1)
    
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
