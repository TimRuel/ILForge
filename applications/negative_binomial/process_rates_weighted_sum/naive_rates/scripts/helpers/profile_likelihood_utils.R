# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/profile_likelihood_utils.R

# -------------------------------------------------------------------------
# Profile Likelihood Utilities (Negative Binomial, Weighted-Sum ψ)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Full profile likelihood (left + right) and merging
# -------------------------------------------------------------------------

get_profile_LL <- function(config, data, weights) {
  pl_cfg <- config$optimization_specs$PL
  pl_cfg <- sanitize_config(pl_cfg)
  invisible(list2env(pl_cfg, environment()))
  
  alpha          <- min(alpha_levels)
  Y             <- data$Y
  t             <- data$t
  n_per_process <- table(data$process)
  J             <- length(n_per_process)
  
  model <- fit_model(data)
  
  theta_MLE      <- get_theta_MLE(model)
  phi_MLE        <- get_phi_MLE(model)
  theta_phi_MLE  <- list(theta = theta_MLE, phi = phi_MLE)
  psi_MLE        <- get_psi(theta_MLE, weights)
  
  crit <- 0.5 * stats::qchisq(1 - alpha, df = 1)
  
  data_args <- list(
    psi_MLE       = psi_MLE,
    Y             = Y,
    t             = t,
    n_per_process = n_per_process,
    weights       = weights
  )
  
  likelihood_args <- list(
    theta_phi_0 = theta_phi_MLE,
    p_cutoff    = p_cutoff,
    max_y_cap   = max_y_cap
  )
  
  solver_args <- list(
    max_retries = max_retries,
    localsolver = localsolver,
    control     = list(
      xtol_rel = xtol_rel,
      maxeval  = maxeval
    )
  )
  
  eval_psi_fun <- .make_eval_psi_fun(
    data_args       = data_args,
    likelihood_args = likelihood_args,
    solver_args     = solver_args
  )
  
  result <- foreach::foreach(
    i               = seq_len(2),
    .combine        = "list",
    .multicombine   = TRUE,
    .errorhandling  = "remove",
    .options.future = list(
      seed       = TRUE,
      chunk.size = 1,
      packages   = c("nloptr", "dplyr")
    )
  ) %dofuture% {
    
    PL_branch <- .generate_branch(branch_args, eval_psi_fun)
    
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
