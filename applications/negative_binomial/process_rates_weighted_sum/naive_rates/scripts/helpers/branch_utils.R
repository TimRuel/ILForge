# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/branch_utils.R

# ======================================================================
# One-sided branch extension (left or right)
# ======================================================================

.run_branch_side <- function(increment,
                             start,
                             branch_cutoff,
                             init_guess,
                             eval_psi_fun,
                             max_retries) {
  current_psi <- start
  current_val <- Inf
  current_par <- init_guess
  
  branch_df <- tibble::tibble(psi = numeric(), value = numeric())
  
  while (is.finite(current_val) && current_val >= branch_cutoff) {
    retry <- 0
    
    repeat {
      eval <- eval_psi_fun(current_psi, current_par)
      
      if (eval$branch_val <= current_val || retry >= max_retries) break
      
      retry <- retry + 1
      jitter_scale <- 0.1 * retry
      
      current_par <- current_par + stats::rnorm(
        n  = length(current_par),
        sd = jitter_scale
      )
    }
    
    if (eval$branch_val > current_val) {
      warning(
        sprintf(
          "Monotonicity violation at psi = %.4f after %d retries; resetting init_guess.",
          current_psi, retry
        ),
        call. = FALSE
      )
      
      eval <- eval_psi_fun(current_psi, current_par)
    }
    
    current_val <- eval$branch_val
    current_par <- eval$theta_phi_hat
    
    branch_df <- dplyr::bind_rows(
      branch_df,
      tibble::tibble(psi = current_psi, value = current_val)
    )
    
    current_psi <- current_psi + increment
  }
  
  branch_df |>
    dplyr::distinct() |>
    dplyr::arrange(.data$psi)
}

# ======================================================================
# Full branch computation (left + right expansion)
# ======================================================================

.generate_branch <- function(branch_args, eval_psi_fun) {
  left <- .run_branch_side(
    increment     = -branch_args$increment,
    start         = branch_args$left_start,
    branch_cutoff = branch_args$branch_cutoff,
    init_guess    = branch_args$init_guess,
    eval_psi_fun  = eval_psi_fun,
    max_retries   = branch_args$max_retries
  )
  
  right <- .run_branch_side(
    increment     = branch_args$increment,
    start         = branch_args$right_start,
    branch_cutoff = branch_args$branch_cutoff,
    init_guess    = branch_args$init_guess,
    eval_psi_fun  = eval_psi_fun,
    max_retries   = branch_args$max_retries
  )
  
  dplyr::bind_rows(left, right) |>
    dplyr::arrange(.data$psi) |>
    dplyr::mutate(
      k = round((psi - branch_args$psi_MLE) / branch_args$increment)
    )
}
