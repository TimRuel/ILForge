# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/psi_utils.R

suppressPackageStartupMessages({
  library(dplyr)
})

#' Compute weighted sum (psi)
#'
#' @param theta Numeric vector of process rates.
#' @param weights Numeric vector of process weights.
#' @return Scalar, the weighted sum of theta.
get_psi <- function(theta, weights) {
  sum(theta * weights)
}

#' Compute standard error of psi (weighted sum) from data
#'
#' This function calculates the variance and standard error of psi = sum(theta * weights)
#' using theta and phi MLEs and exposure times in the data.
#'
#' @param theta_MLE Named numeric vector of theta MLEs (per process).
#' @param phi_MLE Named numeric vector of phi MLEs (per process).
#' @param weights Named numeric vector of process weights.
#' @param data Data frame with columns `process` and `t` (exposure per observation).
#' @return Numeric scalar, the standard error of psi.
get_psi_MLE_SE <- function(theta_MLE, phi_MLE, weights, data) {
  se_terms <- data %>%
    group_by(process) %>%
    summarise(
      S1 = sum(t),
      S2 = sum(t^2),
      .groups = "drop"
    ) %>% as.data.frame()
  
  theta_vec <- theta_MLE[match(se_terms$process, names(theta_MLE))]
  phi_vec   <- phi_MLE[match(se_terms$process, names(phi_MLE))]
  
  var_theta <- theta_vec / se_terms$S1 + theta_vec^2 / phi_vec * se_terms$S2 / (se_terms$S1^2)
  sqrt(sum((weights^2) * var_theta))
}

#' Generate psi grid for profile likelihood
#'
#' Constructs a combined coarse and fine grid for psi, centered at the MLE
#' and covering a specified number of standard errors.
#'
#' @param theta_MLE Named numeric vector of theta MLEs.
#' @param phi_MLE Named numeric vector of phi MLEs.
#' @param weights Named numeric vector of process weights.
#' @param data Data frame with columns `process` and `t`.
#' @param num_std_errors Numeric, how many standard errors to span around MLE (default 3).
#' @param step_size Numeric, step size for coarse grid (default 0.1).
#' @param fine_step_size Numeric, step size for fine grid (default 0.01).
#' @param fine_window Numeric, half-width around MLE for fine grid (default 1).
#' @return Numeric vector of psi values (sorted, unique).
get_psi_grid <- function(theta_MLE, phi_MLE, weights, data,
                         num_std_errors = 3,
                         step_size = 0.1,
                         fine_step_size = 0.01,
                         fine_window = 1) {
  
  psi_MLE <- get_psi(theta_MLE, weights)
  psi_MLE_SE <- get_psi_MLE_SE(theta_MLE, phi_MLE, weights, data)
  endpoints <- psi_MLE + c(-1, 1) * num_std_errors * psi_MLE_SE
  
  coarse_grid <- seq(endpoints[1], endpoints[2], by = step_size)
  
  lower <- max(coarse_grid[coarse_grid <= (psi_MLE - fine_window)], min(coarse_grid))
  upper <- min(coarse_grid[coarse_grid >= (psi_MLE + fine_window)], max(coarse_grid))
  
  fine_grid <- seq(lower, upper, by = fine_step_size)
  
  sort(unique(c(coarse_grid, fine_grid, psi_MLE)))
}

get_search_interval <- function(model,
                                data,
                                weights, 
                                num_std_errors) {
  
  theta_MLE <- get_theta_MLE(model)
  
  phi_MLE <- get_phi_MLE(model)
  
  psi_MLE <- get_psi(theta_MLE, weights)
  
  psi_MLE_SE <- get_psi_MLE_SE(theta_MLE, phi_MLE, weights, data)
  
  psi_MLE + c(-1, 1) * num_std_errors * psi_MLE_SE
}
