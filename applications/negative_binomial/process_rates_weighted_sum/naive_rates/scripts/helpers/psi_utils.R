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
.compute_psi_MLE_SE <- function(theta_MLE, phi_MLE, weights, data) {
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
