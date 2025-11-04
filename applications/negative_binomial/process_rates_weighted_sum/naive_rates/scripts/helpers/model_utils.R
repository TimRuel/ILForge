# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/model_utils.R

suppressPackageStartupMessages({
  library(glmmTMB)
})

#' Fit a Negative Binomial model with process-specific rates and dispersions
#'
#' @param data A tibble or data frame containing:
#'   - Y: observed counts
#'   - t: exposure times
#'   - process: factor identifying each process
#'
#' @return A fitted glmmTMB model with process-level fixed effects and
#'         process-specific dispersion parameters.
#'
#' @details
#' The fitted model assumes:
#'   - Y_ij ~ NegBinomial(mu_ij, phi_j)
#'   - log(mu_ij) = log(t_ij) + log(theta_j)
#'
#' This corresponds to a regression with no intercept, process-level fixed
#' effects, and an offset for log-exposure.
#'
fit_model <- function(data) {
  stopifnot(all(c("Y", "t", "process") %in% names(data)))
  
  glmmTMB::glmmTMB(
    formula = Y ~ 0 + process + offset(log(t)),
    family = glmmTMB::nbinom2(),           # Variance = mu + phi * mu^2
    dispformula = ~ 0 + process,           # phi varies by process
    data = data
  )
}

#' Extract MLEs for process-specific rate parameters (theta_j)
#'
#' @param model A fitted glmmTMB object.
#' @return Named numeric vector of exp(fixed effects), one per process.
#'
get_theta_MLE <- function(model) {
  theta <- exp(glmmTMB::fixef(model)$cond)
  names(theta) <- levels(model$frame$process)
  theta
}

#' Extract MLEs for process-specific dispersion parameters (phi_j)
#'
#' @param model A fitted glmmTMB object.
#' @return Named numeric vector of exp(dispersion fixed effects), one per process.
#'
get_phi_MLE <- function(model) {
  phi <- exp(glmmTMB::fixef(model)$disp)
  names(phi) <- levels(model$frame$process)
  phi
}

#' Compute fitted mean counts per process (mu_hat)
#'
#' @param model A fitted glmmTMB object.
#' @param data The same data used for fitting, containing columns `t` and `process`.
#'
#' @return Numeric vector of fitted means `mu_hat = t * theta_hat` per observation.
#'         The vector is ordered the same as `data`.
#'
#' @details
#' This helper multiplies each observation’s exposure time by the corresponding
#' process-specific rate parameter estimated from the fitted model.
#'
get_mu_hat <- function(model, data) {
  stopifnot(all(c("t", "process") %in% names(data)))
  
  theta_hat <- get_theta_MLE(model)
  mu_hat <- data$t * theta_hat[data$process]
  mu_hat
}
