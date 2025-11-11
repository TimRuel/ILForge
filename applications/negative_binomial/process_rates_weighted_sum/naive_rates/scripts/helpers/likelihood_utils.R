# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/likelihood_utils.R

# -------------------------------------------------------------------------
# LOG-LIKELIHOOD
# -------------------------------------------------------------------------

#' Log-likelihood for Negative Binomial
#'
#' Computes the log-likelihood of a negative binomial model with per-process
#' rates `theta` and dispersion `phi`.
#'
#' @param theta Numeric vector of process rates.
#' @param phi Numeric vector of dispersion parameters.
#' @param Y Numeric vector of observed counts.
#' @param t Numeric vector of exposure times.
#' @param n_per_process Integer vector with number of observations per process.
#' @return Numeric scalar: log-likelihood.
#' @export
log_likelihood <- function(theta, phi, Y, t, n_per_process) {
  dnbinom(
    x = Y,
    size = rep(phi, times = n_per_process),
    mu = t * rep(theta, times = n_per_process),
    log = TRUE
  ) |> sum()
}

#' Likelihood for Negative Binomial
#'
#' Exponentiates the log-likelihood.
#'
#' @inheritParams log_likelihood
#' @return Numeric scalar: likelihood.
#' @export
likelihood <- function(theta, phi, Y, t, n_per_process) {
  exp(log_likelihood(theta, phi, Y, t, n_per_process))
}

# -------------------------------------------------------------------------
# Exact expectation of gamma/digamma-related functions under NB
# -------------------------------------------------------------------------

#' Exact expected value of F(Y + φ) under NB(θ, φ), per observation
#'
#' Computes the exact expectation
#' \deqn{E_\omega \left[ F(Y_{ij} + \phi_j) \right]}
#' for each observation, where
#' \eqn{Y_{ij} \sim \text{NB}(\mu_{ij} = t_{ij} \theta_{\omega j}, \phi_{\omega j})}.
#'
#' The expectation is evaluated exactly by summing over the Negative Binomial
#' pmf truncated at a negligible tail mass (controlled by `p_cutoff`). This
#' function returns one expected value per observation and is typically used
#' inside integrated-likelihood components.
#'
#' @param omega List with elements `theta` and `phi`, representing the parameters
#'   defining the expectation distribution (the “with respect to” NB parameters).
#'   Each will be recycled to match total observations implied by `n_per_process`.
#' @param phi Numeric vector of φ values used inside `F(Y + φ)`; length equals total observations.
#' @param t Numeric vector of exposure times (one per observation).
#' @param n_per_process Integer vector giving the number of observations per process.
#' @param FUN Function to apply to \eqn{Y + φ}, typically `lgamma` or `digamma`.
#' @param p_cutoff Numeric; truncate NB tail when cumulative mass exceeds `1 - p_cutoff`.
#' @param max_y_cap Integer; maximum `y` value to include for exact computation (default 1e6).
#'
#' @return Numeric vector of expected values, one per observation.
#'
#' @examples
#' omega <- list(theta = c(2, 3), phi = c(5, 8))
#' phi   <- rep(c(5, 8), each = 10)
#' t     <- runif(20, 1, 3)
#' n_per_process <- c(10, 10)
#' E_gamma_vec(
#'   omega = omega,
#'   phi = phi,
#'   t = t,
#'   n_per_process = n_per_process,
#'   FUN = lgamma
#' )
#'
#' @export
E_gamma <- function(omega,
                    phi,
                    t,
                    n_per_process,
                    FUN = lgamma,
                    p_cutoff = 1e-12,
                    max_y_cap = 1e6) {
  
  n <- sum(n_per_process)
  
  # --- Expand process-level params to observation-level ---
  omega_mu_vec <- t * rep(omega$theta, times = n_per_process)
  omega_phi_vec <- rep(omega$phi,   times = n_per_process)
  
  # --- Initialize output ---
  E_gamma_vals <- numeric(n)
  
  # --- Loop over observations ---
  for (i in seq_len(n)) {
    # Build NB support
    y_max <- min(max_y_cap, ceiling(qnbinom(1 - p_cutoff, size = omega_phi_vec[i], mu = omega_mu_vec[i])))
    y_seq <- 0:y_max
    p_y   <- dnbinom(y_seq, size = omega_phi_vec[i], mu = omega_mu_vec[i])
    p_y   <- p_y / sum(p_y)  # normalize for safety
    
    # Compute expected value for this observation
    E_gamma_vals[i] <- sum(FUN(y_seq + phi[i]) * p_y)
  }
  
  E_gamma_vals
}

#' Expected log-gamma 
#'
#' @inheritParams E_gamma
#' @export
E_log_gamma <- function(omega, phi, t, n_per_process, ...) {
  E_gamma(omega, phi, t, n_per_process, FUN = lgamma, ...)
}

#' Expected digamma
#'
#' @inheritParams E_gamma
#' @export
E_digamma <- function(omega, phi, t, n_per_process, ...) {
  E_gamma(omega, phi, t, n_per_process, FUN = digamma, ...)
}

# -------------------------------------------------------------------------
# EXPECTED LOG-LIKELIHOOD CLOSURE
# -------------------------------------------------------------------------

#' Create expected log-likelihood closure for (theta, phi)
#'
#' Returns a closure computing the expected log-likelihood of the model
#' given parameters theta, phi. Expectation is w.r.t. omega.
#'
#' @param omega List with components `theta` and `phi`, evaluation point.
#' @param t Numeric vector of exposure times.
#' @param n_per_process Integer vector of number of observations per process.
#' @param p_cutoff Numeric, tail cutoff for exact summation.
#' @param max_y_cap Numeric, maximum y to sum to in exact expectation.
#' @return Function closure: takes theta_phi vector and returns expected log-likelihood scalar.
#' @export
expected_log_likelihood_closure <- function(omega, t, n_per_process,
                                            p_cutoff = 1e-12, max_y_cap = 1e6) {
  
  J <- length(omega$theta)
  omega_mu_vec <- t * rep(omega$theta, times = n_per_process)
  
  function(theta_phi) {
    theta <- theta_phi[1:J]
    phi   <- theta_phi[(J+1):(2*J)]
    
    theta_vec <- rep(theta, times = n_per_process)
    phi_vec   <- rep(phi, times = n_per_process)
    mu_vec    <- t * theta_vec
    
    lgamma_phi  <- lgamma(phi_vec)
    phi_log_phi <- phi_vec * log(phi_vec)
    
    E_log_gamma_vals <- E_log_gamma(
      omega, phi_vec, t, n_per_process,
      p_cutoff = p_cutoff,
      max_y_cap = max_y_cap
    )
    
    theta_term <- omega_mu_vec * log(mu_vec) - (omega_mu_vec + phi_vec) * log(mu_vec + phi_vec)
    
    sum(E_log_gamma_vals - lgamma_phi + phi_log_phi + theta_term)
  }
}

# -------------------------------------------------------------------------
# PARTIAL GRADIENTS
# -------------------------------------------------------------------------

#' Gradient w.r.t. theta
#' @param phi_vec Numeric vector of phi values repeated to observation level
#' @param theta_vec Numeric vector of theta values repeated to observation level
#' @param mu_vec Numeric vector of expected counts (t * theta_vec)
#' @param omega_mu_vec Numeric vector of omega_mu (t * omega$theta)
#' @param obs_to_process Integer vector mapping each observation to its process
#' @return Numeric vector length J: gradient w.r.t. theta
#' @export
expected_gradient_theta <- function(phi_vec, theta_vec, mu_vec, omega_mu_vec, t, obs_to_process) {
  grad_obs <- omega_mu_vec / theta_vec - t * (omega_mu_vec + phi_vec) / (mu_vec + phi_vec)
  as.numeric(rowsum(grad_obs, group = obs_to_process))
}

#' Gradient w.r.t. phi
#' @param phi_vec Numeric vector of phi values repeated to observation level
#' @param mu_vec Numeric vector of expected counts
#' @param omega_mu_vec Numeric vector of omega_mu
#' @param E_digamma_vals Expected digamma values per observation
#' @param obs_to_process Observation-to-process mapping
#' @return Numeric vector length J: gradient w.r.t. phi
#' @export
expected_gradient_phi <- function(phi_vec, mu_vec, omega_mu_vec, E_digamma_vals, obs_to_process) {
  grad_obs <- E_digamma_vals - digamma(phi_vec) + log(phi_vec) + 1 -
    log(phi_vec + mu_vec) - (omega_mu_vec + phi_vec) / (mu_vec + phi_vec)
  as.numeric(rowsum(grad_obs, group = obs_to_process))
}

# -------------------------------------------------------------------------
# FULL GRADIENT CLOSURE
# -------------------------------------------------------------------------

#' Create closure for expected gradient of the log-likelihood
#'
#' Returns a function closure computing the gradient w.r.t. theta and phi
#' for a given evaluation point omega.
#'
#' @param omega List with components `theta` and `phi`
#' @param t Numeric vector of exposure times
#' @param n_per_process Integer vector of observations per process
#' @param p_cutoff Numeric, tail cutoff
#' @param max_y_cap Numeric, maximum y for exact expectation
#' @return Function closure taking theta_phi vector, returning gradient vector
#' @export
expected_gradient_closure <- function(omega, t, n_per_process,
                                      p_cutoff = 1e-12, max_y_cap = 1e6) {
  
  J <- length(omega$theta)
  omega_mu_vec <- t * rep(omega$theta, times = n_per_process)
  obs_to_process <- rep(seq_len(J), times = n_per_process)
  
  function(theta_phi) {
    theta <- theta_phi[1:J]
    phi   <- theta_phi[(J+1):(2*J)]
    
    theta_vec <- rep(theta, times = n_per_process)
    phi_vec   <- rep(phi, times = n_per_process)
    mu_vec    <- t * theta_vec
    
    E_digamma_vals <- E_digamma(
      omega, phi_vec, t, n_per_process,
      p_cutoff = p_cutoff,
      max_y_cap = max_y_cap
    )
    
    grad_theta <- expected_gradient_theta(phi_vec, theta_vec, mu_vec, omega_mu_vec, t, obs_to_process)
    grad_phi   <- expected_gradient_phi(phi_vec, mu_vec, omega_mu_vec, E_digamma_vals, obs_to_process)
    
    c(grad_theta, grad_phi)
  }
}
