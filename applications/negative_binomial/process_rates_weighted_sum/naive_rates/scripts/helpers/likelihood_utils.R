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
# HYBRID EXPECTED GAMMA (EXACT + MONTE CARLO)
# -------------------------------------------------------------------------

#' Hybrid expected gamma function (exact + Monte Carlo)
#'
#' Computes E[F(Y_ij + phi_j) | omega_j] using exact summation for moderate overdispersion
#' and Monte Carlo sampling for extreme overdispersion.
#'
#' @param omega List with components `theta` and `phi`, evaluation point.
#' @param theta Numeric vector of true theta_j parameters.
#' @param phi Numeric vector of true phi_j parameters.
#' @param t Numeric vector of exposure times for each observation.
#' @param n_per_process Integer vector giving number of observations per process.
#' @param FUN Function to apply: typically `lgamma` or `digamma`.
#' @param ratio_threshold Numeric, mu/phi ratio threshold for Monte Carlo (default 100).
#' @param n_mc Integer, number of Monte Carlo samples (default 1e5).
#' @param p_cutoff Numeric, tail cutoff for exact summation (default 1e-12).
#' @param max_y_cap Numeric, maximum y to sum to for exact computation (default 1e6).
#' @return Numeric vector of expected values per observation.
#' @export
E_gamma_hybrid <- function(omega, theta, phi, t, n_per_process,
                           FUN = lgamma, ratio_threshold = 100, n_mc = 1e5,
                           p_cutoff = 1e-12, max_y_cap = 1e6) {
  
  J <- length(theta)
  n <- sum(n_per_process)
  obs_to_process <- rep(seq_len(J), times = n_per_process)
  
  theta_vec <- rep(theta, times = n_per_process)
  phi_vec   <- rep(phi, times = n_per_process)
  mu_vec    <- t * theta_vec
  
  omega_theta_vec <- rep(omega$theta, times = n_per_process)
  omega_phi_vec   <- rep(omega$phi, times = n_per_process)
  omega_mu_vec    <- t * omega_theta_vec
  
  ratio <- omega_mu_vec / omega_phi_vec
  E_vals <- numeric(n)
  
  # Exact computation
  idx_exact <- which(ratio <= ratio_threshold)
  if (length(idx_exact) > 0) {
    E_vals[idx_exact] <- E_gamma_vec(
      omega = list(theta = omega_theta_vec[idx_exact], phi = omega_phi_vec[idx_exact]),
      theta = theta_vec[idx_exact],
      phi   = phi_vec[idx_exact],
      t     = t[idx_exact],
      n_per_process = rep(1, length(idx_exact)),
      FUN   = FUN,
      p_cutoff = p_cutoff,
      max_y_cap = max_y_cap
    )
  }
  
  # Monte Carlo for high ratios
  idx_mc <- which(ratio > ratio_threshold)
  if (length(idx_mc) > 0) {
    y_samp <- matrix(
      rnbinom(n_mc * length(idx_mc),
              mu = omega_mu_vec[idx_mc],
              size = omega_phi_vec[idx_mc]),
      nrow = n_mc, ncol = length(idx_mc)
    )
    phi_mat <- matrix(rep(phi_vec[idx_mc], each = n_mc), nrow = n_mc)
    E_vals[idx_mc] <- colMeans(FUN(y_samp + phi_mat))
  }
  
  E_vals
}

#' Expected log-gamma using hybrid method
#'
#' @inheritParams E_gamma_hybrid
#' @export
E_log_gamma_hybrid <- function(omega, theta, phi, t, n_per_process, ...) {
  E_gamma_hybrid(omega, theta, phi, t, n_per_process, FUN = lgamma, ...)
}

#' Expected digamma using hybrid method
#'
#' @inheritParams E_gamma_hybrid
#' @export
E_digamma_hybrid <- function(omega, theta, phi, t, n_per_process, ...) {
  E_gamma_hybrid(omega, theta, phi, t, n_per_process, FUN = digamma, ...)
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
#' @param ratio_threshold Numeric, threshold for hybrid expectation.
#' @param n_mc Integer, number of Monte Carlo samples.
#' @param p_cutoff Numeric, tail cutoff for exact summation.
#' @param max_y_cap Numeric, maximum y to sum to in exact expectation.
#' @return Function closure: takes theta_phi vector and returns expected log-likelihood scalar.
#' @export
expected_log_likelihood_closure <- function(omega, t, n_per_process,
                                            ratio_threshold = 100, n_mc = 1e5,
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
    
    E_log_gamma_vals <- E_log_gamma_hybrid(
      omega, theta, phi, t, n_per_process,
      ratio_threshold = ratio_threshold,
      n_mc = n_mc,
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
expected_gradient_theta <- function(phi_vec, theta_vec, mu_vec, omega_mu_vec, obs_to_process) {
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
#' @param ratio_threshold Numeric, hybrid expectation threshold
#' @param n_mc Integer, Monte Carlo samples
#' @param p_cutoff Numeric, tail cutoff
#' @param max_y_cap Numeric, maximum y for exact expectation
#' @return Function closure taking theta_phi vector, returning gradient vector
#' @export
expected_gradient_closure <- function(omega, t, n_per_process,
                                      ratio_threshold = 100, n_mc = 1e5,
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
    
    E_digamma_vals <- E_digamma_hybrid(
      omega, theta, phi, t, n_per_process,
      ratio_threshold = ratio_threshold,
      n_mc = n_mc,
      p_cutoff = p_cutoff,
      max_y_cap = max_y_cap
    )
    
    grad_theta <- expected_gradient_theta(phi_vec, theta_vec, mu_vec, omega_mu_vec, obs_to_process)
    grad_phi   <- expected_gradient_phi(phi_vec, mu_vec, omega_mu_vec, E_digamma_vals, obs_to_process)
    
    c(grad_theta, grad_phi)
  }
}
