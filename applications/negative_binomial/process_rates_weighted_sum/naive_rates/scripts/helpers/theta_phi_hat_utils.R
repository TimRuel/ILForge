# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/theta_phi_hat_utils.R

# -------------------------------------------------------------------------
# OBJECTIVE, GRADIENT, AND CONSTRAINT CLOSURES
# -------------------------------------------------------------------------

#' Create closures for constrained estimation of (theta, phi)
#'
#' Returns closures for the negative expected log-likelihood objective,
#' its analytic gradient, and the weighted-sum constraint and Jacobian.
#'
#' The expected log-likelihood and gradient are computed using precomputation
#' based on \code{omega_hat}. The constraint enforces
#' \eqn{\sum_j w_j \theta_j = \psi}, while allowing \eqn{\phi_j} to vary freely.
#'
#' @param omega_hat List with components `theta` and `phi`, evaluation point for omega.
#' @param t Numeric vector of exposure times for each observation.
#' @param n_per_process Integer vector giving number of observations per process.
#' @param weights Numeric vector of process weights.
#' @param psi Numeric target value for weighted sum of theta.
#' @param ratio_threshold Numeric, threshold for hybrid expectation (default 100).
#' @param n_mc Integer, number of Monte Carlo draws for hybrid expectation (default 1e5).
#' @param p_cutoff Numeric, probability cutoff for truncating exact expectation (default 1e-12).
#' @param max_y_cap Numeric, maximum y value for truncating expectation (default 1e6).
#'
#' @return List with closures:
#'   \describe{
#'     \item{\code{objective}}{Negative expected log-likelihood function of \code{theta_phi}.}
#'     \item{\code{gradient}}{Analytic gradient of the objective.}
#'     \item{\code{constraint}}{Equality constraint enforcing \eqn{\sum_j w_j \theta_j = \psi}.}
#'     \item{\code{jacobian}}{Jacobian (row vector) of the constraint.}
#'   }
#' @seealso \code{\link{expected_log_likelihood_closure}},
#'   \code{\link{expected_gradient_closure}}
#' @export
make_theta_phi_hat_closures <- function(omega_hat, t, n_per_process,
                                        weights, psi,
                                        ratio_threshold = 100, n_mc = 1e5,
                                        p_cutoff = 1e-12, max_y_cap = 1e6) {
  
  J <- length(weights)
  
  # Precompute expected log-likelihood and gradient closures
  ell_closure <- expected_log_likelihood_closure(
    omega = omega_hat,
    t = t,
    n_per_process = n_per_process,
    ratio_threshold = ratio_threshold,
    n_mc = n_mc,
    p_cutoff = p_cutoff,
    max_y_cap = max_y_cap
  )
  
  grad_closure <- expected_gradient_closure(
    omega = omega_hat,
    t = t,
    n_per_process = n_per_process,
    ratio_threshold = ratio_threshold,
    n_mc = n_mc,
    p_cutoff = p_cutoff,
    max_y_cap = max_y_cap
  )
  
  # Objective and gradient (negated for minimization)
  objective <- function(theta_phi) -ell_closure(theta_phi)
  gradient  <- function(theta_phi) -grad_closure(theta_phi)
  
  # Constraint closure: enforce sum_j w_j * theta_j = psi
  constraint <- function(theta_phi) {
    theta <- theta_phi[1:J]
    sum(weights * theta) - psi
  }
  
  # Jacobian closure: derivatives w.r.t. theta only
  jacobian <- function(theta_phi) {
    cbind(matrix(weights, nrow = 1), matrix(0, nrow = 1, ncol = J))
  }
  
  list(
    objective  = objective,
    gradient   = gradient,
    constraint = constraint,
    jacobian   = jacobian
  )
}

# -------------------------------------------------------------------------
# CONSTRAINED OPTIMIZATION WRAPPER
# -------------------------------------------------------------------------

#' Estimate theta_phi_hat = (theta_hat, phi_hat) subject to weighted-sum constraint
#'
#' Solves the constrained optimization problem
#' \deqn{\sum_j w_j \theta_j = \psi}
#' jointly estimating \eqn{\theta} and \eqn{\phi}.
#'
#' @param omega_hat List with reference values (`theta`, `phi`) for omega precomputation.
#' @param t Numeric vector of exposure times.
#' @param n_per_process Integer vector of number of observations per process.
#' @param init_guess Numeric vector of length 2*J (theta + phi), initial guess for optimization.
#' @param weights Numeric vector of process weights.
#' @param psi Numeric target for weighted sum constraint.
#' @param ratio_threshold,n_mc,p_cutoff,max_y_cap Control parameters for hybrid expectation.
#' @param localsolver,xtol_rel,maxeval Optimization control arguments for auglag.
#' @param return_full Logical; if TRUE, also returns the full `nloptr` result.
#'
#' @return List with:
#'   \describe{
#'     \item{\code{theta_hat}}{Estimated theta values.}
#'     \item{\code{phi_hat}}{Estimated phi values.}
#'     \item{\code{result}}{(Optional) Full `nloptr` object if \code{return_full = TRUE}.}
#'   }
#' @seealso \code{\link{make_theta_phi_hat_closures}}
#' @export
get_theta_phi_hat <- function(omega_hat, t, n_per_process,
                              init_guess, weights, psi,
                              ratio_threshold = 100, n_mc = 1e5,
                              p_cutoff = 1e-12, max_y_cap = 1e6,
                              localsolver = "SLSQP",
                              xtol_rel = 1e-8,
                              maxeval = 1000,
                              return_full = FALSE) {
  
  closures <- make_theta_phi_hat_closures(
    omega_hat = omega_hat,
    t = t,
    n_per_process = n_per_process,
    weights = weights,
    psi = psi,
    ratio_threshold = ratio_threshold,
    n_mc = n_mc,
    p_cutoff = p_cutoff,
    max_y_cap = max_y_cap
  )
  
  res <- nloptr::auglag(
    x0 = init_guess,
    fn = closures$objective,
    gr = closures$gradient,
    heq = closures$constraint,
    heqjac = closures$jacobian,
    lower = rep(1e-12, length(init_guess)),
    localsolver = localsolver,
    control = list(xtol_rel = xtol_rel, maxeval = maxeval),
    deprecatedBehavior = FALSE
  )
  
  J <- length(weights)
  out <- list(
    theta_hat = res$par[1:J],
    phi_hat   = res$par[(J + 1):(2 * J)]
  )
  if (return_full) out$result <- res
  out
}
