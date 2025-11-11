# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/omega_hat_utils.R

# -------------------------------------------------------------------------
# OMEGA_HAT UTILS (FULL THETA + PHI)
# -------------------------------------------------------------------------

#' @title Constraint closure for omega_hat (theta + phi)
#' @description Creates an equality constraint closure: sum(weights * theta) - psi_MLE.
#'   Only the theta components of omega are constrained; phi components are unconstrained.
#' @param weights Numeric vector of process weights (length J)
#' @param psi_MLE Numeric, target weighted sum
#' @param J Integer, number of processes
#' @return Function closure taking omega_hat vector (length 2*J, first J = theta, last J = phi)
.omega_hat_con_fn_closure <- function(weights, psi_MLE) {
  
  J <- length(weights)
  
  function(omega_hat) {
    theta_hat <- omega_hat[1:J]  # first J elements are theta
    sum(weights * theta_hat) - psi_MLE
  }
}

#' @title Optimize omega_hat given a constraint (theta + phi)
#' @description Wrapper for constrained optimization using auglag.
#'   omega_hat is now a vector of length 2*J: first J = theta, last J = phi.
#'   The objective is zero; the constraint is applied to theta only.
#' @param omega_hat_con_fn Function closure defining the equality constraint
#' @param init_guess Numeric vector of length 2*J, initial guess for omega_hat
#' @param localsolver Character, solver to use (default "SLSQP")
#' @param xtol_rel Numeric, relative tolerance for solver
#' @param maxeval Integer, maximum number of iterations
#' @return List with components:
#'   - theta: optimized theta vector
#'   - phi: optimized phi vector
.draw_omega_hat <- function(
    omega_hat_con_fn, 
    init_guess,
    ...
    ) {
  
  res <- nloptr::auglag(
    x0 = init_guess,
    fn = function(omega_hat) 0,      # Dummy objective, we only need the constraint
    heq = omega_hat_con_fn,
    lower = rep(1e-12, length(init_guess)),
    deprecatedBehavior = FALSE,
    ...
  )
  
  J <- length(init_guess) / 2
  list(
    theta = res$par[1:J],
    phi   = res$par[(J+1):(2*J)]
  )
}
