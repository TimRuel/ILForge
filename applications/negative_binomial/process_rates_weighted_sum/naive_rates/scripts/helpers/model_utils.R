# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/model_utils.R

suppressPackageStartupMessages({
  library(glmmTMB)
})

fit_model <- function(data) {
  
  glmmTMB(
    Y ~ 0 + process + offset(log(t)),
    family = nbinom2(),        # variance = mu + phi*mu^2
    dispformula = ~ 0 + process, # φ varies by process
    data = data
  )
}

get_theta_MLE <- function(model) {
  
  model |> 
    fixef() |> 
    (\(x) x$cond)() |> 
    exp() |> 
    set_names(levels(model$frame$process))
}

get_phi_MLE <- function(model) {
  
  model |> 
    fixef() |> 
    (\(x) x$disp)() |> 
    exp() |> 
    set_names(levels(model$frame$process))
}

get_psi_MLE <- function(theta_MLE, weights) sum(theta_MLE * weights)

get_psi_MLE_SE <- function(theta_MLE, phi_MLE, weights, data) {
  
  t <- data |> 
    group_by(process) |> 
    summarise(S1 = sum(t),
              S2 = sum(t^2),
              .groups = "drop")
  
  theta_MLE <- theta_MLE[match(t$process, names(theta_MLE))]
  phi_MLE <- phi_MLE[match(t$process, names(phi_MLE))]
  
  var_theta <- theta_MLE / t$S1 + theta_MLE^2 / phi_MLE * t$S2 / (t$S1)^2
  
  sqrt(sum((weights^2) * var_theta))
}

log_likelihood <- function(theta, phi, Y, t, n_per_process) {
  
  dnbinom(x = Y, 
          size = rep(phi, times = n_per_process), 
          mu = t * rep(theta, times = n_per_process),
          log = TRUE) |>
    sum()
}

likelihood <- function(theta, phi, Y, t, n_per_process) exp(log_likelihood(theta, phi, Y, t, n_per_process))
