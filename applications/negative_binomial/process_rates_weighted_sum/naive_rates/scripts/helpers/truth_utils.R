# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/truth_utils.R

#' Generate process labels automatically
#'
#' @param n_processes Number of processes to label.
#' @param style Label style: "auto_upper", "auto_lower", "auto_num", or a custom character vector.
#' @return Character vector of labels.
generate_process_labels <- function(n_processes, style = "auto_upper") {
  if (is.character(style) && length(style) == 1) {
    return(
      switch(style,
             auto_upper = {
               labels <- character(n_processes)
               for (i in seq_len(n_processes)) {
                 num <- i - 1
                 name <- ""
                 repeat {
                   name <- paste0(LETTERS[(num %% 26) + 1], name)
                   num <- num %/% 26 - 1
                   if (num < 0) break
                 }
                 labels[i] <- name
               }
               labels
             },
             auto_lower = tolower(generate_process_labels(n_processes, "auto_upper")),
             auto_num   = as.character(seq_len(n_processes)),
             stop(sprintf("Unknown label style: '%s'", style))
      )
    )
  }
  
  if (length(style) == n_processes) return(as.character(style))
  stop("Labels must be a known style keyword or a vector of length n_processes.")
}

#' Generate a parameter vector (e.g., theta, phi)
#'
#' @param cfg Full configuration list.
#' @param name Name of the parameter entry ("theta", "phi", etc.).
#' @return Numeric vector of parameter values.
generate_param_values <- function(cfg, name) {
  n_processes <- cfg$processes$n_processes
  param_cfg <- cfg[[name]]
  
  # Case 1: explicit vector of values provided
  if (!is.null(param_cfg$values)) {
    validate_param_vector(param_cfg$values, n_processes, name)
    return(param_cfg$values)
  }
  
  # Case 2: distribution-based generation
  if (is.null(param_cfg$distribution) || is.null(param_cfg$args)) {
    stop(sprintf(
      "%s must include either 'values' or ('distribution' and 'args') in the config.",
      name
    ))
  }
  
  dist_fun  <- match.fun(param_cfg$distribution)
  dist_args <- param_cfg$args
  do.call(dist_fun, c(list(n = n_processes), as.list(dist_args)))
}

#' Validate a parameter vector
#'
#' @param param Numeric vector of parameter values.
#' @param expected_length Expected length.
#' @param name Parameter name (for messaging).
#' @return Invisibly TRUE if valid; otherwise throws error.
validate_param_vector <- function(param, expected_length, name) {
  if (length(param) != expected_length) {
    stop(sprintf(
      "%s vector length (%d) does not match n_processes (%d).",
      name, length(param), expected_length
    ))
  }
  if (any(param < 0)) stop(sprintf("%s values must be non-negative.", name))
  invisible(TRUE)
}

#' Generate process weights
#'
#' @param cfg Full configuration list.
#' @return Named numeric vector of weights.
generate_process_weights <- function(cfg) {
  weight_cfg <- cfg$weights
  n_processes <- cfg$processes$n_processes
  
  if (is.null(weight_cfg$distribution) || is.null(weight_cfg$args)) {
    stop("Weight distribution and arguments must be specified in the config.")
  }
  
  dist_fun  <- match.fun(weight_cfg$distribution)
  dist_args <- weight_cfg$args
  raw_weights <- do.call(dist_fun, c(list(n = n_processes), as.list(dist_args)))
  
  weights <- raw_weights
  if (!is.null(weight_cfg$normalize_mean_to)) {
    weights <- raw_weights / mean(raw_weights) * weight_cfg$normalize_mean_to
  } else if (!is.null(weight_cfg$normalize_sum_to)) {
    weights <- raw_weights / sum(raw_weights) * weight_cfg$normalize_sum_to
  }
  
  labels <- generate_process_labels(n_processes, cfg$processes$labels)
  names(weights) <- labels
  weights
}

#' Expand n_per_process to match number of processes
#'
#' @param n_per_process Single integer or vector.
#' @param n_processes Total number of processes.
#' @return Integer vector of length n_processes.
expand_n_per_process <- function(n_per_process, n_processes) {
  if (length(n_per_process) == 1) {
    return(rep(n_per_process, n_processes))
  }
  if (length(n_per_process) < n_processes) {
    return(rep_len(n_per_process, n_processes))
  }
  if (length(n_per_process) == n_processes) {
    return(n_per_process)
  }
  stop("n_per_process must be length 1, a divisor of n_processes, or exactly n_processes.")
}

#' Generate true parameter values for all processes
#'
#' @param cfg Full configuration list for an experiment.
#' @return List containing theta_0, phi_0, weights, and n_per_process.
generate_true_parameters <- function(cfg) {
  n_processes <- cfg$processes$n_processes
  labels <- generate_process_labels(n_processes, cfg$processes$labels)
  
  theta_vals <- generate_param_values(cfg, "theta")
  phi_vals   <- generate_param_values(cfg, "phi")
  
  validate_param_vector(theta_vals, n_processes, "theta")
  validate_param_vector(phi_vals,   n_processes, "phi")
  
  weights <- generate_process_weights(cfg)
  n_per_process <- expand_n_per_process(cfg$processes$n_per_process, n_processes)
  
  names(theta_vals)     <- labels
  names(phi_vals)       <- labels
  names(weights)        <- labels
  names(n_per_process)  <- labels
  
  list(
    theta_0       = theta_vals,
    phi_0         = phi_vals,
    weights       = weights,
    n_per_process = n_per_process
  )
}
