# applications/multinomial_logistic_regression/entropy/scripts/helpers/truth_utils.R

get_num_predictors <- function(X1_level_names, J) {
  
  ncol(model.matrix( ~ factor(X1_level_names)[1] * J - 1))
}

get_entropy_ranges <- function(X1_level_names, J, entropy_range_specs) {
  
  list2env(entropy_range_specs, environment())

  num_ranges <- length(X1_level_names)

  lower <- 0 + offset[1]

  upper <- log(J) - offset[2]

  entropy_ranges <- seq(lower, upper, length.out = num_ranges + 1) |>
    (\(x) mapply(c, x[-length(x)], x[-1], SIMPLIFY = FALSE))() |>
    map(\(x) {
      midpoint <- mean(x)
      desired_length <- (x[2] - x[1]) * padding
      range <- (midpoint + c(-1, 1) * desired_length / 2) |> 
        round(2)
      return(range)
    }
    )

  names(entropy_ranges) <- rev(X1_level_names)

  return(rev(entropy_ranges))
}

compute_probabilities <- function(intercept, slope, X2_vals) {

  eta <- matrix(intercept, nrow = length(X2_vals), ncol = length(intercept), byrow = TRUE) + outer(X2_vals, slope, `*`)
  exp_eta <- exp(eta)

  denom <- 1 + rowSums(exp_eta)
  probs <- cbind(1 / denom, exp_eta / denom)

  return(probs)
}

Beta_0_objective_fn <- function(params, X2_vals, Beta2) {

  Jm1 <- length(params) / 2
  intercept <- params[1:Jm1]
  slope <- params[(Jm1+1):(2*Jm1)] + Beta2

  prob_matrix <- compute_probabilities(intercept, slope, X2_vals)
  entropies <- apply(prob_matrix, 1, entropy)

  return(-diff(range(entropies)))
}

Beta_0_constraint_fn <- function(params, X2_vals, entropy_range, Beta2) {

  Jm1 <- length(params) / 2
  intercept <- params[1:Jm1]
  slope <- params[(Jm1+1):(2*Jm1)] + Beta2

  prob_matrix <- compute_probabilities(intercept, slope, X2_vals)
  entropies <- apply(prob_matrix, 1, entropy)

  return(c(entropy_range[1] - min(entropies), max(entropies) - entropy_range[2]))
}

optimize_Beta_0 <- function(J, entropy_range, X2_support, reps, Beta2) {

  X2_vals <- seq(X2_support[1], X2_support[2], length.out = reps)

  Jm1 <- J - 1
  init_params <- rnorm(2 * Jm1)

  result <- nloptr::auglag(
    x0 = init_params,
    fn = function(params) Beta_0_objective_fn(params, X2_vals, Beta2),
    lower = rep(-10, length(init_params)),
    upper = rep(10, length(init_params)),
    hin = function(params) Beta_0_constraint_fn(params, X2_vals, entropy_range, Beta2),
    deprecatedBehavior = FALSE)

  params <- result$par

  return(list(intercept = params[1:Jm1], slope = params[(Jm1+1):(2*Jm1)]))
}

get_Beta_0 <- function(X1_levels, J, p, reps) {

  X1_level_names <- names(X1_levels)

  X1_ref_level <- X1_levels |>
    map_lgl(\(x) x$ref_level) |>
    which() |>
    names()

  X1_main_effect_names <- paste0("X1", X1_level_names)
  
  X1X2_interaction_names <- paste0("X1", setdiff(X1_level_names, X1_ref_level)) |>
    paste0(":X2")

  Beta_0 <- matrix(NA,
                   nrow = p,
                   ncol = J - 1)
  
  rownames(Beta_0) <- c(X1_main_effect_names, "X2", X1X2_interaction_names)
  
  for (h in X1_level_names) {

    list2env(X1_levels[[h]]$X2, environment())

    if (h == X1_ref_level) {

      Beta2 <- 0

      params <- optimize_Beta_0(J, entropy_range, support, reps, Beta2)

      Beta_0[paste0("X1", h), ] <- params$intercept
      Beta_0["X2", ] <- params$slope
    } else {

      Beta2 <- Beta_0["X2", ]

      params <- optimize_Beta_0(J, entropy_range, support, reps, Beta2)

      Beta_0[grepl(h, rownames(Beta_0)), ] <- unname(rbind(params$intercept, params$slope))
    }
  }
  
  return(Beta_0)
}

get_experiment_parameters <- function(X1_levels, model_specs) {
  
  model_specs$p <- get_num_predictors(names(X1_levels), model_specs$J)
  
  list2env(model_specs, environment())
  
  entropy_ranges <- X1_levels |> 
    names() |> 
    get_entropy_ranges(J, entropy_range_specs)
  
  m <- map_int(X1_levels, \(x) x$m)
  
  model_specs$n <- sum(m)
  
  X1_levels <- X1_levels |> 
    imap(\(level, h) {
      level$X2$entropy_range <- entropy_ranges[[h]]
      level$m <- reps$pY_0
      level
    })
  
  Beta_0 <- get_Beta_0(X1_levels, J, p, reps$Beta_0)
  
  X1 <- get_X1(X1_levels)
  
  X2 <- get_X2(X1_levels)
  
  mm_formula <- substring(formula, 2)
  
  X_design <- get_X_design(mm_formula, X1, X2)
  
  pY_0 <- get_Y_probs(X_design, Beta_0)
  
  theta_0 <- pY_0 |> 
    select(-X2) |> 
    group_by(X1) |> 
    summarise(across(everything(), mean)) |> 
    data.frame()
  
  H_0 <- theta_0 |> 
    group_by(X1) |> 
    rowwise() |> 
    mutate(entropy = entropy(c_across(everything()))) |> 
    select(X1, entropy) |> 
    data.frame()
  
  experiment_parameters <- list(true_params = list(Beta_0 = Beta_0,
                                                   theta_0 = theta_0,
                                                   H_0 = H_0),
                                pY_0 = pY_0,
                                model_specs = model_specs)
  
  return(experiment_parameters)
}

