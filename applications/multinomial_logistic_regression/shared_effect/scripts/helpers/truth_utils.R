# applications/multinomial_logistic_regression/shared_effect/scripts/helpers/truth_utils.R

generate_true_parameters <- function(config) {
  
  set.seed(config$data_generation$seed)
  
  J <- config$model$response$num_classes
  num_effects <- J - 1
  
  # Generate intercepts
  intercept_dist <- config$model$intercepts$distribution
  intercept_args <- config$model$intercepts$parameters
  intercept_args$n <- num_effects
  intercepts <- list("(Intercepts)" = do.call(intercept_dist, intercept_args))
  
  # Read shared effect (psi_0)
  shared_effect <- config$model$predictors$shared$effect
  
  # Generate class-specific effects
  class_effects <- list()
  for (pred in config$model$predictors$class_specific) {
    dist_fn <- pred$effects$distribution
    dist_args <- pred$effects$parameters
    dist_args$n <- num_effects
    vec <- do.call(dist_fn, dist_args)
    class_effects[[pred$name]] <- vec
  }
  
  # Combine into Beta_0 matrix
  Beta_0 <- do.call(rbind, c(intercepts, class_effects)) |> 
    round(2)
  colnames(Beta_0) <- paste0("Y", seq_len(num_effects))
  
  list(psi_0 = shared_effect, Beta_0 = Beta_0)
}


