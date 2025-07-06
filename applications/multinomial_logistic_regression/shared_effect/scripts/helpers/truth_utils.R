# applications/multinomial_logistic_regression/shared_effect/scripts/helpers/truth_utils.R

get_Beta_0_from_config <- function(config) {
  
  model <- config$model
  num_effects <- model$response$num_effects
  
  # Initialize list to hold rows of Beta
  beta_rows <- list()
  row_names <- c()
  
  # Intercepts
  intercepts <- model$intercepts
  stopifnot(length(intercepts) == num_effects)
  beta_rows[["Intercept"]] <- intercepts
  
  # Shared predictors
  if (!is.null(model$predictors$shared)) {
    for (pred in model$predictors$shared) {
      name <- pred$name
      effect <- pred$effect
      beta_rows[[name]] <- rep(effect, num_effects)
    }
  }
  
  # Class-specific predictors
  if (!is.null(model$predictors$class_specific)) {
    for (pred in model$predictors$class_specific) {
      name <- pred$name
      effects <- pred$effects
      stopifnot(length(effects) == num_effects)
      beta_rows[[name]] <- effects
    }
  }
  
  # Combine into matrix
  Beta_0 <- do.call(rbind, beta_rows)
  colnames(Beta_0) <- paste0("Class", seq_len(num_effects))
  rownames(Beta_0) <- names(beta_rows)
  
  return(Beta_0)
}
