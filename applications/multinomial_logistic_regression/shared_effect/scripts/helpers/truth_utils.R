# applications/multinomial_logistic_regression/shared_effect/scripts/helpers/truth_utils.R

generate_true_parameters <- function(config) {
  set.seed(config$data_generation$seed)
  
  # Extract formula terms
  form <- as.formula(config$model$formula)
  terms_obj <- terms(form)
  all_terms <- attr(terms_obj, "term.labels")
  
  J <- config$model$response$num_classes
  num_effects <- J - 1
  
  # Prepare intercepts
  intercept_dist <- config$model$intercepts$distribution
  intercept_args <- config$model$intercepts$parameters
  intercept_args$n <- num_effects
  intercepts <- list("(Intercept)" = do.call(intercept_dist, intercept_args))
  
  # Determine shared vs. class-specific terms in formula
  shared_name <- config$model$predictors$shared$name
  psi_0 <- NULL
  if (shared_name %in% all_terms) {
    psi_0 <- config$model$predictors$shared$effect
  }
  
  # Class-specific effects
  class_effects <- list()
  cs_config <- config$model$predictors$class_specific
  
  # Create a lookup for class-specific definitions
  cs_lookup <- setNames(cs_config, vapply(cs_config, `[[`, "", "name"))
  
  # Loop over all formula terms and generate appropriate coefficients
  for (term in all_terms) {
    if (term == shared_name) next  # Already handled as psi_0
    
    # Get distribution settings from the matching class-specific predictor
    dist_fn <- cs_lookup[[term]]$effects$distribution
    dist_args <- cs_lookup[[term]]$effects$parameters
    dist_args$n <- num_effects
    
    # Generate vector of coefficients
    class_effects[[term]] <- do.call(dist_fn, dist_args)
  }
  
  # Combine intercept and class-specific effects into matrix
  Beta_0 <- do.call(rbind, c(intercepts, class_effects)) |> round(2)
  colnames(Beta_0) <- paste0("Y", seq_len(num_effects))
  
  list(psi_0 = psi_0, Beta_0 = Beta_0)
}




