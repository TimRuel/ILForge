# applications/multinomial_logistic_regression/shared_effect/scripts/helpers/truth_utils.R

generate_true_parameters <- function(config) {
  set.seed(config$data_generation$seed)
  
  form <- as.formula(config$model$formula)
  terms_obj <- terms(form)
  all_terms <- attr(terms_obj, "term.labels")
  
  J <- config$model$response$num_classes
  num_effects <- J - 1
  
  # Intercept
  intercept_dist <- config$model$intercepts$distribution
  intercept_args <- config$model$intercepts$parameters
  intercept_args$n <- num_effects
  intercepts <- list("(Intercept)" = do.call(intercept_dist, intercept_args))
  
  # Shared predictor (e.g., Z)
  shared_name <- config$model$predictors$shared$name
  psi_0 <- NULL
  if (shared_name %in% all_terms) {
    psi_0 <- config$model$predictors$shared$effect
  }
  
  # Class-specific predictors
  cs_effects <- list()
  cs_config <- config$model$predictors$class_specific
  cs_lookup <- setNames(cs_config, vapply(cs_config, `[[`, "", "name"))
  
  for (term in all_terms) {
    if (term == shared_name || grepl(":", term)) next
    if (!term %in% names(cs_lookup)) {
      warning(sprintf("Term '%s' not found in class_specific config", term))
      next
    }
    dist_fn <- cs_lookup[[term]]$effects$distribution
    dist_args <- cs_lookup[[term]]$effects$parameters
    dist_args$n <- num_effects
    cs_effects[[term]] <- do.call(dist_fn, dist_args)
  }
  
  # Interactions
  interaction_effects <- list()
  interaction_config <- config$model$predictors$interactions %||% list()
  interaction_lookup <- setNames(interaction_config, vapply(interaction_config, `[[`, "", "name"))
  
  for (term in all_terms) {
    if (!grepl(":", term)) next
    if (!term %in% names(interaction_lookup)) {
      warning(sprintf("Interaction '%s' not found in config$model$interactions", term))
      next
    }
    dist_fn <- interaction_lookup[[term]]$effects$distribution
    dist_args <- interaction_lookup[[term]]$effects$parameters
    dist_args$n <- num_effects
    interaction_effects[[term]] <- do.call(dist_fn, dist_args)
  }
  
  # Combine all effects
  Beta_0 <- do.call(rbind, c(intercepts, cs_effects, interaction_effects)) |> round(2)
  colnames(Beta_0) <- paste0("Y", seq_len(num_effects))
  
  list(psi_0 = psi_0, Beta_0 = Beta_0)
}





