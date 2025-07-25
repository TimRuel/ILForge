# applications/multinomial_logistic_regression/shared_effect/scripts/helpers/data_utils.R

get_eta <- function(Z, X, psi, Beta) psi * Z + X %*% Beta

softmax <- function(x) exp(x) / sum(exp(x))

softmax_adj <- function(x) exp(x) / (1 + sum(exp(x)))

generate_data <- function(config, theta_0) {
  
  n <- config$data_generation$n_obs
  model <- config$model
  num_classes <- model$response$num_classes
  ref_class <- model$response$reference_class
  num_effects <- num_classes - 1
  min_obs <- config$data_generation$min_obs
  max_tries <- config$data_generation$max_tries %||% 100
  
  formula_obj <- as.formula(model$formula)
  
  # Remove the shared predictor Z from the RHS
  rhs_terms <- attr(terms(formula_obj), "term.labels")
  rhs_terms_no_Z <- setdiff(rhs_terms, model$predictors$shared$name)
  rhs_formula <- reformulate(rhs_terms_no_Z)
  
  for (i in 1:max_tries) {
    
    # Generate shared predictor Z
    shared <- model$predictors$shared
    Z_args <- shared$observations$parameters
    Z_args$n <- n
    Z <- do.call(shared$observations$distribution, Z_args)
    
    # Generate class-specific predictors (including Intercept)
    X_list <- list(Intercept = rep(1, n))
    for (pred in model$predictors$class_specific) {
      dist_fn <- pred$observations$distribution
      dist_args <- pred$observations$parameters
      dist_args$n <- n
      X_list[[pred$name]] <- do.call(dist_fn, dist_args)
    }
    
    # Assemble full data frame for model.matrix()
    sim_df <- as.data.frame(X_list[names(X_list) != "Intercept"])
    
    # Get class-specific design matrix X
    X <- model.matrix(rhs_formula, data = sim_df)
    
    sim_df[[shared$name]] <- Z
    
    # Extract true parameters
    psi_0 <- theta_0$psi_0
    Beta_0 <- theta_0$Beta_0
    
    # Compute linear predictors
    eta <- get_eta(Z, X, psi_0, Beta_0)
    
    # Expand to full class set
    eta_full <- matrix(0, nrow = n, ncol = num_classes)
    eta_full[, -ref_class] <- eta
    colnames(eta_full) <- paste0("Y", 1:num_classes)
    
    # Convert to probabilities
    Y_probs <- t(vapply(1:n, function(i) softmax(eta_full[i, ]), numeric(num_classes)))
    
    # Sample outcomes
    Y_numeric <- apply(Y_probs, 1, function(p) sample(1:num_classes, 1, prob = p))
    Y_factor <- factor(Y_numeric, levels = 1:num_classes, labels = paste0("Y", 1:num_classes))
    
    class_counts <- table(Y_factor)
    
    if (length(class_counts) == num_classes && all(class_counts >= min_obs)) {
      message(sprintf("✅ Balanced dataset generated on try %d", i))
      
      model_df <- sim_df
      model_df[[model$response$name]] <- Y_factor
      
      return(list(model_df = model_df, Y_probs = Y_probs))
    }
  }
  
  stop(sprintf("❌ Failed to generate a balanced dataset after %d tries. Consider relaxing `min_obs` or adjusting your true parameters.", max_tries))
}



