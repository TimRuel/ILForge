# applications/multinomial_logistic_regression/shared_effect/scripts/helpers/data_utils.R

generate_data <- function(config, Beta_0) {
  
  n <- config$data_generation$n_obs
  model <- config$model
  num_classes <- model$response$num_classes
  num_effects <- model$response$num_effects
  ref_class <- model$response$reference_class
  
  # Generate each predictor as a vector of length n
  X_list <- list(Intercept = rep(1, n))  # Intercept row is all 1s
  
  # Shared predictors
  if (!is.null(model$predictors$shared)) {
    for (pred in model$predictors$shared) {
      dist <- pred$distribution
      name <- pred$name
      if (dist$type == "normal") {
        x <- rnorm(n, mean = dist$parameters$mean, sd = dist$parameters$sd)
      } else {
        stop("Unsupported distribution type for shared predictor: ", dist$type)
      }
      X_list[[name]] <- x
    }
  }
  
  # Class-specific predictors
  if (!is.null(model$predictors$class_specific)) {
    for (pred in model$predictors$class_specific) {
      dist <- pred$distribution
      name <- pred$name
      if (dist$type == "normal") {
        x <- rnorm(n, mean = dist$parameters$mean, sd = dist$parameters$sd)
      } else {
        stop("Unsupported distribution type for class-specific predictor: ", dist$type)
      }
      X_list[[name]] <- x
    }
  }
  
  # Construct full design matrix X (n x p)
  X_design <- do.call(cbind, X_list)
  
  # Compute eta = X_design %*% Beta_0  (n x (J-1))
  eta <- X_design %*% Beta_0
  
  # Append baseline logits (0s) and apply softmax
  eta_full <- cbind(eta, rep(0, n))  # class 6 = reference
  Y_probs <- exp(eta_full)
  Y_probs <- Y_probs / rowSums(Y_probs)
  
  # Sample Y based on probs
  Y_numeric <- apply(Y_probs, 1, function(p) sample(1:num_classes, 1, prob = p))
  Y_factor <- factor(Y_numeric, levels = 1:num_classes)
  
  model_df <- as.data.frame(X_list[-which(names(X_list) == "Intercept")])
  model_df$Y <- Y_factor
  
  return(list(
    model_df = model_df,
    X_design = X_design,
    Y_probs = Y_probs
  ))
}
