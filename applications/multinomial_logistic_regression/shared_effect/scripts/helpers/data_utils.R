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
  min_prop <- config$data_generation$min_prop
  min_required <- ceiling(n * min_prop)
  max_tries <- config$data_generation$max_tries %||% 100
  
  for (i in 1:max_tries) {
    
    # Shared predictor (Z)
    shared <- model$predictors$shared
    Z_dist <- shared$observations
    Z_args <- Z_dist$parameters
    Z_args$n <- n
    Z <- do.call(Z_dist$distribution, Z_args)
    
    # Class-specific predictors (X-)
    X_list <- list(Intercept = rep(1, n))
    for (pred in model$predictors$class_specific) {
      dist_fn <- pred$observations$distribution
      dist_args <- pred$observations$parameters
      dist_args$n <- n
      X_list[[pred$name]] <- do.call(dist_fn, dist_args)
    }
    
    # Construct design matrix (includes intercept)
    X <- do.call(cbind, X_list)
    colnames(X) <- names(X_list)
    
    # Pull true parameters
    psi_0 <- theta_0$psi_0
    Beta_0 <- theta_0$Beta_0
    
    # Linear predictor for each class (excluding reference)
    eta <- get_eta(Z, X, psi_0, Beta_0)
    
    # Expand to full K-class system with reference category
    eta_full <- matrix(0, nrow = n, ncol = num_classes)
    eta_full[, -ref_class] <- eta
    colnames(eta_full) <- paste0("Y", 1:num_classes)
    
    # Convert logits to probabilities
    Y_probs <- t(vapply(1:n, function(i) softmax(eta_full[i, ]), numeric(num_classes)))
    
    # Sample Y using probabilities
    Y_numeric <- apply(Y_probs, 1, function(p) sample(1:num_classes, 1, prob = p))
    Y_factor <- factor(Y_numeric, levels = 1:num_classes, labels = paste0("Y", 1:num_classes))
    
    class_counts <- table(Y_factor)

    if (length(class_counts) == num_classes) {
      if (all(class_counts >= min_required)) {
        message(sprintf("✅ Balanced dataset generated on try %d", i))

        # Create dataframe with predictors + Y (exclude Intercept from model_df)
        model_df <- as.data.frame(X_list[names(X_list) != "Intercept"])
        model_df[[shared$name]] <- Z
        model_df[[model$response$name]] <- Y_factor
        data <- list(model_df = model_df,
                     Y_probs = Y_probs)
        return(data)
      }
    }
  }
  stop(sprintf("❌ Failed to generate a balanced dataset after %d tries. Consider relaxing `min_prop` or adjusting your true parameters.", max_tries))
}

