# applications/multinomial_logistic_regression/shared_effect/scripts/helpers/data_utils.R

get_eta <- function(Z, X, psi, Beta) psi * Z + X %*% Beta

softmax <- function(x) exp(x) / sum(exp(x))

softmax_adj <- function(x) exp(x) / (1 + sum(exp(x)))

generate_data <- function(config, theta_0) {
  
  n <- config$data_generation$n_obs
  model <- config$model
  num_classes <- model$response$num_classes
  num_effects <- model$response$num_effects
  ref_class <- model$response$reference_class
  
  # Shared predictor
  dist <- model$predictors$shared$distribution
  print(model)
  if (dist$type == "normal") {
    Z <- rnorm(n, mean = dist$parameters$mean, sd = dist$parameters$sd)
  } else {
    stop("Unsupported distribution type for shared predictor: ", dist$type)
  }
  
  # Class-specific predictors
  # Generate each class-specific predictor as a vector of length n
  X_list <- list(Intercept = rep(1, n))  # Intercept row is all 1s
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

  X <- do.call(cbind, X_list)
  
  psi_0 <- theta_0$psi_0
  Beta_0 <- theta_0$Beta_0
  
  eta <- get_eta(Z, X, psi_0, Beta_0)
  
  eta_full <- cbind(eta, rep(0, n))  
  Y_probs <- softmax(eta_full)
  
  Y_numeric <- apply(Y_probs, 1, function(p) sample(1:num_classes, 1, prob = p))
  Y_factor <- factor(Y_numeric, levels = 1:num_classes)
  
  model_df <- as.data.frame(X_list[-which(names(X_list) == "Intercept")])
  model_df$Z <- Z
  model_df$Y <- Y_factor
  
  return(list(
    model_df = model_df,
    Y_probs = Y_probs
  ))
}
