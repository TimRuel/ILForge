# applications/poisson/group_rates_weighted_sum/fixed_effects_regression/scripts/helpers/model_utils.R

make_formula <- function(config) {
  # Build RHS parts
  rhs_parts <- c("0", "group")
  
  covs <- config$model$covariates
  homo_covs <- Filter(\(c) c$type == "homogeneous", covs)
  hetero_covs <- Filter(\(c) c$type == "heterogeneous", covs)
  
  # Homogeneous covariates (same slope across groups)
  if (length(homo_covs) > 0) {
    homo_covs <- sapply(homo_covs, \(x) x$variable$symbol)
    rhs_parts <- c(rhs_parts, homo_covs)
  }
  
  # Heterogeneous covariates (different slope across groups, can be interpreted as interaction with group)
  if (length(hetero_covs) > 0) {
    hetero_covs <- sapply(hetero_covs, \(x) x$variable$symbol)
    rhs_parts <- c(rhs_parts, paste0(hetero_covs, ":group"))
  }
  
  # Add group main effect (intercepts per group)
  rhs <- paste(rhs_parts, collapse = " + ")
  
  # Full formula
  f <- as.formula(paste("Y ~", rhs))
  
  return(f)
}

fit_model <- function(config, data) {
  
  formula <- make_formula(config)
  
  glm(formula, offset = log(t), family = poisson(), data = data)
}

rename_coefs <- function(coef_names) {
  # Separate intercepts, main effects, and interactions
  intercepts   <- coef_names[!grepl("X", coef_names)]              # group names
  main_effects <- coef_names[grepl("^X\\d+$", coef_names)]         # X1, X2, ...
  interactions <- coef_names[grepl(":", coef_names)]               # group:X2, ...
  
  # Figure out which X's are fixed vs varying
  # Fixed slopes: main effects only
  fixed_Xs <- gsub("^X", "", main_effects)
  # Varying slopes: from interactions
  varying_Xs <- unique(gsub(".*:X", "", interactions))
  
  # Replacements
  new_names <- character(length(coef_names))
  
  # Intercepts → α_group
  new_names[match(intercepts, coef_names)] <- paste0("α_", sub("^group", "", intercepts))
  
  # Fixed slopes → γ#
  new_names[match(main_effects, coef_names)] <- 
    paste0("γ", gsub("^X", "", main_effects))
  
  # Varying slopes → ζ#_group
  for (int in interactions) {
    parts <- strsplit(int, ":")[[1]]   # e.g., c("groupA", "X2")
    group <- sub("^group", "", parts[1])
    xnum  <- sub("^X", "", parts[2])
    new_names[match(int, coef_names)] <- paste0("ζ", xnum, "_", group)
  }
  
  return(new_names)
}

get_Beta_MLE <- function(model) {
  
  coefs <- coef(model)
  Beta_MLE <- as.matrix(coefs, ncol = 1)  # stack as single-column matrix
  
  rownames(Beta_MLE) <- rename_coefs(names(coefs))
  
  Beta_MLE
}

get_eta <- function(Beta, X) X %*% Beta

get_theta <- function(Beta, X) exp(get_eta(Beta, X))

get_theta_g <- function(Beta, X) {
  
  data.frame(theta = get_theta(Beta, X),
             group = factor(rownames(X))) |> 
    group_by(group) |> 
    summarise(theta = mean(theta)) |> 
    deframe()
}

get_mu <- function(eta, t) t * exp(eta)

get_psi <- function(Beta, X, weights) {
  
  theta <- get_theta_g(Beta, X)
  
  sum(theta * weights)
}

# General delta-method version of SE(psi_hat)
get_psi_MLE_SE <- function(model, weights, X) {
  # Extract Beta_MLE as numeric vector
  Beta_MLE <- drop(get_Beta_MLE(model))
  Beta_cov <- vcov(model)
  
  # Group membership for each observation
  groups <- model$data$group
  group_labels <- levels(groups)
  
  # Ensure weights are a named vector (one per group)
  if (is.null(names(weights))) {
    if (length(weights) != length(group_labels)) {
      stop("If 'weights' is unnamed, it must have length equal to number of groups.")
    }
    weights <- setNames(weights, group_labels)
  }
  
  # Linear predictors and exp
  eta <- drop(get_eta(Beta_MLE, X))
  exp_eta <- exp(eta)
  
  # Group sizes
  n_g <- table(groups)
  
  # Per-observation effective weights = w_g / n_g
  obs_w <- weights[as.character(groups)] / n_g[as.character(groups)]
  
  # Gradient: sum_i obs_w[i] * exp(eta_i) * x_i
  grad <- as.numeric(t(X) %*% (obs_w * exp_eta))
  names(grad) <- colnames(X)
  
  # Delta-method SE
  se <- sqrt(as.numeric(t(grad) %*% Beta_cov %*% grad))
  
  return(se)
}

log_likelihood <- function(Beta, X, Y, t) {
  
  eta <- get_eta(Beta, X)
  
  mu <- get_mu(eta, t)
  
  sum(Y * (log(t) + eta) - mu - lgamma(Y+1))
}

likelihood <- function(Beta, X, Y, t) exp(log_likelihood(Beta, X, Y, t))
