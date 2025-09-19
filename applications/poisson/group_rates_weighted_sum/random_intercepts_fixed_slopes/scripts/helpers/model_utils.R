# applications/poisson/group_rates_weighted_sum/random_intercepts_fixed_slopes/scripts/helpers/model_utils.R

library(Matrix)

make_formula <- function(config) {
  # Build RHS parts
  rhs_parts <- c()
  
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
  
  rhs_parts <- c(rhs_parts, "(1 | group)")
  
  # Add group main effect (intercepts per group)
  rhs <- paste(rhs_parts, collapse = " + ")
  
  # Full formula
  f <- as.formula(paste("Y ~", rhs))
  
  return(f)
}

fit_model <- function(config, data) {
  
  formula <- make_formula(config)
  
  glmer(formula, offset = log(t), data = data, family = poisson)
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
  new_names[match(intercepts, coef_names)] <- "μ_α"
  
  # Fixed slopes → γ#
  new_names[match(main_effects, coef_names)] <- 
    paste0("γ", gsub("^X", "", main_effects))
  
  # Varying slopes → ζ#_group
  for (int in interactions) {
    parts <- strsplit(int, ":")[[1]]   # e.g., c("groupA", "X2")
    xnum  <- sub("^X", "", parts[1])
    group <- sub("^group", "", parts[2])
    new_names[match(int, coef_names)] <- paste0("ζ", xnum, "_", group)
  }
  
  return(new_names)
}

get_Beta_MLE <- function(model) {
  
  fixed_effects <- fixef(model)
  Beta_MLE <- as.matrix(fixed_effects, ncol = 1)  # stack as single-column matrix
  
  rownames(Beta_MLE) <- rename_coefs(names(fixed_effects))
  
  Beta_MLE
}

get_sigma_alpha_MLE <- function(model) unname(attr(VarCorr(model)[["group"]], "stddev"))

get_eta <- function(Beta, X, alpha_g) X %*% Beta + alpha_g

get_theta <- function(Beta, X, alpha_g) exp(get_eta(Beta, X, alpha_g))

get_theta_g <- function(Beta, sigma_alpha, X) {
  
  data.frame(exp_eta = exp(X %*% Beta),
             group = factor(rownames(X))) |> 
    group_by(group) |> 
    summarise(theta = mean(exp_eta) * exp(sigma_alpha^2 / 2)) |> 
    deframe()
}

get_mu <- function(eta, t) t * exp(eta)

get_psi <- function(Beta, sigma_alpha, X, weights) {
  
  theta_g <- get_theta_g(Beta, sigma_alpha, X)
  
  sum(theta_g * weights)
}

# General delta-method version of SE(psi_hat)
get_psi_MLE_SE <- function(model, weights) {
  # Extract fixed-effects estimates
  Beta_MLE <- fixef(model)
  Beta_cov <- vcov(model)
  
  # Extract sigma_alpha (random intercept std dev)
  sigma_alpha <- as.numeric(attr(VarCorr(model)[[1]], "stddev"))
  
  # Design matrix and group labels
  X <- model.matrix(model)
  groups <- model@frame$group
  group_labels <- levels(groups)
  
  # Ensure weights are a named vector
  if (is.null(names(weights))) {
    if (length(weights) != length(group_labels)) {
      stop("weights must have length equal to number of groups")
    }
    weights <- setNames(weights, group_labels)
  }
  
  # Linear predictor and exponentiated values
  eta <- drop(X %*% Beta_MLE)
  exp_eta <- exp(eta)
  
  # Compute per-group means
  n_g <- table(groups)
  theta_g <- tapply(exp_eta, groups, mean)
  
  # Compute psi_hat
  psi_hat <- sum(weights * theta_g) * exp(sigma_alpha^2 / 2)
  
  # Gradient w.r.t Beta
  obs_w <- weights[as.character(groups)] / n_g[as.character(groups)]
  grad_beta <- as.numeric(t(X) %*% (obs_w * exp_eta)) * exp(sigma_alpha^2 / 2)
  
  # Gradient w.r.t sigma_alpha^2
  grad_sigma <- 0.5 * psi_hat
  
  # Combine gradients
  grad <- c(grad_beta, grad_sigma)
  
  # Variance-covariance matrix for (Beta, sigma_alpha^2)
  # lme4 does not provide Var(sigma_alpha^2), so we assume independence and zero for simplicity
  V <- bdiag(Beta_cov, Matrix(0, 1, 1))
  
  # Delta-method SE
  SE <- sqrt(as.numeric(t(grad) %*% V %*% grad))
  
  return(SE)
}

# Define log integrand Q_g(α)
Qg <- function(Beta, sigma_alpha, alpha_g, Xg, Yg, tg) {
  
  eta_g <- get_eta(Beta, Xg, alpha_g)
  mu_g <- get_mu(eta_g, tg)
  
  ll_cond <- sum(Yg * (log(tg) + eta_g) - mu_g - lgamma(Yg + 1))
  prior <- -0.5 * (alpha_g^2 / sigma_alpha^2) - 0.5 * log(2 * pi * sigma_alpha^2)
  ll_cond + prior
}

# First derivative of Q_g
Qg_prime <- function(Beta, sigma_alpha, alpha_g, Xg, Yg, tg) {
  sum(Yg - tg * exp(as.numeric(Xg %*% Beta + alpha_g))) - alpha_g / sigma_alpha^2
}

log_likelihood <- function(Beta, sigma_alpha, X, Y, t) {
  
  groups <- factor(rownames(X))
  group_levels <- levels(groups)
  
  logLik_val <- 0
  
  for (g in group_levels) {
    # Subset data for this group
    idx <- which(groups == g)
    Xg <- X[idx, , drop = FALSE]
    Yg <- Y[idx]
    tg <- t[idx]
    
    # Find conditional mode α̂_g by root-finding
    alpha_g_hat <- uniroot(function(alpha_g) Qg_prime(Beta, sigma_alpha, alpha_g, Xg, Yg, tg), interval = c(-10, 10))$root
    
    # Second derivative at mode
    Hg <- sum(tg * exp(as.numeric(Xg %*% Beta + alpha_g_hat))) + 1 / sigma_alpha^2
    
    # Laplace approx for group contribution
    logLik_val <- logLik_val + Qg(Beta, sigma_alpha, alpha_g_hat, Xg, Yg, tg) + 0.5 * log(2 * pi / Hg)
  }
  
  return(logLik_val)
}

likelihood <- function(Beta, sigma_alpha, X, Y, t) exp(log_likelihood(Beta, sigma_alpha, X, Y, t))
