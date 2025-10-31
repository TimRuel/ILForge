args <- c(
  
  "poisson_regression",
  "weighted_sum",
  "exp_v3.0.0",
  "sim_01",
  "iter_01"
)

Beta_0 <- readRDS("C:/Northwestern/ILForge/experiments/exp_v3.0.0/true_params/Beta_0.rds")

covariates <- readRDS("C:/Northwestern/ILForge/experiments/exp_v3.0.0/true_params/covariates.rds")

weights <- readRDS("C:/Northwestern/ILForge/experiments/exp_v3.0.0/true_params/weights.rds")

lambda_0 <- readRDS("C:/Northwestern/ILForge/experiments/exp_v3.0.0/true_params/lambda_0.rds")

psi_0 <- readRDS("C:/Northwestern/ILForge/experiments/exp_v3.0.0/true_params/psi_0.rds")

data <- readRDS("C:/Northwestern/ILForge/experiments/exp_v3.0.0/simulations/sim_01/iter_01/data/data.rds")

config <- yaml::read_yaml("C:/Northwestern/ILForge/experiments/exp_v3.0.0/simulations/sim_01/iter_02/config_snapshot.yml")

compute_lambda_j <- function(alpha_j, beta_j, mu_j = NULL, Sigma_j = NULL) {
  # Ensure beta_j is a matrix (groups in rows)
  beta_j <- as.matrix(beta_j)
  n_groups <- nrow(beta_j)
  p <- ncol(beta_j)
  
  # Handle defaults for covariate distribution
  if (is.null(mu_j)) {
    mu_j <- matrix(0, n_groups, p)
  } else {
    mu_j <- as.matrix(mu_j)
  }
  
  if (is.null(Sigma_j)) {
    # Assume identity matrix for each group
    Sigma_j_list <- replicate(n_groups, diag(p), simplify = FALSE)
  } else if (is.list(Sigma_j)) {
    Sigma_j_list <- Sigma_j
  } else {
    # Shared covariance matrix for all groups
    Sigma_j_list <- replicate(n_groups, Sigma_j, simplify = FALSE)
  }
  
  # Compute lambda_j for each group
  lambda_j <- numeric(n_groups)
  for (j in seq_len(n_groups)) {
    mean_part <- sum(mu_j[j, ] * beta_j[j, ])
    var_part <- 0.5 * t(beta_j[j, ]) %*% Sigma_j_list[[j]] %*% beta_j[j, ]
    lambda_j[j] <- exp(alpha_j[j] + mean_part + var_part)
  }
  
  return(lambda_j)
}


# Example: 3 groups, 2 covariates
alpha_j <- c(0.1, 0.3, -0.2)
beta_j <- matrix(c(0.4, 0.5,
                   0.6, -0.3,
                   -0.1, 0.7), nrow = 3, byrow = TRUE)

# Assume N(0, I) covariates
lambda_shared <- compute_lambda_j(alpha_j, beta_j)

# Group-specific means and shared covariance
mu_j <- matrix(c(0, 0,
                 1, 1,
                 -1, 0), nrow = 3, byrow = TRUE)
Sigma <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)

lambda_custom <- compute_lambda_j(alpha_j, beta_j, mu_j = mu_j, Sigma_j = Sigma)

library(tidyverse)
library(MASS)

# --- Step 1: Global Parameters ---
J <- 5                 # Number of groups
n_per_group <- 100      # Observations per group
p <- 2                  # Number of covariates

# --- Step 2: Generate Group-Level Parameters ---
set.seed(123)

alpha_j <- rnorm(J, mean = 0, sd = 0.5)
beta_j  <- replicate(J, rnorm(p, mean = 0.5, sd = 0.2), simplify = FALSE)
mu_j     <- replicate(J, runif(p, -1, 1), simplify = FALSE)
Sigma_j  <- replicate(J, diag(runif(p, 0.2, 1.0)), simplify = FALSE)

# --- Step 3: Compute True Group-Level Rates (λ_j) ---
lambda_j <- map_dbl(seq_len(J), \(j) {
  beta <- beta_j[[j]]
  mu   <- mu_j[[j]]
  Sigma <- Sigma_j[[j]]
  eta_bar <- alpha_j[j] + sum(beta * mu) + 0.5 * t(beta) %*% Sigma %*% beta
  as.numeric(exp(eta_bar))
})

# --- Step 4: Simulate Data per Group ---
sim_data <- map_dfr(seq_len(J), \(j) {
  n <- n_per_group
  X <- MASS::mvrnorm(n, mu = mu_j[[j]], Sigma = Sigma_j[[j]])
  exposure <- runif(n, 0.5, 2.0)
  eta <- alpha_j[j] + X %*% beta_j[[j]]
  lambda_ij <- exp(eta)
  y <- rpois(n, lambda = lambda_ij * exposure)
  
  tibble(
    group = j,
    y = y,
    exposure = exposure,
    !!!set_names(as.data.frame(X), paste0("X", 1:p))
  )
})

library(glmmTMB)

# Convert group to factor if not already
sim_data$group <- factor(sim_data$group)

# Build the formula: allow for random intercept and slopes for all X vars
formula <- as.formula("y ~ X1 + X2 + (1 + X1 + X2 | group) + offset(log(exposure))")

# Fit the Poisson mixed model
fit <- glmmTMB(formula, family = poisson(), data = sim_data)

# View summary
coef(fit)

logLik(fit)



