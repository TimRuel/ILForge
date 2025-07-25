library(VGAM)
library(tidyverse)

set.seed(1)
n <- 500
X1 <- rnorm(n)
X2 <- rnorm(n)
X <- matrix(c(X1, X2), ncol = 2, byrow = FALSE)
Z  <- rnorm(n)

# True parameters
psi_0   <- 1.2
Beta_0  <- matrix(c(0.5, -0.5, 
                    1, -1,
                    0.5, 0.5), nrow = 3, byrow = TRUE)

eta_fn <- function(psi, beta_vec) {
  beta <- matrix(beta_vec, nrow = 3, byrow = TRUE)  # 3 rows: intercept, X1, X2
  eta <- cbind(1, X) %*% beta + psi * Z
  eta  # n × 2 matrix
}

# Linear predictors
eta <- eta_fn(psi_0, Beta_0)
eta1 <- eta[,1]
eta2 <- eta[,2]

# Probabilities (softmax with baseline class 3)
denom <- 1 + exp(eta1) + exp(eta2)
p1 <- exp(eta1) / denom
p2 <- exp(eta2) / denom
p3 <- 1 - p1 - p2

Y <- apply(cbind(p1, p2, p3), 1, function(p) sample(1:3, 1, prob = p))
Y <- factor(Y, levels = 1:3, labels = c("A", "B", "C"))
df <- data.frame(Y, Z, X1, X2)

Y_mat <- model.matrix(~ Y - 1, data = df)
Y_design <- Y_mat[, 1:2]

log_likelihood <- function(psi, lambda) {
  eta <- eta_fn(psi, lambda)
  sum(rowSums(Y_design * eta) - log(1 + rowSums(exp(eta))))
}

# Define the formula
formula <- Y ~ Z + X1 + X2

# Construct constraints matrix
# 3 response levels → 2 logits (reference = class 3)
# We will constrain the Z coefficient to be the same across both logits
# Other coefficients (X1, X2, intercept) remain unconstrained

C <- list("(Intercept)" = diag(2),           # Intercepts free
          "Z" = matrix(1, nrow = 2, ncol = 1), # Shared effect
          "X1" = diag(2),                     # Separate
          "X2" = diag(2)                      # Separate
          )

# Fit the model
fit <- vglm(formula, family = multinomial(refLevel = 3), 
            data = df, constraints = C)

summary(fit)

psi_hat <- coef(fit)["Z"]

lambda_hat <- coef(fit)[setdiff(names(coef(fit)), "Z")]

vcov_lambda <- vcov(fit)

Sigma <- vcov_lambda[-3, -3]

set.seed(NULL)

psi1_grid <- seq(0.5, 1.3, 0.05)

softmax <- function(x) exp(x) / sum(exp(x))

R <- 100

l_tilde_mat <- matrix(1,
                      nrow = R,
                      ncol = length(psi1_grid))

for (j in seq_along(psi1_grid)) {
  
  psi1 <- psi1_grid[j]
  
  print(psi1)
  
  for (i in 1:R) {
    
    phi <- MASS::mvrnorm(1, lambda_hat, Sigma) |>
      matrix(nrow = nrow(Beta_0),
             ncol = ncol(Beta_0),
             byrow = TRUE)
    
    nu <- (psi_hat * Z + cbind(1, X) %*% phi) |>
      cbind(0) |>
      apply(1, softmax) |>
      t()
    
    obj_fn <- function(lambda) -sum(rowSums(nu[, c(1,2)] * eta_fn(psi1, lambda)) - log(1 + rowSums(exp(eta_fn(psi1, lambda)))))
    
    lambda_tilde <- nloptr::auglag(
      x0 = rep(1, 6),
      fn = obj_fn,
      localsolver = "SLSQP",
      deprecatedBehavior = FALSE)$par
    
    l_tilde_mat[i,j] <- log_likelihood(psi1, lambda_tilde)
  }
}

log_L_bar <- matrixStats::colLogSumExps(l_tilde_mat, na.rm = TRUE) - log(R)

plot(psi1_grid, log_L_bar)

spline_mod <- smooth.spline(psi1_grid, log_L_bar)

MLE <- optimize(\(psi) predict(spline_mod, psi)$y,
               lower = head(psi1_grid, 1),
               upper = tail(psi1_grid, 1),
               maximum = TRUE)

log_L_bar_mod <- function(psi) predict(spline_mod, psi)$y - MLE$objective

alpha <- 0.05

crit <- qchisq(1 - alpha, df = 1) / 2

lower_bound <- uniroot(\(psi) log_L_bar_mod(psi) + crit,
                       interval = c(head(psi1_grid, 1), MLE$maximum))$root |> 
  round(3)

upper_bound <- uniroot(\(psi) log_L_bar_mod(psi) + crit,
                       interval = c(MLE$maximum, tail(psi1_grid, 1)))$root |> 
  round(3)

CI1 <- c(lower_bound, upper_bound)

psi_hat_se <- coef(summary(fit))["Z", "Std. Error"]

num_se <- qnorm(1-alpha/2)

CI2 <- (psi_hat + c(-1, 1) * psi_hat_se * num_se) |> 
  round(3)

CI1
CI2

diff(CI1)
diff(CI2)

###########################################################################
###########################################################################
###########################################################################

formula <- Y ~ Z + X1 + X2

C <- list("(Intercept)" = diag(5),           
          "Z" = matrix(1, nrow = 5, ncol = 1), 
          "X1" = diag(5),                     
          "X2" = diag(5))

ref_level <- 6

model <- fit_model(model_df, formula, C, ref_level)

Beta_MLE <- get_Beta_MLE(model)

psi_hat <- get_psi_hat(model)

X <- model_df |> 
  select(starts_with("X")) |> 
  as.matrix()
X <- cbind(1, X)

Z <- model_df$Z

Y_one_hot <- model.matrix(~ Y - 1, data = model_df)
Y_one_hot <- Y_one_hot[, 1:(ref_level-1)]

log_likelihood(X, Y_one_hot, Z, psi_hat, Beta_MLE)

mu <- get_E_Y(Z, X, psi_hat, Beta_MLE)

log_likelihood(X, mu, Z, psi_hat + 1, Beta_MLE)

alpha <- 0.05

psi_hat_se <- VGAM::summaryvglm(model)@coef3["Z", "Std. Error"]

num_se <- qnorm(1-alpha/2)

interval <- psi_hat + c(-2, 2) * psi_hat_se * num_se

vcov_Beta <- vcov(model)

Sigma <- vcov_Beta[-6, -6]

vcov_Beta[!rownames(vcov_Beta) %in% "Z", !colnames(vcov_Beta) %in% "Z"]

vcov_Beta[-"Z"]

phi <- MASS::mvrnorm(1, gdata::unmatrix(Beta_MLE, byrow=TRUE), Sigma) |>
  matrix(nrow = nrow(Beta_MLE),
         ncol = ncol(Beta_MLE),
         byrow = TRUE)

mu <- get_E_Y(Z, X, psi_hat, phi)

branch <- function(psi) {
  
  Beta_hat_obj_fn <- function(Beta) {
    
    Beta <- Beta |> 
      matrix(nrow = nrow(phi),
             ncol = ncol(phi),
             byrow = TRUE)
    
    return(-log_likelihood(X, mu, Z, psi, Beta))
  }
  
  Beta_hat <- get_Beta_hat(Beta_hat_obj_fn, gdata::unmatrix(phi, byrow=TRUE)) |> 
    matrix(nrow = nrow(phi),
           ncol = ncol(phi),
           byrow = TRUE)
  
  return(log_likelihood(X, Y_one_hot, Z, psi, Beta_hat))
}

psi_grid <- seq(interval[1], interval[2], 0.05)

log_L_tilde <- purrr::map_dbl(psi_grid, branch)

plot(psi_grid, log_L_tilde)

get_branch_mode(phi, psi_hat, Z, Y_one_hot, X, interval) 

library(gdata)
gdata::unmatrix(Beta_MLE, byrow=TRUE) |> unname()

library(foreach)
library(doFuture)
plan(multisession)
registerDoFuture()

foreach(i = 1:3) %dofuture% {
  i^2
}


config <- yaml::read_yaml("C:/Northwestern/ILForge/experiments/exp_v2.0.0/simulations/sim_01/iter_01/config_snapshot.yml")
model_df <- readRDS("C:/Northwestern/ILForge/experiments/exp_v2.0.0/simulations/sim_01/iter_01/data/model_df.rds")

