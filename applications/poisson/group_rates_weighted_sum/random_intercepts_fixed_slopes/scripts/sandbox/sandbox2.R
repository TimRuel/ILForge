# Load package
library(lme4)

set.seed(123)

# Simulation settings
n_groups <- 30          # number of groups
n_per_group <- 20       # observations per group
beta0 <- 1.0            # fixed intercept
beta1 <- 0.5            # slope for x
sigma_b <- 0.7          # SD of random intercepts

# Generate group-level random intercepts
b <- rnorm(n_groups, mean = 0, sd = sigma_b)

# Data frame
group <- rep(1:n_groups, each = n_per_group)
x <- rnorm(n_groups * n_per_group, mean = 0, sd = 1)
eta <- beta0 + beta1 * x + b[group]
mu <- exp(eta)
y <- rpois(length(mu), lambda = mu)

dat <- data.frame(y, x, group = factor(group))

# Fit random intercepts Poisson regression
fit <- glmer(y ~ x + (1 | group), data = dat, family = poisson)

summary(fit)

library(lme4)

set.seed(123)

# ----------------------------
# Simulation setup
# ----------------------------
G <- 5             # number of groups
n_per_group <- 200
N <- G * n_per_group

alpha_mean <- 0.5
alpha_sd   <- 0.3

gamma1     <- 0.4   # homogeneous slope (fixed, shared across groups)
gamma2     <- -0.2  # homogeneous slope (fixed, shared across groups)

# heterogeneous slopes for X3, one per group
gamma3 <- rnorm(G, mean = 0.6, sd = 0.1)

# group IDs
group <- factor(rep(1:G, each = n_per_group))

# random intercepts for groups
alpha <- rnorm(G, mean = alpha_mean, sd = alpha_sd)

# covariates
X1 <- rnorm(N, 0, 1)   # homogeneous slope
X2 <- runif(N, -1, 1)  # homogeneous slope
X3 <- rnorm(N, 0, 1)   # heterogeneous slope (different coefficient by group)

# linear predictor
eta <- alpha[group] + gamma1 * X1 + gamma2 * X2 + gamma3[group] * X3

# Poisson mean and response
lambda <- exp(eta)
Y <- rpois(N, lambda)

# ----------------------------
# Fit model
# ----------------------------
# Trick: interaction of X3 with group creates group-specific slopes
dat <- data.frame(Y, X1, X2, X3, group)

fit <- glmer(Y ~ X1 + X2 + X3:group + (1 | group),
             family = poisson,
             data = dat)

summary(fit)

# ----------------------------
# Extract estimates
# ----------------------------
fixef(fit)        # γ1, γ2, plus γ3 for each group
ranef(fit)$group  # α_g deviations

library(lme4)

set.seed(123)

# ------------------------
# Simulation setup
# ------------------------
G <- 10                     # number of groups
n_per_group <- 50           # observations per group
N <- G * n_per_group

# True random intercept distribution
mu_alpha <- 1.5
sigma_alpha <- 0.5
alpha_g <- rnorm(G, mean = mu_alpha, sd = sigma_alpha)  # true group intercepts

# One fixed slope
beta <- 0.3

# Covariate
x <- rnorm(N, 0, 1)

# Group labels
group <- rep(1:G, each = n_per_group)

# Linear predictor and Poisson outcome
eta <- alpha_g[group] + beta * x
lambda <- exp(eta)
y <- rpois(N, lambda)

dat <- data.frame(y, x, group = factor(group))

# ------------------------
# Fit mixed model
# ------------------------
fit <- glmer(y ~ x + (1 | group), data = dat, family = poisson)

fit <- glmer(Y ~ X1 + X2:group + (1 | group), offset = log(t), data = data, family = poisson)

sigma_alpha_hat <- unname(attr(VarCorr(fit)[["group"]], "stddev"))

data |> 
  mutate(mu_hat = predict(fit, type = "response", re.form = NA),
         theta_hat = mu_hat / t) |> 
  group_by(group) |> 
  summarise(theta_hat = mean(theta_hat)) |> 
  deframe() |> 
  (\(x) x * exp(sigma_alpha_hat^2 / 2))()

# ------------------------
# Compare estimates
# ------------------------
cat("True mean of group intercepts (mu_alpha):", mu_alpha, "\n")
cat("Sample mean of simulated alpha_g:", mean(alpha_g), "\n")
cat("Estimated fixed intercept:", fixef(fit)["(Intercept)"], "\n")
cat("Estimated slope:", fixef(fit)["x"], "\n")

# Inspect random intercept deviations
head(ranef(fit)$group)

| Class                                                        | Intercept | Slopes                                | Random effects to integrate        | Covariate integration | Notes                                                                                                     |
  | ------------------------------------------------------------ | --------- | ------------------------------------- | ---------------------------------- | --------------------- | --------------------------------------------------------------------------------------------------------- |
  | **A — No random effects**                                    | Fixed     | All fixed                             | None                               | Yes                   | Profile likelihood = marginal likelihood; integrate only over covariates for true marginal rate           |
  | **B — Random intercepts only**                               | Random    | All fixed                             | $\alpha_g$                         | Yes                   | Standard random intercept model; marginalize over intercepts and covariates                               |
  | **C — Random intercept + fixed slopes (some heterogeneous)** | Random    | Some fixed slopes vary by group       | $\alpha_g$                         | Yes                   | Heterogeneous slopes are fixed; integrate over random intercepts and covariates                           |
  | **D — Random slopes only**                                   | Fixed     | All slopes random                     | $b_{gk}$ for all slopes            | Yes                   | No random intercept; integrate over slopes and covariates                                                 |
  | **E — Random intercept + random slopes (all random)**        | Random    | All slopes random                     | $\alpha_g, b_{gk}$                 | Yes                   | Full hierarchical model; integrate over intercept, slopes, and covariates                                 |
  | **F — Random intercept + mixed slopes**                      | Random    | Some slopes fixed, some slopes random | $\alpha_g, b_{gk}^{\text{random}}$ | Yes                   | Hybrid: integrate over random intercept, random slopes, and covariates; fixed slopes treated as constants |
  | **G — Fixed intercept + mixed slopes**                       | Fixed     | Some slopes fixed, some slopes random | $b_{gk}^{\text{random}}$           | Yes                   | Hybrid without random intercept; integrate over random slopes and covariates                              |
  
  
  | Class                                                 | Intercept | Slopes                                                  | Random effects to integrate        | Covariate integration | Notes                                                                                            |
  | ----------------------------------------------------- | --------- | ------------------------------------------------------- | ---------------------------------- | --------------------- | ------------------------------------------------------------------------------------------------ |
  | **A — No random effects**                             | Fixed     | All fixed (shared or group-specific)                    | None                               | Yes                   | Standard regression; only integrate over covariates for true marginal rate                       |
  | **B — Random intercepts with fixed slopes**           | Random    | All slopes fixed (some homogeneous, some heterogeneous) | $\alpha_g$                         | Yes                   | Heterogeneous fixed slopes treated as constants; integrate over random intercepts and covariates |
  | **E — Random intercept + random slopes (all random)** | Random    | All slopes random                                       | $\alpha_g, b_{gk}$                 | Yes                   | Full hierarchical model; integrate over intercept, slopes, and covariates                        |
  | **F — Random intercept + mixed slopes**               | Random    | Some slopes fixed, some slopes random                   | $\alpha_g, b_{gk}^{\text{random}}$ | Yes                   | Hybrid model; integrate over random intercept + random slopes; fixed slopes are constants        |
  

