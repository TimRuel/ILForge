args <- c(
  
  "negative_binomial",
  "process_rates_weighted_sum",
  "naive_rates",
  "exp_v7.0.1",
  "sim_01",
  "iter_01"
)

phi_0 <- readRDS("C:/Northwestern/ILForge/experiments/exp_v7.0.1/true_params/phi_0.rds")

theta_0 <- readRDS("C:/Northwestern/ILForge/experiments/exp_v7.0.1/true_params/theta_0.rds")

weights <- readRDS("C:/Northwestern/ILForge/experiments/exp_v7.0.1/true_params/weights.rds")

n_per_process <- readRDS("C:/Northwestern/ILForge/experiments/exp_v7.0.1/true_params/n_per_process.rds")

data <- readRDS("C:/Northwestern/ILForge/experiments/exp_v7.0.1/simulations/sim_01/iter_01/data/data.rds")

set.seed(123)

# Common mean rate
lambda <- 10

# Sample size
n <- 10000

# Generate Poisson data
y_pois <- rpois(n, lambda = lambda)

# Generate Negative Binomial data with same mean
# The NB mean = mu, variance = mu + mu^2 / size
# A smaller 'size' increases dispersion
size <- 2  # high dispersion (smaller = more overdispersed)
y_nb <- rnbinom(n, size = size, mu = lambda)

# Compare summary stats
cat("Poisson variance:", var(y_pois), "\n")
cat("NegBin variance :", var(y_nb), "\n")

# Plot histograms side by side
hist(
  y_pois, breaks = 40, col = rgb(0.2, 0.4, 0.8, 0.5),
  xlim = c(0, 30), main = "Poisson vs Negative Binomial (same mean)",
  xlab = "Count", freq = FALSE
)
hist(
  y_nb, breaks = 40, col = rgb(0.8, 0.3, 0.3, 0.4),
  add = TRUE, freq = FALSE
)

legend(
  "topright", legend = c("Poisson", "Neg. Binomial"),
  fill = c(rgb(0.2, 0.4, 0.8, 0.5), rgb(0.8, 0.3, 0.3, 0.4))
)
