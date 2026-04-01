# =====================================================
# Application 4.1 – Combining Mean and Median for Symmetric Distributions
# Replication of Table 1 (MSEs) and Table 2 (Coverage rates)
# =====================================================

# -------------------------
# Required libraries
# -------------------------
set.seed(123)
library(MASS)       # For simulation and basic distributions
library(stats)      # For density, mean, median


# -------------------------
# Simulation parameters
# -------------------------
nsim <- 10000           # Number of Monte Carlo simulations
B <- 1000               # Bootstrap resamples
sample_sizes <- c(30, 50, 100)  # Sample sizes for the study

# -------------------------
# Bandwidth estimation using Silverman's rule-of-thumb
# -------------------------
silverman_bandwidth <- function(x) {
  1.06 * sd(x) * length(x)^(-1/5)
}

# -------------------------
# Kernel density estimate at a specific point x0
# -------------------------
kernel_density <- function(x, x0, h) {
  n <- length(x)
  sum(dnorm((x - x0) / h)) / (n * h)
}

# -------------------------
# Laplace-based averaging estimator (AV)
# Combines mean and median using theoretical optimal weights
# -------------------------
estimate_AV <- function(x) {
  n <- length(x)
  x0 <- median(x)
  s2 <- var(x)
  m_hat <- mean(abs(x - x0))                  # Estimate of E|X - theta|
  h <- silverman_bandwidth(x)
  f_hat <- kernel_density(x, x0, h)           # Density estimate at x0 = median
  p1 <- 1 / (4 * f_hat) - m_hat / 2
  p2 <- s2 * f_hat - m_hat / 2
  alpha1 <- p1 / (p1 + p2)
  alpha2 <- 1 - alpha1
  return(alpha1 * mean(x) + alpha2 * x0)
}

# -------------------------
# Bootstrap-based averaging estimator (AVB)
# Combines mean and median using optimal weights from estimated covariance
# -------------------------
estimate_AVB <- function(x, B = 1000) {
  n <- length(x)
  mean_x <- mean(x)
  median_x <- median(x)
  means <- numeric(B)
  medians <- numeric(B)
  for (b in 1:B) {
    x_b <- sample(x, size = n, replace = TRUE)
    means[b] <- mean(x_b)
    medians[b] <- median(x_b)
  }
  v11 <- mean((means - mean_x)^2)
  v12 <- mean((means - mean_x) * (medians - median_x))
  v22 <- mean((medians - median_x)^2)
  Sigma_hat <- matrix(c(v11, v12, v12, v22), 2, 2)
  w <- solve(Sigma_hat, rep(1, 2))
  w <- w / sum(w)                               # Normalized weights
  return(as.numeric(w[1] * mean_x + w[2] * median_x))
}

# -------------------------
# Mixture distribution: 50% N(-2,1) and 50% N(2,1)
# -------------------------
rmix <- function(n) {
  ifelse(runif(n) < 0.5, rnorm(n, -2, 1), rnorm(n, 2, 1))
}

# Simulation wrapper
simulate_distribution <- function(rdist, name) {
  results <- list()
  for (n in sample_sizes) {
    mse_mean <- mse_median <- mse_AV <- mse_AVB <- numeric(nsim)
    for (i in 1:nsim) {
      x <- rdist(n)
      m1 <- mean(x)
      m2 <- median(x)
      av <- estimate_AV(x)
      avb <- estimate_AVB(x, B)
      theta <- 0 # true center of symmetry
      mse_mean[i] <- (m1 - theta)^2
      mse_median[i] <- (m2 - theta)^2
      mse_AV[i] <- (av - theta)^2
      mse_AVB[i] <- (avb - theta)^2
    }
    results[[paste0("n", n)]] <- list(
      mean = 100 * mean(mse_mean),
      median = 100 * mean(mse_median),
      av = 100 * mean(mse_AV),
      avb = 100 * mean(mse_AVB),
      sd_mean = 100 * sd(mse_mean) / sqrt(nsim),
      sd_median = 100 * sd(mse_median) / sqrt(nsim),
      sd_av = 100 * sd(mse_AV) / sqrt(nsim),
      sd_avb = 100 * sd(mse_AVB) / sqrt(nsim)
    )
  }
  return(results)
}


# -------------------------
# List of symmetric distributions to simulate
# -------------------------
distributions <- list(
  Cauchy = function(n) rcauchy(n),
  Student4 = function(n) rt(n, df = 4),
  Student7 = function(n) rt(n, df = 7),
  Logistic = function(n) rlogis(n),
  Gaussian = function(n) rnorm(n),
  Mixture = function(n) rmix(n)
)

# -------------------------
# Run simulations for Table 1 (MSE and Std)
# -------------------------
all_results <- list()
for (distname in names(distributions)) {
  cat("Simulating:", distname, "\n")
  all_results[[distname]] <- simulate_distribution(distributions[[distname]], distname)
}

# -------------------------
# Display results for Table 1
# -------------------------
for (n in sample_sizes) {
  cat("====== n =", n, "======\n")
  for (distname in names(all_results)) {
    res <- all_results[[distname]][[paste0("n", n)]]
    cat(sprintf("%-10s | Mean: %6.2f (%5.2f) | Median: %6.2f (%5.2f) | AV: %6.2f (%5.2f) | AVB: %6.2f (%5.2f)\n",
                distname,
                res$mean, res$sd_mean,
                res$median, res$sd_median,
                res$av, res$sd_av,
                res$avb, res$sd_avb))
  }
  cat("\n")
}

# =====================================================
# Table 2: Confidence Intervals and Empirical Coverage Rates
# =====================================================

# Computes the empirical coverage of 95% CI for AV and AVB
compute_coverage <- function(n, rdist, nsim = 10000, B = 1000, alpha = 0.05) {
  cover_AV <- cover_AVB <- numeric(nsim)
  z_alpha <- qnorm(1 - alpha / 2)
  
  for (i in 1:nsim) {
    x <- rdist(n)
    x0 <- median(x)
    
    # AV estimator and CI
    h <- silverman_bandwidth(x)
    f_hat <- kernel_density(x, x0, h)
    s2 <- var(x)
    m_hat <- mean(abs(x - x0))
    p1 <- 1 / (4 * f_hat) - m_hat / 2
    p2 <- s2 * f_hat - m_hat / 2
    w1 <- p1 / (p1 + p2)
    w2 <- 1 - w1
    av <- w1 * mean(x) + w2 * x0
    var_AV <- (w1^2) * var(x) / n + (w2^2) * (pi / 2) * var(x) / n
    lower_AV <- av - z_alpha * sqrt(var_AV)
    upper_AV <- av + z_alpha * sqrt(var_AV)
    cover_AV[i] <- as.numeric(0 >= lower_AV & 0 <= upper_AV)
    
    # AVB estimator and CI via bootstrap
    means <- medians <- numeric(B)
    for (b in 1:B) {
      xb <- sample(x, n, replace = TRUE)
      means[b] <- mean(xb)
      medians[b] <- median(xb)
    }
    v11 <- mean((means - mean(x))^2)
    v12 <- mean((means - mean(x)) * (medians - median(x)))
    v22 <- mean((medians - median(x))^2)
    Sigma <- matrix(c(v11, v12, v12, v22), 2, 2)
    w <- solve(Sigma, rep(1, 2))
    w <- w / sum(w)
    avb <- w[1] * mean(x) + w[2] * median(x)
    var_AVB <- t(w) %*% Sigma %*% w
    lower_AVB <- avb - z_alpha * sqrt(var_AVB)
    upper_AVB <- avb + z_alpha * sqrt(var_AVB)
    cover_AVB[i] <- as.numeric(0 >= lower_AVB & 0 <= upper_AVB)
  }
  
  return(c(mean(cover_AV) * 100, mean(cover_AVB) * 100))
}

# -------------------------
# Compute coverage results (Table 2)
# -------------------------
coverage_results <- matrix(NA, nrow = length(distributions), ncol = length(sample_sizes) * 2)
rownames(coverage_results) <- names(distributions)
colnames(coverage_results) <- c("AV_n30", "AVB_n30", "AV_n50", "AVB_n50", "AV_n100", "AVB_n100")

for (dist_name in names(distributions)) {
  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]
    cat("Coverage for", dist_name, "n =", n, "\n")
    coverage <- compute_coverage(n, distributions[[dist_name]], nsim = nsim, B = B)
    coverage_results[dist_name, (2*j - 1):(2*j)] <- round(coverage, 2)
  }
}

# -------------------------
# Print Table 2 coverage rates
# -------------------------
cat("\n===== Empirical Coverage Rates (in %) – Table 2 replication =====\n\n")
print(coverage_results)
