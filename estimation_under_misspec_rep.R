# === Libraries ===
library(fitdistrplus)  # Parametric distribution fitting
library(actuar)        # For Burr distribution

# === Simulation parameters ===
p_quantile <- 0.99          # quantile of interest
n_replications <- 200       # Monte Carlo replications
n_bootstrap_reps <- 200     # bootstrap replicates per sample

# ------------------------------------------------------------------
# Solve convex quadratic program
# 
# Uses subset enumeration (active set) strategy.
# ------------------------------------------------------------------
solve_convex_by_subset <- function(Sigma, tol = 1e-10) {
  k <- ncol(Sigma)
  Sigma <- Sigma + 1e-12 * diag(k)  # tiny ridge for numerical stability
  
  best_val <- Inf
  best_lambda <- rep(0, k)
  
  for (s in 1:k) {
    cmbs <- combn(k, s, simplify = FALSE)  # all subsets of size s
    for (idx in cmbs) {
      SigS <- Sigma[idx, idx, drop = FALSE]
      if (rcond(SigS) < 1e-14) next
      invSigS <- tryCatch(solve(SigS), error = function(e) NULL)
      if (is.null(invSigS)) next
      
      oneS <- rep(1, length(idx))
      lamS <- invSigS %*% oneS
      denom <- as.numeric(t(oneS) %*% lamS)
      if (abs(denom) < .Machine$double.eps) next
      lamS <- as.numeric(lamS / denom)  # normalize weights
      
      if (all(lamS > -tol)) {
        lam <- rep(0, k)
        lam[idx] <- pmax(lamS, 0)
        val <- as.numeric(t(lam) %*% Sigma %*% lam)
        if (val < best_val) {
          best_val <- val
          best_lambda <- lam
        }
      }
    }
  }
  
  # fallback: equal weights if nothing works
  if (sum(best_lambda) <= 0) best_lambda <- rep(1/k, k)
  best_lambda <- best_lambda / sum(best_lambda)
  
  return(list(lambda = best_lambda, obj = best_val))
}

# ------------------------------------------------------------------
# Compute all quantile estimators for one sample
#   - Nonparametric sample quantile
#   - Parametric quantiles: Weibull, Gamma, Burr
# Returns a list of estimates
# ------------------------------------------------------------------
compute_all_estimators <- function(x, p_quantile) {
  out <- list(qNP = NA, qW = NA, qG = NA, qB = NA)
  tiny <- 1e-12
  x_pos <- pmax(x[is.finite(x)], tiny)  # clamp away from 0, drop non-finite
  
  # --- Nonparametric (empirical) quantile ---
  out$qNP <- tryCatch(
    quantile(x_pos, p_quantile, type = 7, names = FALSE),
    error = function(e) NA_real_
  )
  
  # --- Weibull MLE quantile ---
  out$qW <- tryCatch({
    stW <- tryCatch(as.list(mmedist(x_pos, "weibull")$estimate),
                    error = function(e) list(shape = 1.5, scale = median(x_pos)))
    fw <- suppressWarnings(
      fitdist(x_pos, "weibull", method = "mle",
              start = stW,
              lower = c(shape = tiny, scale = tiny),
              optim.method = "L-BFGS-B",
              silent = TRUE)
    )
    qweibull(p_quantile, fw$estimate["shape"], fw$estimate["scale"])
  }, error = function(e) NA_real_)
  
  # --- Gamma MLE quantile ---
  out$qG <- tryCatch({
    stG <- tryCatch(as.list(mmedist(x_pos, "gamma")$estimate),
                    error = function(e) list(shape = 1, rate = 1 / mean(x_pos)))
    fg <- suppressWarnings(
      fitdist(x_pos, "gamma", method = "mle",
              start = stG,
              lower = c(shape = tiny, rate = tiny),
              optim.method = "L-BFGS-B",
              silent = TRUE)
    )
    qgamma(p_quantile, fg$estimate["shape"], fg$estimate["rate"])
  }, error = function(e) NA_real_)
  
  # --- Burr MLE quantile ---
  out$qB <- tryCatch({
    fb <- suppressWarnings(
      fitdist(x_pos, "burr", method = "mle",
              start = list(shape1 = 1, shape2 = 1),
              lower = c(shape1 = tiny, shape2 = tiny),
              optim.method = "L-BFGS-B",
              silent = TRUE)
    )
    qburr(p_quantile, fb$estimate["shape1"], fb$estimate["shape2"])
  }, error = function(e) NA_real_)
  
  return(out)
}

# ------------------------------------------------------------------
# Main estimator function
#   1. Compute original estimators
#   2. Bootstrap to estimate covariance matrix of deviations
#   3. Solve QP for optimal convex weights
#   4. Form averaged estimator
# ------------------------------------------------------------------
get_estimators <- function(x, p_quantile, true_q, n_bootstrap_reps) {
  # Step 1: original estimators
  estimators <- compute_all_estimators(x, p_quantile)
  
  # Step 2: bootstrap replicates
  B <- n_bootstrap_reps
  boots <- matrix(NA, nrow = B, ncol = 4)
  colnames(boots) <- c("qNP", "qW", "qG", "qB")
  for (b in 1:B) {
    xb <- sample(x, size = length(x), replace = TRUE)
    eb <- compute_all_estimators(xb, p_quantile)
    boots[b, ] <- c(eb$qNP, eb$qW, eb$qG, eb$qB)
  }
  boots_clean <- boots[complete.cases(boots), , drop = FALSE]
  
  # Step 3: compute covariance and solve QP
  weights <- rep(NA, 4)
  names(weights) <- c("qNP", "qW", "qG", "qB")
  if (nrow(boots_clean) >= 5) {
    center <- estimators$qNP
    deviations <- sweep(boots_clean, 2, center, "-")
    Sigma_hat <- cov(deviations)
    sol <- tryCatch(solve_convex_by_subset(Sigma_hat), error = function(e) NULL)
    
    if (!is.null(sol)) {
      weights <- sol$lambda
    } else {
      # fallback: inverse variance weights
      diag_m <- diag(Sigma_hat)
      diag_m[diag_m <= 0 | is.na(diag_m)] <- max(diag_m, na.rm = TRUE)
      w <- (1 / diag_m) / sum(1 / diag_m)
      weights <- setNames(as.numeric(w), names(weights))
      warning("Subset solver failed; fell back to inverse-variance weights.")
    }
  } else {
    # fallback: equal weights
    weights <- rep(1/4, 4)
    warning("Too few successful bootstrap fits; using equal weights.")
  }
  
  # Step 4: averaged estimator
  orig_vals <- c(estimators$qNP, estimators$qW, estimators$qG, estimators$qB)
  na_mask <- is.na(orig_vals)
  if (any(na_mask)) {
    weights[na_mask] <- 0
    if (sum(weights) <= 0) {
      weights <- rep(1 / sum(!na_mask), 4) * (!na_mask)
    } else {
      weights <- weights / sum(weights)
    }
  }
  estimators$qAV <- sum(weights * orig_vals, na.rm = TRUE)
  
  # record extra info
  estimators$weights <- weights
  if (exists("Sigma_hat")) estimators$Sigma_hat <- Sigma_hat else estimators$Sigma_hat <- NULL
  estimators$boots_used <- nrow(boots_clean)
  
  return(estimators)
}

# ------------------------------------------------------------------
# Run Monte Carlo simulation for a given distribution and sample size
#   - Generate n samples repeatedly
#   - Compute all estimators
#   - Evaluate MSE across replications
# ------------------------------------------------------------------
run_simulation <- function(dist_name, n, n_replications, p_quantile, n_bootstrap_reps) {
  # Select distribution
  if (dist_name == "Weibull") {
    shape <- 3; scale <- 2
    true_q <- qweibull(p_quantile, shape = shape, scale = scale)
    r_func <- function(n) rweibull(n, shape, scale)
  } else if (dist_name == "Gamma") {
    shape <- 3; rate <- 2
    true_q <- qgamma(p_quantile, shape = shape, rate = rate)
    r_func <- function(n) rgamma(n, shape, rate)
  } else if (dist_name == "Burr") {
    shape1 <- 2; shape2 <- 1
    true_q <- qburr(p_quantile, shape1 = shape1, shape2 = shape2)
    r_func <- function(n) rburr(n, shape1, shape2)
  } else if (dist_name == "Lognormal") {
    meanlog <- 0; sdlog <- 1
    true_q <- qlnorm(p_quantile, meanlog = meanlog, sdlog = sdlog)
    r_func <- function(n) rlnorm(n, meanlog, sdlog)
  }
  
  # Store results for each replication
  results <- matrix(NA, nrow = n_replications, ncol = 5,
                    dimnames = list(NULL, c("qW", "qG", "qB", "qNP", "qAV")))
  
  for (i in 1:n_replications) {
    x <- r_func(n)
    est <- get_estimators(x, p_quantile, true_q, n_bootstrap_reps)
    results[i, ] <- c(est$qW, est$qG, est$qB, est$qNP, est$qAV)
  }
  
  # Monte Carlo MSE and its Monte Carlo SD
  mse_results <- colMeans((results - true_q)^2, na.rm = TRUE)
  sd_results <- apply((results - true_q)^2, 2, sd, na.rm = TRUE) / sqrt(n_replications)
  
  return(list(mse = mse_results, sd = sd_results))
}

# ------------------------------------------------------------------
# Loop over distributions and sample sizes, collect results
# ------------------------------------------------------------------
results_all <- list()
distributions <- c("Weibull", "Gamma", "Burr", "Lognormal")
n_sizes <- c(100, 1000)
estimators <- c("qW", "qG", "qB", "qNP", "qAV")

total_jobs <- length(distributions) * length(n_sizes)
job_idx <- 0

for (dist_name in distributions) {
  for (n in n_sizes) {
    job_idx <- job_idx + 1
    cat(sprintf("[%02d/%02d] Running simulation for %-9s with n = %4d ... ",
                job_idx, total_jobs, dist_name, n))
    flush.console()
    
    sim_result <- run_simulation(dist_name, n, n_replications, p_quantile, n_bootstrap_reps)
    results_all[[paste0(dist_name, "_n", n)]] <- list(
      dist_name = dist_name,
      n = n,
      mse = sim_result$mse,
      sd  = sim_result$sd
    )
    
    cat("done\n")
  }
}
# ------------------------------------------------------------------
# Build wide-format results table
# Rows = distributions, columns = Estimator_n
# ------------------------------------------------------------------
make_colname <- function(est, n) paste0(est, "_n", n)

final_table_data <- matrix("", 
                           nrow = length(distributions), 
                           ncol = length(n_sizes) * length(estimators),
                           dimnames = list(distributions, as.character(
                             unlist(lapply(n_sizes, function(nn) sapply(estimators, make_colname, n = nn)))
                           )))

for (key in names(results_all)) {
  rec <- results_all[[key]]
  dist_name <- rec$dist_name
  n         <- rec$n
  mse_vec   <- rec$mse
  sd_vec    <- rec$sd
  
  for (est in estimators) {
    mse_val <- if (!is.null(mse_vec[est])) mse_vec[est] else NA
    sd_val  <- if (!is.null(sd_vec[est])) sd_vec[est] else NA
    if (dist_name == "Weibull" && !is.na(mse_val)) {
      mse_val <- mse_val * 1000
      sd_val  <- sd_val  * 1000
    }
    formatted_val <- if (is.na(mse_val)) "NA" else sprintf("%.2f (%.2f)", mse_val, sd_val)
    final_table_data[dist_name, paste0(est, "_n", n)] <- formatted_val
  }
}

final_table_df <- as.data.frame(final_table_data, stringsAsFactors = FALSE)

# ------------------------------------------------------------------
# Print results
# ------------------------------------------------------------------
cat("\n=== Simulation summary (wide format: Estimator_n) ===\n")
print(final_table_df, right = FALSE)
cat("\nNote: Weibull MSE/SD values have been multiplied by 1000 for display.\n")
