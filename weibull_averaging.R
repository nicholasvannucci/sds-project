############################################################
# A general procedure to combine estimators - Sec. 4.2 (Weibull)
# R implementation with didactic, line-by-line English comments
############################################################

## ---------------------------------------------------------
## 0) Setup & utilities
## ---------------------------------------------------------

# Reproducibility
set.seed(123)

# (Optional) For nicer plots; if unavailable, base R will be used.
.use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

# Small clipping helper to avoid numerical issues
clip_num <- function(x, lo = -Inf, hi = Inf) pmin(pmax(x, lo), hi)

# Simple progress bar
.make_pb <- function(total, label="") utils::txtProgressBar(min=0, max=total, style=3, width=40)


## ---------------------------------------------------------
## 1) Weibull data generation
## ---------------------------------------------------------
# X ~ Weibull(beta, eta): density f(x) = (beta/eta^beta) x^(beta-1) exp(-(x/eta)^beta), x>0
# In R: rweibull(n, shape=beta, scale=eta)
rweibull_data <- function(n, beta, eta) {
  if (beta <= 0 || eta <= 0) stop("Parametri Weibull non validi.")
  stats::rweibull(n, shape = beta, scale = eta)
}


## ---------------------------------------------------------
## 2) Classical estimators for beta and eta
##    2.1) MM: Method of Moments for beta (+ eta from beta)
##    2.2) OLS: Weibull plot for beta
##    2.3) ML: Maximum Likelihood for beta (+ eta)
## ---------------------------------------------------------

## 2.1) METHOD OF MOMENTS (MM)
# Equation: s^2 / mean^2 = Γ(1+2/b) / Γ(1+1/b)^2 - 1
# Solve for beta with uniroot() on a wide interval, first searching for a sign change.
beta_mm <- function(x, use_unbiased_var = TRUE) {
  x <- x[is.finite(x) & x > 0]
  n <- length(x)
  if (n < 3) return(NA_real_)
  mean_x <- mean(x)
  # Variance: you may choose "unbiased" (R's default var) or "biased" for stability
  var_x <- if (use_unbiased_var) stats::var(x) else mean((x-mean_x)^2)
  if (!is.finite(mean_x) || mean_x <= 0) return(NA_real_)
  r_obs <- var_x / (mean_x^2)
  
  # Root function: f(beta) = Γ(1+2/b)/Γ(1+1/b)^2 - 1 - r_obs
  f <- function(b) {
    if (!is.finite(b) || b <= 0) return(NA_real_)
    a <- 1 + 2/b
    c <- 1 + 1/b
    val <- gamma(a) / (gamma(c)^2) - 1 - r_obs
    if (!is.finite(val)) return(NA_real_)
    val
  }
  
  # Find a bracket with sign change (exploratory grid)
  grid <- c(seq(0.1, 2, length.out = 50), seq(2, 10, length.out = 50))
  fg   <- sapply(grid, f)
  # look for an interval with sign change
  idx <- which(is.finite(fg[-length(fg)]) & is.finite(fg[-1]) & (fg[-length(fg)] * fg[-1] < 0))
  if (length(idx) > 0) {
    a <- grid[idx[1]]; b <- grid[idx[1]+1]
    out <- try(stats::uniroot(f, lower = a, upper = b), silent = TRUE)
    if (!inherits(out, "try-error")) return(out$root)
  }
  
  # fallback: try the whole bracket [0.1, 10]
  out <- try(stats::uniroot(f, lower = 0.1, upper = 10), silent = TRUE)
  if (!inherits(out, "try-error")) return(out$root)
  
  # last fallback: NA
  NA_real_
}

# Eta from MM, given beta_MM: eta_hat = mean / Γ(1+1/beta)
eta_from_beta_mm <- function(x, beta_hat) {
  if (!is.finite(beta_hat) || beta_hat <= 0) return(NA_real_)
  mean(x) / gamma(1 + 1/beta_hat)
}


## 2.2) OLS (Weibull plot) for beta
# Linearization: log(-log(1 - F(x))) = beta*log(x) - beta*log(eta)
# With mean ranks F(x_(i)) ≈ i/(n+1), the OLS slope estimates beta.
beta_ols <- function(x) {
  x <- sort(x[is.finite(x) & x > 0])
  n <- length(x)
  if (n < 3) return(NA_real_)
  u <- (1:n)/(n+1)
  # avoid exactly 0 or 1
  u <- clip_num(u, 1e-10, 1-1e-10)
  Y <- log(-log(1 - u))
  X <- log(x)
  if (!all(is.finite(X)) || !all(is.finite(Y))) return(NA_real_)
  fit <- stats::lm(Y ~ X)
  as.numeric(coef(fit)[["X"]])
}


## 2.3) ML: Maximum likelihood for beta and eta
# Score equation for beta:
# g(beta) = n/beta + sum(log x) - n * (sum x^beta log x) / (sum x^beta) = 0
# Then eta_ML = ( (1/n) sum x^beta )^(1/beta)
mle_beta_eta <- function(x, init_beta = NA_real_) {
  x <- x[is.finite(x) & x > 0]
  n <- length(x)
  if (n < 3) return(c(beta = NA_real_, eta = NA_real_))
  logs <- log(x); s_logx <- sum(logs)
  
  g <- function(b) {
    if (!is.finite(b) || b <= 0) return(NA_real_)
    xb <- x^b
    s_xb <- sum(xb)
    s_xb_logx <- sum(xb * logs)
    n/b + s_logx - n * (s_xb_logx / s_xb)
  }
  
  # Robust initialization: OLS or MM -> then uniroot on a bracket with sign change
  if (!is.finite(init_beta) || init_beta <= 0) {
    init_beta <- beta_ols(x)
    if (!is.finite(init_beta) || init_beta <= 0) {
      init_beta <- beta_mm(x)
      if (!is.finite(init_beta) || init_beta <= 0) init_beta <- 1.0
    }
  }
  
  # Find a bracket where g changes sign
  grid <- c(seq(0.1, 2, length.out = 40), seq(2, 10, length.out = 40))
  gg <- sapply(grid, g)
  idx <- which(is.finite(gg[-length(gg)]) & is.finite(gg[-1]) & (gg[-length(gg)]*gg[-1] < 0))
  beta_hat <- NA_real_
  
  if (length(idx) > 0) {
    a <- grid[idx[1]]; b <- grid[idx[1]+1]
    out <- try(stats::uniroot(g, lower = a, upper = b), silent = TRUE)
    if (!inherits(out, "try-error")) beta_hat <- out$root
  }
  
  # fallback: safeguarded Newton
  if (!is.finite(beta_hat)) {
    b <- init_beta
    for (iter in 1:80) {
      gb <- g(b)
      if (!is.finite(gb)) break
      h  <- 1e-4 * max(1, b)
      gbp <- g(b + h)
      gbn <- g(max(1e-6, b - h))
      dg  <- (gbp - gbn) / ( (if (b - h <= 0) h else 2*h) )
      if (!is.finite(dg) || abs(dg) < 1e-12) {
        # small step in the descent direction
        bnew <- b - sign(gb)*0.1
      } else {
        step <- gb/dg
        bnew <- b - step
      }
      if (bnew <= 0.05) bnew <- 0.5*(b + 0.05)
      if (abs(bnew - b) < 1e-8) { b <- bnew; break }
      b <- bnew
    }
    beta_hat <- b
  }
  
  if (!is.finite(beta_hat) || beta_hat <= 0) return(c(beta = NA_real_, eta = NA_real_))
  
  # Eta ML
  xb <- x^beta_hat
  eta_hat <- (mean(xb))^(1/beta_hat)
  
  c(beta = beta_hat, eta = eta_hat)
}


## ---------------------------------------------------------
## 3) Estimators vector T: (beta_ML, beta_MM, beta_OLS, eta_ML)
## ---------------------------------------------------------
compute_estimators_T <- function(x) {
  # ML (beta, eta)
  ml <- mle_beta_eta(x)
  beta_ML <- ml[["beta"]]; eta_ML <- ml[["eta"]]
  
  # MM (beta only; eta-MM is used only as a check, the paper uses eta_ML)
  beta_MM <- beta_mm(x)
  
  # OLS (beta)
  beta_OLS <- beta_ols(x)
  
  # T vector in fixed order (as in the paper)
  T <- c(beta_ML, beta_MM, beta_OLS, eta_ML)
  names(T) <- c("beta_ML","beta_MM","beta_OLS","eta_ML")
  T
}


## ---------------------------------------------------------
## 4) Σ estimation via parametric bootstrap (as in the paper)
##     - Starting values: beta0 = mean( beta_ML, beta_MM, beta_OLS ), eta0 = eta_ML
##     - J: constraint matrix (4x2) as in Sec. 4.2
##     - Σ ≈ (1/B) * sum_b (T^(b) - J theta0) (T^(b) - J theta0)^T
## ---------------------------------------------------------
estimate_Sigma_param_boot <- function(beta0, eta0, n, B = 1000, ridge = 1e-8, verbose = FALSE) {
  # J matrix (4x2): column 1 for beta, column 2 for eta
  # Constraints: for beta, sum of weights on (beta_ML, beta_MM, beta_OLS) = 1; eta_ML does not contribute
  #              for eta, weight 1 on eta_ML; sum of weights on betas = 0 (column 2).
  J <- matrix(c(
    1, 0,  # beta_ML
    1, 0,  # beta_MM
    1, 0,  # beta_OLS
    0, 1   # eta_ML
  ), nrow = 4, ncol = 2, byrow = TRUE)
  
  theta0 <- c(beta0, eta0)
  JT0    <- as.numeric(J %*% theta0)
  
  # Collect the 4 components of the estimators over B bootstrap replications
  Ts <- matrix(NA_real_, nrow = 4, ncol = B)
  if (verbose) pb <- .make_pb(B, "bootstrap Σ")
  
  for (b in 1:B) {
    xb <- rweibull_data(n, beta = beta0, eta = eta0)
    T  <- compute_estimators_T(xb)
    Ts[, b] <- T
    if (verbose) utils::setTxtProgressBar(pb, b)
  }
  if (verbose) close(pb)
  
  # Center with respect to J theta0
  D <- Ts - JT0
  # Σ_hat = (1/B) * D %*% t(D)
  Sigma_hat <- (D %*% t(D)) / B
  
  # Small ridge for numerical stability
  Sigma_hat <- Sigma_hat + diag(ridge, 4)
  list(Sigma = Sigma_hat, J = J)
}


## ---------------------------------------------------------
## 5) Averaging: paper formula (Sec. 2, Λmax case)
##     theta_AV = (J' Σ^{-1} J)^{-1} J' Σ^{-1} T
## ---------------------------------------------------------
averaging_theta <- function(T, Sigma_hat, J) {
  # T: vector (4), Sigma_hat: 4x4, J: 4x2
  Sinv <- solve(Sigma_hat)
  M    <- solve(t(J) %*% Sinv %*% J) %*% (t(J) %*% Sinv %*% matrix(T, ncol = 1))
  as.numeric(M)  # c(beta_AV, eta_AV)
}


## ---------------------------------------------------------
## 6) One-shot demo on a single dataset (useful for debugging)
## ---------------------------------------------------------

one_shot <- function(beta_true = 2, eta_true = 10, n = 40, B = 300, verbose = TRUE) {
  x  <- rweibull_data(n, beta_true, eta_true)
  T  <- compute_estimators_T(x)
  
  # Bootstrap starting values, as in the paper:
  beta0 <- mean(T[1:3], na.rm = TRUE) # average of (beta_ML, beta_MM, beta_OLS)
  if (!is.finite(beta0) || beta0 <= 0) beta0 <- T[["beta_ML"]]
  eta0  <- T[["eta_ML"]]
  
  est   <- estimate_Sigma_param_boot(beta0 = beta0, eta0 = eta0, n = n, B = B, verbose = verbose)
  th_av <- averaging_theta(T, est$Sigma, est$J)
  
  out <- list(
    true      = c(beta_true, eta_true),
    T         = T,
    theta_AV  = setNames(th_av, c("beta_AV", "eta_AV")),
    Sigma_hat = est$Sigma,
    J         = est$J
  )
  if (verbose) print(out)
  invisible(out)
}


## ---------------------------------------------------------
## 7) Outer Monte Carlo: MSE and (optional) coverage
##     - R: number of outer replications (e.g., 1e4 in the paper)
##     - B: inner bootstrap for Σ (e.g., 1000 in the paper)
## ---------------------------------------------------------
mc_weibull <- function(beta_true, eta_true = 10, n, R = 500, B = 300,
                       compute_coverage = FALSE, verbose = TRUE) {
  # Results container
  RES <- data.frame(
    beta_ML  = numeric(R),
    beta_MM  = numeric(R),
    beta_OLS = numeric(R),
    beta_AV  = numeric(R),
    eta_ML   = numeric(R),
    eta_AV   = numeric(R)
  )
  
  # For coverage (optional): asymptotic variance estimated by (J' Σ^{-1} J)^{-1}
  beta_cover <- eta_cover <- rep(NA, R)
  
  if (verbose) pb <- .make_pb(R, "Monte Carlo")
  
  for (r in 1:R) {
    # 1) Simulate data
    x <- rweibull_data(n, beta_true, eta_true)
    
    # 2) Single-sample estimators
    T <- compute_estimators_T(x)
    
    # 3) Starting values for Σ
    beta0 <- mean(T[1:3], na.rm = TRUE)
    if (!is.finite(beta0) || beta0 <= 0) {
      # robust fallback
      beta0 <- T[["beta_ML"]]
      if (!is.finite(beta0) || beta0 <= 0) beta0 <- 1.0
    }
    eta0  <- T[["eta_ML"]]
    if (!is.finite(eta0) || eta0 <= 0) eta0 <- stats::median(x)
    
    # 4) Σ via parametric bootstrap and averaging
    est   <- estimate_Sigma_param_boot(beta0, eta0, n, B = B, ridge = 1e-8, verbose = FALSE)
    th_av <- averaging_theta(T, est$Sigma, est$J)
    
    # 5) Save results
    RES$beta_ML[r]  <- T[["beta_ML"]]
    RES$beta_MM[r]  <- T[["beta_MM"]]
    RES$beta_OLS[r] <- T[["beta_OLS"]]
    RES$beta_AV[r]  <- th_av[1]
    RES$eta_ML[r]   <- T[["eta_ML"]]
    RES$eta_AV[r]   <- th_av[2]
    
    # 6) (Optional) 95% coverage using asymptotic variance ≈ (J' Σ^{-1} J)^{-1}
    if (compute_coverage) {
      Sinv  <- solve(est$Sigma)
      Vhat  <- solve(t(est$J) %*% Sinv %*% est$J) # 2x2
      se_b  <- sqrt(Vhat[1,1])
      se_e  <- sqrt(Vhat[2,2])
      beta_cover[r] <- as.integer((th_av[1] - 1.96*se_b <= beta_true) && (beta_true <= th_av[1] + 1.96*se_b))
      eta_cover[r]  <- as.integer((th_av[2] - 1.96*se_e <= eta_true) && (eta_true <= th_av[2] + 1.96*se_e))
    }
    
    if (verbose) utils::setTxtProgressBar(pb, r)
  }
  if (verbose) close(pb)
  
  # MSE
  MSE <- c(
    beta_ML  = mean((RES$beta_ML  - beta_true)^2, na.rm = TRUE),
    beta_MM  = mean((RES$beta_MM  - beta_true)^2, na.rm = TRUE),
    beta_OLS = mean((RES$beta_OLS - beta_true)^2, na.rm = TRUE),
    beta_AV  = mean((RES$beta_AV  - beta_true)^2, na.rm = TRUE),
    eta_ML   = mean((RES$eta_ML   - eta_true)^2, na.rm = TRUE),
    eta_AV   = mean((RES$eta_AV   - eta_true)^2, na.rm = TRUE)
  )
  
  out <- list(results = RES, MSE = MSE)
  if (compute_coverage) {
    out$coverage <- c(
      beta_AV = mean(beta_cover, na.rm = TRUE),
      eta_AV  = mean(eta_cover,  na.rm = TRUE)
    )
  }
  out
}

## ---------------------------------------------------------
## 8) Figure 1: distributions of beta estimators (ML, MM, OLS, AV)
##     - Scenario A: beta=0.5, eta=10, n=20
##     - Scenario B: beta=3,   eta=10, n=20
## ---------------------------------------------------------
figure1_like <- function(R = 1000, B = 300, n = 20, eta_true = 10,
                         betas = c(0.5, 3), verbose = TRUE) {
  # For each beta in 'betas' collect beta_ML, beta_MM, beta_OLS, beta_AV over R replications
  PLOTS <- list()
  for (btrue in betas) {
    mc <- mc_weibull(beta_true = btrue, eta_true = eta_true, n = n, R = R, B = B,
                     compute_coverage = FALSE, verbose = verbose)
    DF <- mc$results
    # long format
    long <- data.frame(
      estimator = factor(rep(c("ML","MM","OLS","AV"), each = nrow(DF)),
                         levels = c("ML","MM","OLS","AV")),
      value     = c(DF$beta_ML, DF$beta_MM, DF$beta_OLS, DF$beta_AV),
      beta_true = btrue
    )
    
    if (.use_ggplot) {
      library(ggplot2)
      p <- ggplot2::ggplot(long, ggplot2::aes(x = estimator, y = value)) +
        ggplot2::geom_violin(fill = "grey85", color = "grey40", trim = TRUE) +
        ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.7) +
        ggplot2::geom_hline(yintercept = btrue, linetype = 2, color = "red") +
        ggplot2::labs(title = paste0("Figura 1-like: n=", n, ", eta=10, beta=", btrue),
                      y = expression(hat(beta)), x = NULL) +
        ggplot2::theme_minimal(base_size = 13)
      print(p)
      PLOTS[[as.character(btrue)]] <- p
    } else {
      # Base R fallback
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar), add = TRUE)
      boxplot(value ~ estimator, data = long,
              main = paste0("Figura 1-like (base): n=", n, ", eta=10, beta=", btrue),
              ylab = expression(hat(beta)), col = "lightgrey", border = "grey25")
      abline(h = btrue, lty = 2, col = "red")
    }
  }
  invisible(PLOTS)
}


############################################################
# FAST BLOCK
# Quick simulation params: B = 300, R = 300
############################################################

# --- A) One-shot on a single dataset
out1 <- one_shot(beta_true = 2, eta_true = 10, n = 40, B = 300, verbose = TRUE)

# --- B) Monte Carlo for MSE (small values for a test; then scale to R=1e4, B=1000)
mc_out <- mc_weibull(beta_true = 2, eta_true = 10, n = 20, R = 300, B = 300,
                     compute_coverage = TRUE, verbose = TRUE)
mc_out$MSE
mc_out$coverage  # coverage ~ 95% (as in the paper), with shorter intervals

# --- C) Figure 1-like (reduced)
figure1_like(R = 500, B = 300, n = 20, eta_true = 10, betas = c(0.5, 3), verbose = TRUE)

############################################################
# PAPER BLOCK: reproduces Sec. 4.2 results
# Paper-like parameters: R = 10000, B = 1000
############################################################

cat("\n================= DEMO PAPER =================\n")

# One-shot on a single dataset
out1 <- one_shot(beta_true = 2, eta_true = 10, n = 40, B = 1000, verbose = TRUE)
cat("Stimatori classici (T):\n"); print(out1$T)
cat("Averaging (beta_AV, eta_AV):\n"); print(out1$theta_AV)

# Outer Monte Carlo to estimate MSE and coverage
mc_out <- mc_weibull(beta_true = 2, eta_true = 10, n = 20,
                     R = 10000, B = 1000, compute_coverage = TRUE, verbose = TRUE)
cat("MSE stimatori:\n"); print(mc_out$MSE)
cat("Coverage IC 95% (averaging):\n"); print(mc_out$coverage)

# Figure 1 (beta = 0.5 and beta = 3)
figure1_like(R = 10000, B = 1000, n = 20, eta_true = 10,
             betas = c(0.5, 3), verbose = TRUE)

############################################################
# FINAL NOTES
# - Expectations:
#    MSE(beta_AV) << MSE(beta_ML/MM/OLS) in many scenarios (Table 3 in the paper).
#    MSE(eta_AV) ≈ MSE(eta_ML), i.e., it does not worsen (Table 4).
#    Coverage ~ 95% for beta_AV and eta_AV with asymptotic CIs, and shorter lengths (Table 5).
# - Numerical stability: if solve(Sigma) is unstable, increase 'ridge' in estimate_Sigma_param_boot.
# - Performance: for large R and B=1000, consider parallelization (future.apply / parallel).
############################################################
