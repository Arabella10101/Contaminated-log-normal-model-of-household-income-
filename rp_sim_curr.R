library(moments)
library(writexl)
set.seed(123)

# Bounded so lambda/epsilon can't run off toward the degenerate region
# (lambda->Inf, epsilon->0) that was corrupting fits in the hardest cells.
eps_lo <- 0.5
lam_hi <- 50

inv_lam <- function(l) log((l - 1) / (lam_hi - l))
inv_eps <- function(e) log((e - eps_lo) / (1 - e))

# =============================================================================
# Model log-likelihoods
# =============================================================================

# (A/B) Single lognormal -- "regular" (mu,sigma) and "mode-parameterized"
# (m,sigma) are the SAME distribution under a bijective reparameterization
# (m = exp(mu - sigma^2)), so they are fit once in closed form below and
# will always share one AIC value. No optim() needed: the MLE is exact.

# (C) Contaminated (mixture) lognormal, REGULAR parameterization.
# par = c(mu, log(sigma), logit(lambda), logit(epsilon)). This is your
# original dlnL2, renamed for clarity.
dlnL_mix_reg <- function(par, x) {
  mu      <- par[1]
  sigma   <- exp(par[2])
  lambda  <- 1 + (lam_hi - 1) / (1 + exp(-par[3]))
  epsilon <- eps_lo + (1 - eps_lo) / (1 + exp(-par[4]))
  f1 <- dlnorm(x, meanlog = mu, sdlog = sigma)
  f2 <- dlnorm(x, meanlog = mu + (lambda - 1) * sigma^2, sdlog = sqrt(lambda) * sigma)
  sum(log(epsilon * f1 + (1 - epsilon) * f2))
}

# (D) Contaminated (mixture) lognormal, MODE parameterization.
# Same mixture density as (C) -- the baseline component's location is just
# searched over in mode-space (par[1] = log(m)) instead of mu-space, with
# mu recovered via mu = log(m) + sigma^2 before being plugged into the same
# two-component mixture. Mathematically the surface has the same global
# maximum as (C); in finite samples the optimizer can land on a different
# value here because it walks a differently-shaped/scaled surface.
dlnL_mix_mode <- function(par, x) {
  m       <- exp(par[1])
  sigma   <- exp(par[2])
  mu      <- log(m) + sigma^2
  lambda  <- 1 + (lam_hi - 1) / (1 + exp(-par[3]))
  epsilon <- eps_lo + (1 - eps_lo) / (1 + exp(-par[4]))
  f1 <- dlnorm(x, meanlog = mu, sdlog = sigma)
  f2 <- dlnorm(x, meanlog = mu + (lambda - 1) * sigma^2, sdlog = sqrt(lambda) * sigma)
  sum(log(epsilon * f1 + (1 - epsilon) * f2))
}

# =============================================================================
# Generic multi-start path search for a mixture log-likelihood.
# Identical algorithm to your original loop; par1_init is in whatever space
# the supplied loglik_fn expects for its first parameter (mu or log(m)).
# =============================================================================
fit_mixture_path <- function(x, loglik_fn, par1_init, par2_init,
                             n_steps = 15, base_lambda = 1, base_eps = 0.99,
                             lambda_step = 0.5, eps_step = -0.032) {
  path <- matrix(NA, nrow = n_steps, ncol = 5)
  colnames(path) <- c("par1", "sigma", "lambda", "epsilon", "logLik")
  
  for (i in 1:n_steps) {
    curr_L <- base_lambda + i * lambda_step
    curr_E <- base_eps    + i * eps_step
    curr_L <- max(1.001, min(lam_hi - 0.001, curr_L))
    curr_E <- max(eps_lo + 0.001, min(0.99, curr_E))
    
    init_pars <- c(par1_init, par2_init, inv_lam(curr_L), inv_eps(curr_E))
    est <- try(optim(par = init_pars, fn = loglik_fn, x = x, control = list(fnscale = -1)), silent = TRUE)
    
    if (inherits(est, "try-error") || !is.finite(est$value) || est$convergence != 0) next
    
    sigma_i   <- exp(est$par[2])
    lambda_i  <- 1 + (lam_hi - 1) / (1 + exp(-est$par[3]))
    epsilon_i <- eps_lo + (1 - eps_lo) / (1 + exp(-est$par[4]))
    path[i, ] <- c(est$par[1], sigma_i, lambda_i, epsilon_i, est$value)
  }
  
  ok <- is.finite(path[, "logLik"])
  if (any(ok)) {
    best_idx <- which(ok)[which.max(path[ok, "logLik"])]
    list(par1 = path[best_idx, "par1"], sigma = path[best_idx, "sigma"],
         lambda = path[best_idx, "lambda"], epsilon = path[best_idx, "epsilon"],
         logLik = path[best_idx, "logLik"], n_converged = sum(ok))
  } else {
    list(par1 = NA, sigma = NA, lambda = NA, epsilon = NA, logLik = NA, n_converged = 0)
  }
}

m     <- 2
sigma <- 0.5

n_vals       <- c(100, 1000, 10000)
lambda_vals  <- c(2, 5)
epsilon_vals <- c(0.6, 0.9)
n_reps       <- 500

n_steps     <- 15
base_lambda <- 1
base_eps    <- 0.99
lambda_step <- 0.5
eps_step    <- -0.032

k_single <- 2
k_mix    <- 4

results      <- list()
summary_rows <- list()

for (n in n_vals) {
  for (lambda in lambda_vals) {
    for (epsilon in epsilon_vals) {
      
      n1 <- round(n * epsilon)
      n2 <- n - n1
      
      rep_results <- matrix(NA, nrow = n_reps, ncol = 16)
      colnames(rep_results) <- c(
        "m_hat", "sigma_hat", "lambda_hat", "epsilon_hat", "logLik", "n_converged",
        "logLik_single", "AIC_single",
        "AIC_mix_reg",
        "m_hat_mode", "sigma_hat_mode", "lambda_hat_mode", "epsilon_hat_mode",
        "logLik_mix_mode", "AIC_mix_mode", "n_converged_mode"
      )
      
      for (r in seq_len(n_reps)) {
        x1 <- rlnorm(n1, meanlog = log(m) + sigma^2, sdlog = sigma)
        x2 <- rlnorm(n2, meanlog = log(m) + lambda*sigma^2, sdlog = sqrt(lambda)*sigma)
        x  <- c(x1, x2)
        
        mu0    <- mean(log(x))
        sigma0 <- sd(log(x))
        m0     <- exp(mu0 - sigma0^2)   # mode-space starting value
        
        # ---- (A/B) single lognormal, closed-form MLE ----
        mu_s     <- mean(log(x))
        sigma_s  <- sqrt(mean((log(x) - mu_s)^2))
        logLik_s <- sum(dlnorm(x, meanlog = mu_s, sdlog = sigma_s, log = TRUE))
        AIC_s    <- 2 * k_single - 2 * logLik_s
        
        # ---- (C) contaminated lognormal, regular (mu) parameterization ----
        fit_reg <- fit_mixture_path(x, dlnL_mix_reg, par1_init = mu0, par2_init = log(sigma0),
                                    n_steps, base_lambda, base_eps, lambda_step, eps_step)
        m_hat_reg   <- if (is.na(fit_reg$par1)) NA else exp(fit_reg$par1 - fit_reg$sigma^2)
        AIC_mix_reg <- if (is.na(fit_reg$logLik)) NA else 2 * k_mix - 2 * fit_reg$logLik
        
        # ---- (D) contaminated lognormal, mode parameterization ----
        fit_mode <- fit_mixture_path(x, dlnL_mix_mode, par1_init = log(m0), par2_init = log(sigma0),
                                     n_steps, base_lambda, base_eps, lambda_step, eps_step)
        m_hat_mode   <- if (is.na(fit_mode$par1)) NA else exp(fit_mode$par1)
        AIC_mix_mode <- if (is.na(fit_mode$logLik)) NA else 2 * k_mix - 2 * fit_mode$logLik
        
        rep_results[r, ] <- c(
          m_hat_reg, fit_reg$sigma, fit_reg$lambda, fit_reg$epsilon, fit_reg$logLik, fit_reg$n_converged,
          logLik_s, AIC_s,
          AIC_mix_reg,
          m_hat_mode, fit_mode$sigma, fit_mode$lambda, fit_mode$epsilon,
          fit_mode$logLik, AIC_mix_mode, fit_mode$n_converged
        )
      }
      
      key <- sprintf("n=%d_lambda=%g_epsilon=%g", n, lambda, epsilon)
      results[[key]] <- rep_results
      cat("Done:", key, "\n")
      
      best_rep_idx <- which.max(rep_results[, "logLik"])
      cat(sprintf("  best mixture (reg) fit for this scenario (replicate %d/%d): m=%.3f sigma=%.3f lambda=%.3f epsilon=%.3f logLik=%.2f\n",
                  best_rep_idx, n_reps,
                  rep_results[best_rep_idx, "m_hat"], rep_results[best_rep_idx, "sigma_hat"],
                  rep_results[best_rep_idx, "lambda_hat"], rep_results[best_rep_idx, "epsilon_hat"],
                  rep_results[best_rep_idx, "logLik"]))
      
      m_hat       <- rep_results[, "m_hat"]
      sigma_hat   <- rep_results[, "sigma_hat"]
      lambda_hat  <- rep_results[, "lambda_hat"]
      epsilon_hat <- rep_results[, "epsilon_hat"]
      
      AIC_single_mean  <- mean(rep_results[, "AIC_single"], na.rm = TRUE)
      AIC_mix_reg_mean <- mean(rep_results[, "AIC_mix_reg"], na.rm = TRUE)
      AIC_mix_mode_mean<- mean(rep_results[, "AIC_mix_mode"], na.rm = TRUE)
      
      aic_vec   <- c(single = AIC_single_mean, mix_reg = AIC_mix_reg_mean, mix_mode = AIC_mix_mode_mean)
      best_model <- names(aic_vec)[which.min(aic_vec)]
      
      summary_rows[[key]] <- data.frame(
        scenario      = key,
        n             = n,
        lambda_true   = lambda,
        epsilon_true  = epsilon,
        
        m_mean        = mean(m_hat, na.rm=TRUE),
        sigma_mean    = mean(sigma_hat, na.rm=TRUE),
        lambda_mean   = mean(lambda_hat, na.rm=TRUE),
        epsilon_mean  = mean(epsilon_hat, na.rm=TRUE),
        
        m_sd          = sd(m_hat, na.rm=TRUE),
        sigma_sd      = sd(sigma_hat, na.rm=TRUE),
        lambda_sd     = sd(lambda_hat, na.rm=TRUE),
        epsilon_sd    = sd(epsilon_hat, na.rm=TRUE),
        
        m_bias        = mean(m_hat, na.rm=TRUE)       - m,
        sigma_bias    = mean(sigma_hat, na.rm=TRUE)    - sigma,
        lambda_bias   = mean(lambda_hat, na.rm=TRUE)   - lambda,
        epsilon_bias  = mean(epsilon_hat, na.rm=TRUE)  - epsilon,
        
        m_mse         = mean((m_hat - m)^2, na.rm=TRUE),
        sigma_mse     = mean((sigma_hat - sigma)^2, na.rm=TRUE),
        lambda_mse    = mean((lambda_hat - lambda)^2, na.rm=TRUE),
        epsilon_mse   = mean((epsilon_hat - epsilon)^2, na.rm=TRUE),
        
        mean_n_converged = mean(rep_results[,"n_converged"], na.rm=TRUE),
        
        AIC_single   = AIC_single_mean,
        AIC_mix_reg  = AIC_mix_reg_mean,
        AIC_mix_mode = AIC_mix_mode_mean,
        best_model   = best_model
      )
      
      cat(sprintf("  -> AIC (mean over reps): single=%.2f  mix_reg=%.2f  mix_mode=%.2f  (best: %s)\n\n",
                  AIC_single_mean, AIC_mix_reg_mean, AIC_mix_mode_mean, best_model))
    }
  }
}

summary_df <- do.call(rbind, summary_rows)
print(summary_df)

# =============================================================================
# FINAL SUMMARY -- only the final chosen values per scenario, written to xlsx
# (no per-replicate/per-path-search-step values, just the scenario-level
# recovered means, bias, MSE, and mean AICs used to pick the best model)
# =============================================================================
final_table <- summary_df[, c("n", "lambda_true", "epsilon_true",
                              "m_mean", "m_bias", "m_mse",
                              "sigma_mean", "sigma_bias", "sigma_mse",
                              "lambda_mean", "lambda_bias", "lambda_mse",
                              "epsilon_mean", "epsilon_bias", "epsilon_mse",
                              "AIC_single", "AIC_mix_reg", "AIC_mix_mode", "best_model")]
final_table[, 4:18] <- round(final_table[, 4:18], 4)

out_path <- "contaminated_lognormal_final_summary.xlsx"
write_xlsx(final_table, out_path)
cat("\nFinal summary written to:", out_path, "\n")
