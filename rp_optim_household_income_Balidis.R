#optim of household income
options(scipen = 999) #turn off sci notation

library(readr)
library(dplyr)
library(moments)

households_data <- read_csv("C:/Users/arabe/Documents/Research_Project/Fact_IES2023_Households.csv")

household_income <- households_data %>% 
  select(INCOME)

summary(household_income)

desc_stats <- data.frame(
  n        = length(household_income$INCOME),
  Mean     = mean(household_income$INCOME),
  SD       = sd(household_income$INCOME),
  Skewness = skewness(household_income$INCOME),
  Kurtosis = kurtosis(household_income$INCOME),   # raw/Pearson kurtosis
  Min      = min(household_income$INCOME),
  Max      = max(household_income$INCOME)
)

print(desc_stats)

library(writexl)

# =============================================================================
# Three competing models fit to REAL household income data (no true params
# exist for real data, so per your call: report parameter estimates + AIC/BIC
# only -- no bias/MSE here).
#   (A) Single lognormal                      -- closed-form MLE
#   (B) Contaminated lognormal, REGULAR (mu) parameterization
#   (C) Contaminated lognormal, MODE (m) parameterization -- your dlnL2,
#       kept as the multi-start "nudge path" you already had, just wrapped
#       into a reusable function so (B) and (C) are fit identically.
# =============================================================================

# Bounded so lambda/epsilon can't run off toward the degenerate region
# (same fix used in the simulation study -- keeps both mixture fits stable).
eps_lo <- 0.5
lam_hi <- 50

inv_lam <- function(l) log((l - 1) / (lam_hi - l))
inv_eps <- function(e) log((e - eps_lo) / (1 - e))

# ---- (B) regular (mu) parameterization ----
dlnL_mix_reg <- function(par, x) {
  mu      <- par[1]
  sigma   <- exp(par[2])
  lambda  <- 1 + (lam_hi - 1) / (1 + exp(-par[3]))
  epsilon <- eps_lo + (1 - eps_lo) / (1 + exp(-par[4]))
  f1 <- dlnorm(x, meanlog = mu, sdlog = sigma)
  f2 <- dlnorm(x, meanlog = mu + (lambda - 1) * sigma^2, sdlog = sqrt(lambda) * sigma)
  sum(log(epsilon * f1 + (1 - epsilon) * f2))
}

# ---- (C) mode (m) parameterization -- your dlnL2, same bounded lambda/epsilon
# transform as (B) so the two are directly comparable and equally stable.
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
# Multi-start "nudge path" search -- your original loop, generalized so it
# works for either parameterization.
# =============================================================================
fit_mixture_path <- function(x, loglik_fn, par1_init, par2_init,
                             n_steps = 50, base_lambda = 1, base_eps = 0.99,
                             lambda_step = 0.1, eps_step = -0.01) {
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

# =============================================================================
# Starting values (your originals)
# =============================================================================
get_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
start_m     <- get_mode(household_income$INCOME)
start_sigma <- sd(log(household_income$INCOME[household_income$INCOME < median(household_income$INCOME)]))
start_lambda <- 2
start_eps    <- 0.9

x <- household_income$INCOME
n <- length(x)

k_single <- 2
k_mix    <- 4

n_steps     <- 50
base_lambda <- 1
base_eps    <- 0.99
lambda_step <- 0.1
eps_step    <- -0.01

# =============================================================================
# (A) Single lognormal -- closed-form MLE
# =============================================================================
mu_s     <- mean(log(x))
sigma_s  <- sqrt(mean((log(x) - mu_s)^2))
logLik_s <- sum(dlnorm(x, meanlog = mu_s, sdlog = sigma_s, log = TRUE))
AIC_s    <- 2 * k_single - 2 * logLik_s
BIC_s    <- k_single * log(n) - 2 * logLik_s
m_s      <- exp(mu_s - sigma_s^2)   # mode of the fitted single lognormal, for comparability

cat(sprintf("(A) Single lognormal:      m=%.3f sigma=%.3f logLik=%.2f AIC=%.2f BIC=%.2f\n",
            m_s, sigma_s, logLik_s, AIC_s, BIC_s))

# =============================================================================
# (B) Contaminated lognormal, regular (mu) parameterization
# =============================================================================
fit_reg <- fit_mixture_path(x, dlnL_mix_reg,
                            par1_init = mu_s, par2_init = log(sigma_s),
                            n_steps, base_lambda, base_eps, lambda_step, eps_step)
m_reg    <- if (is.na(fit_reg$par1)) NA else exp(fit_reg$par1 - fit_reg$sigma^2)
AIC_reg  <- if (is.na(fit_reg$logLik)) NA else 2 * k_mix - 2 * fit_reg$logLik
BIC_reg  <- if (is.na(fit_reg$logLik)) NA else k_mix * log(n) - 2 * fit_reg$logLik

cat(sprintf("(B) Mixture (regular):     m=%.3f sigma=%.3f lambda=%.3f epsilon=%.3f logLik=%.2f AIC=%.2f BIC=%.2f (n_converged=%d/%d)\n",
            m_reg, fit_reg$sigma, fit_reg$lambda, fit_reg$epsilon, fit_reg$logLik, AIC_reg, BIC_reg,
            fit_reg$n_converged, n_steps))

# =============================================================================
# (C) Contaminated lognormal, mode (m) parameterization -- your dlnL2
# =============================================================================
fit_mode <- fit_mixture_path(x, dlnL_mix_mode,
                             par1_init = log(start_m), par2_init = log(start_sigma),
                             n_steps, base_lambda, base_eps, lambda_step, eps_step)
m_mode   <- if (is.na(fit_mode$par1)) NA else exp(fit_mode$par1)
AIC_mode <- if (is.na(fit_mode$logLik)) NA else 2 * k_mix - 2 * fit_mode$logLik
BIC_mode <- if (is.na(fit_mode$logLik)) NA else k_mix * log(n) - 2 * fit_mode$logLik

cat(sprintf("(C) Mixture (mode):        m=%.3f sigma=%.3f lambda=%.3f epsilon=%.3f logLik=%.2f AIC=%.2f BIC=%.2f (n_converged=%d/%d)\n",
            m_mode, fit_mode$sigma, fit_mode$lambda, fit_mode$epsilon, fit_mode$logLik, AIC_mode, BIC_mode,
            fit_mode$n_converged, n_steps))

# =============================================================================
# Comparison table -- estimates + AIC/BIC for all three models, no bias/MSE
# (no true parameter exists for real data).
# =============================================================================
comparison_table <- data.frame(
  model      = c("Single lognormal", "Mixture (regular)", "Mixture (mode)"),
  k          = c(k_single, k_mix, k_mix),
  m          = c(m_s, m_reg, m_mode),
  sigma_sq   = c(sigma_s^2, fit_reg$sigma^2, fit_mode$sigma^2),
  lambda     = c(NA, fit_reg$lambda, fit_mode$lambda),
  epsilon    = c(NA, fit_reg$epsilon, fit_mode$epsilon),
  logLik     = c(logLik_s, fit_reg$logLik, fit_mode$logLik),
  AIC        = c(AIC_s, AIC_reg, AIC_mode),
  BIC        = c(BIC_s, BIC_reg, BIC_mode),
  n_converged_of_starts = c(NA, fit_reg$n_converged, fit_mode$n_converged)
)
comparison_table$AIC_rank <- rank(comparison_table$AIC, na.last = "keep")
comparison_table$BIC_rank <- rank(comparison_table$BIC, na.last = "keep")
comparison_table[, c("m", "sigma_sq", "lambda", "epsilon", "logLik", "AIC", "BIC")] <-
  round(comparison_table[, c("m", "sigma_sq", "lambda", "epsilon", "logLik", "AIC", "BIC")], 4)
print(comparison_table)

cat(sprintf("\nAIC prefers: %s\nBIC prefers: %s\n",
            comparison_table$model[which.min(comparison_table$AIC)],
            comparison_table$model[which.min(comparison_table$BIC)]))

# =============================================================================
# Write to a SEPARATE xlsx from the simulation study output
# =============================================================================
out_path <- "household_income_model_comparison.xlsx"
write_xlsx(list(comparison = comparison_table), out_path)
cat("\nModel comparison written to:", out_path, "\n")










