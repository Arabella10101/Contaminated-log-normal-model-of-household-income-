# =============================================================================
# Optimization of household income: model comparison
#   (A) cmLN  - Contaminated lognormal, MODE parameterization (your dlnL2)
#   (B) mG    - Gamma, MODE parameterization
#   (C) mIG   - Inverse Gaussian, MODE parameterization
#
# All three models are fit by maximizing the log-likelihood, each written so
# that the mode m is an explicit parameter (rather than mu, shape/rate, or
# mean/shape in their "textbook" forms). This makes the three models directly
# comparable on m, and keeps the multi-start "nudge path" search consistent
# across models.
#
# Per your last note: NO bias/MSE here (real data has no true parameter to
# compare against). We report parameter estimates, logLik, AIC, BIC, and a
# convergence rate (successful multi-starts / total multi-starts) for each
# model, with bar charts for AIC, BIC, and convergence rate.
# =============================================================================

options(scipen = 999) # turn off sci notation

library(readr)
library(dplyr)
library(moments)
library(writexl)
library(ggplot2)

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

x <- household_income$INCOME
n <- length(x)

# =============================================================================
# Shared multi-start "nudge path" search infrastructure
# =============================================================================
n_steps <- 50   # number of multi-starts per model

get_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
start_m <- get_mode(x)

# Generic driver: given a log-likelihood function of par = c(log(m), disp_raw,
# ...extra...) and a sequence of starting values for the dispersion
# parameter(s), run optim() from each start and keep the best converged fit.
# Returns the best fit's parameters (already back-transformed by `postprocess`)
# plus the convergence rate across all starts attempted.
fit_mode_model <- function(x, loglik_fn, m_init, disp_inits, extra_inits = NULL,
                           postprocess) {
  n_tries <- nrow(disp_inits)
  results <- vector("list", n_tries)
  converged_flags <- logical(n_tries)
  logliks <- rep(NA_real_, n_tries)
  
  for (i in seq_len(n_tries)) {
    init_pars <- c(log(m_init), unlist(disp_inits[i, ]), extra_inits)
    est <- try(optim(par = init_pars, fn = loglik_fn, x = x,
                     control = list(fnscale = -1, maxit = 2000)),
               silent = TRUE)
    
    if (inherits(est, "try-error") || !is.finite(est$value) || est$convergence != 0) {
      next
    }
    converged_flags[i] <- TRUE
    logliks[i] <- est$value
    results[[i]] <- est
  }
  
  n_converged <- sum(converged_flags)
  if (n_converged == 0) {
    return(list(par = NA, logLik = NA, n_converged = 0, n_tries = n_tries))
  }
  
  best_idx <- which.max(logliks)
  best_est <- results[[best_idx]]
  out <- postprocess(best_est$par)
  out$logLik <- best_est$value
  out$n_converged <- n_converged
  out$n_tries <- n_tries
  out
}

# =============================================================================
# (A) cmLN - contaminated lognormal, MODE parameterization
#     par = c(log(m), log(sigma), lambda_raw, epsilon_raw)
# =============================================================================
eps_lo <- 0.5
lam_hi <- 50

inv_lam <- function(l) log((l - 1) / (lam_hi - l))
inv_eps <- function(e) log((e - eps_lo) / (1 - e))

dlnL_cmLN <- function(par, x) {
  m       <- exp(par[1])
  sigma   <- exp(par[2])
  mu      <- log(m) + sigma^2
  lambda  <- 1 + (lam_hi - 1) / (1 + exp(-par[3]))
  epsilon <- eps_lo + (1 - eps_lo) / (1 + exp(-par[4]))
  f1 <- dlnorm(x, meanlog = mu, sdlog = sigma)
  f2 <- dlnorm(x, meanlog = mu + (lambda - 1) * sigma^2, sdlog = sqrt(lambda) * sigma)
  sum(log(epsilon * f1 + (1 - epsilon) * f2))
}

start_sigma  <- sd(log(x[x < median(x)]))
base_lambda  <- 1
base_eps     <- 0.99
lambda_step  <- 0.1
eps_step     <- -0.01

cmLN_disp_inits <- do.call(rbind, lapply(seq_len(n_steps), function(i) {
  curr_L <- max(1.001, min(lam_hi - 0.001, base_lambda + i * lambda_step))
  curr_E <- max(eps_lo + 0.001, min(0.99, base_eps + i * eps_step))
  data.frame(log_sigma = log(start_sigma), lambda_raw = inv_lam(curr_L), eps_raw = inv_eps(curr_E))
}))

fit_cmLN <- fit_mode_model(
  x, dlnL_cmLN, m_init = start_m, disp_inits = cmLN_disp_inits,
  postprocess = function(par) {
    list(
      m       = exp(par[1]),
      sigma   = exp(par[2]),
      lambda  = 1 + (lam_hi - 1) / (1 + exp(-par[3])),
      epsilon = eps_lo + (1 - eps_lo) / (1 + exp(-par[4]))
    )
  }
)
k_cmLN <- 4

# =============================================================================
# (B) mG - Gamma, MODE parameterization
#     mode(Gamma) = (shape - 1) / rate,  requires shape > 1
#     par = c(log(m), log(shape - 1))  =>  rate = (shape - 1) / m
# =============================================================================
dlnL_mG <- function(par, x) {
  m     <- exp(par[1])
  shape <- 1 + exp(par[2])
  rate  <- (shape - 1) / m
  sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
}

mG_disp_inits <- data.frame(log_shape_m1 = log(seq(0.3, 15, length.out = n_steps)))

fit_mG <- fit_mode_model(
  x, dlnL_mG, m_init = start_m, disp_inits = mG_disp_inits,
  postprocess = function(par) {
    m     <- exp(par[1])
    shape <- 1 + exp(par[2])
    rate  <- (shape - 1) / m
    list(m = m, shape = shape, rate = rate)
  }
)
k_mG <- 2

# =============================================================================
# (C) mIG - Inverse Gaussian, MODE parameterization
#     Standard IG(mu, lambda) has mean mu, shape lambda; its mode has no closed
#     form in mu, so we solve for mu numerically given (m, lambda):
#       mode = mu * [ sqrt(1 + 9r^2/4) - 3r/2 ],  r = mu / lambda
#     Substituting mu = r*lambda collapses this to a single-variable equation
#     in r for fixed lambda and target mode m:
#       h(r) = r*lambda*( sqrt(1+2.25 r^2) - 1.5 r ) - m = 0
#     As r -> infinity, h(r) -> lambda/3 - m, so feasibility (existence of a
#     root) is guaranteed whenever lambda > 3m. We enforce that by construction:
#       par = c(log(m), disp_raw),  lambda = 3*m*(1 + exp(disp_raw))
#     which makes lambda/3 - m = m*exp(disp_raw) > 0 always.
# =============================================================================
dinvgauss_manual <- function(x, mu, lambda, log = FALSE) {
  logf <- 0.5 * (log(lambda) - log(2 * pi * x^3)) - lambda * (x - mu)^2 / (2 * mu^2 * x)
  if (log) logf else exp(logf)
}

# Solve r from h(r) = r*lambda*(sqrt(1+2.25r^2) - 1.5r) - m = 0, then mu = r*lambda
solve_mu_from_mode <- function(m, lambda) {
  h <- function(r) r * lambda * (sqrt(1 + 2.25 * r^2) - 1.5 * r) - m
  
  upper <- 10
  # h(upper) must turn positive eventually since h(r) -> lambda/3 - m > 0
  tries <- 0
  while (h(upper) <= 0 && tries < 60) {
    upper <- upper * 2
    tries <- tries + 1
  }
  if (h(upper) <= 0) return(NA_real_)  # numerically could not bracket a root
  
  root <- try(uniroot(h, lower = 1e-8, upper = upper)$root, silent = TRUE)
  if (inherits(root, "try-error")) return(NA_real_)
  root * lambda
}

dlnL_mIG <- function(par, x) {
  m      <- exp(par[1])
  lambda <- 3 * m * (1 + exp(par[2]))
  mu     <- solve_mu_from_mode(m, lambda)
  if (is.na(mu) || mu <= 0) return(-1e10)
  ll <- sum(dinvgauss_manual(x, mu, lambda, log = TRUE))
  if (!is.finite(ll)) return(-1e10)
  ll
}

mIG_disp_inits <- data.frame(disp_raw = seq(-2, 3, length.out = n_steps))

fit_mIG <- fit_mode_model(
  x, dlnL_mIG, m_init = start_m, disp_inits = mIG_disp_inits,
  postprocess = function(par) {
    m      <- exp(par[1])
    lambda <- 3 * m * (1 + exp(par[2]))
    mu     <- solve_mu_from_mode(m, lambda)
    list(m = m, mu = mu, lambda = lambda)
  }
)
k_mIG <- 2

# =============================================================================
# Comparison table
# =============================================================================
AIC_BIC <- function(logLik, k, n) {
  c(AIC = 2 * k - 2 * logLik, BIC = k * log(n) - 2 * logLik)
}

ab_cmLN <- AIC_BIC(fit_cmLN$logLik, k_cmLN, n)
ab_mG   <- AIC_BIC(fit_mG$logLik,   k_mG,   n)
ab_mIG  <- AIC_BIC(fit_mIG$logLik,  k_mIG,  n)

comparison_table <- data.frame(
  model       = c("cmLN", "mG", "mIG"),
  k           = c(k_cmLN, k_mG, k_mIG),
  m           = c(fit_cmLN$m, fit_mG$m, fit_mIG$m),
  param2      = c(fit_cmLN$sigma, fit_mG$shape, fit_mIG$lambda),
  param2_name = c("sigma", "shape", "lambda"),
  lambda      = c(fit_cmLN$lambda, NA, NA),
  epsilon     = c(fit_cmLN$epsilon, NA, NA),
  logLik      = c(fit_cmLN$logLik, fit_mG$logLik, fit_mIG$logLik),
  AIC         = c(ab_cmLN["AIC"], ab_mG["AIC"], ab_mIG["AIC"]),
  BIC         = c(ab_cmLN["BIC"], ab_mG["BIC"], ab_mIG["BIC"]),
  n_converged = c(fit_cmLN$n_converged, fit_mG$n_converged, fit_mIG$n_converged),
  n_tries     = c(fit_cmLN$n_tries, fit_mG$n_tries, fit_mIG$n_tries)
)
comparison_table$convergence_rate <- comparison_table$n_converged / comparison_table$n_tries
comparison_table$AIC_rank <- rank(comparison_table$AIC, na.last = "keep")
comparison_table$BIC_rank <- rank(comparison_table$BIC, na.last = "keep")

num_cols <- c("m", "param2", "lambda", "epsilon", "logLik", "AIC", "BIC", "convergence_rate")
comparison_table[, num_cols] <- round(comparison_table[, num_cols], 4)

print(comparison_table)

cat(sprintf("\nAIC prefers: %s\nBIC prefers: %s\n",
            comparison_table$model[which.min(comparison_table$AIC)],
            comparison_table$model[which.min(comparison_table$BIC)]))

# =============================================================================
# Graphs: AIC, BIC, convergence rate
# =============================================================================
plot_dir <- "model_comparison_plots"
dir.create(plot_dir, showWarnings = FALSE)

model_levels <- c("cmLN", "mG", "mIG")
comparison_table$model <- factor(comparison_table$model, levels = model_levels)

aic_plot <- ggplot(comparison_table, aes(x = model, y = AIC, fill = model)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(AIC, 1)), vjust = -0.4) +
  labs(title = "AIC by model", x = NULL, y = "AIC") +
  theme_minimal() +
  theme(legend.position = "none")

bic_plot <- ggplot(comparison_table, aes(x = model, y = BIC, fill = model)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(BIC, 1)), vjust = -0.4) +
  labs(title = "BIC by model", x = NULL, y = "BIC") +
  theme_minimal() +
  theme(legend.position = "none")

conv_plot <- ggplot(comparison_table, aes(x = model, y = convergence_rate * 100, fill = model)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(convergence_rate * 100, 1), "%")), vjust = -0.4) +
  labs(title = "Convergence rate by model",
       subtitle = paste0("Successful multi-starts out of ", n_steps, " attempts"),
       x = NULL, y = "Convergence rate (%)") +
  ylim(0, 105) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "aic_by_model.png"), aic_plot, width = 6, height = 4.5, dpi = 300)
ggsave(file.path(plot_dir, "bic_by_model.png"), bic_plot, width = 6, height = 4.5, dpi = 300)
ggsave(file.path(plot_dir, "convergence_rate_by_model.png"), conv_plot, width = 6, height = 4.5, dpi = 300)

cat("\nPlots saved to:", normalizePath(plot_dir), "\n")

# =============================================================================
# Write comparison table to xlsx
# =============================================================================
out_path <- "household_income_model_comparison.xlsx"
write_xlsx(list(comparison = comparison_table), out_path)
cat("\nModel comparison written to:", out_path, "\n")




