library(moments)
library(writexl)
library(ggplot2)
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
# will always share one logLik/AIC/BIC value. No optim() needed: the MLE is exact.

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
# two-component mixture.
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
      
      rep_results <- matrix(NA, nrow = n_reps, ncol = 19)
      colnames(rep_results) <- c(
        "m_hat", "sigma_hat", "lambda_hat", "epsilon_hat", "logLik", "n_converged",
        "logLik_single", "AIC_single", "BIC_single",
        "AIC_mix_reg", "BIC_mix_reg",
        "m_hat_mode", "sigma_hat_mode", "lambda_hat_mode", "epsilon_hat_mode",
        "logLik_mix_mode", "AIC_mix_mode", "BIC_mix_mode", "n_converged_mode"
      )
      
      # per-replicate winners under each criterion -- kept as character
      # vectors (not part of the numeric matrix) so we can tally how often
      # each model is picked, instead of reporting an averaged AIC/BIC.
      aic_winner_vec <- character(n_reps)
      bic_winner_vec <- character(n_reps)
      
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
        BIC_s    <- k_single * log(n) - 2 * logLik_s
        
        # ---- (C) contaminated lognormal, regular (mu) parameterization ----
        fit_reg <- fit_mixture_path(x, dlnL_mix_reg, par1_init = mu0, par2_init = log(sigma0),
                                    n_steps, base_lambda, base_eps, lambda_step, eps_step)
        m_hat_reg   <- if (is.na(fit_reg$par1)) NA else exp(fit_reg$par1 - fit_reg$sigma^2)
        AIC_mix_reg <- if (is.na(fit_reg$logLik)) NA else 2 * k_mix - 2 * fit_reg$logLik
        BIC_mix_reg <- if (is.na(fit_reg$logLik)) NA else k_mix * log(n) - 2 * fit_reg$logLik
        
        # ---- (D) contaminated lognormal, mode parameterization ----
        fit_mode <- fit_mixture_path(x, dlnL_mix_mode, par1_init = log(m0), par2_init = log(sigma0),
                                     n_steps, base_lambda, base_eps, lambda_step, eps_step)
        m_hat_mode   <- if (is.na(fit_mode$par1)) NA else exp(fit_mode$par1)
        AIC_mix_mode <- if (is.na(fit_mode$logLik)) NA else 2 * k_mix - 2 * fit_mode$logLik
        BIC_mix_mode <- if (is.na(fit_mode$logLik)) NA else k_mix * log(n) - 2 * fit_mode$logLik
        
        rep_results[r, ] <- c(
          m_hat_reg, fit_reg$sigma, fit_reg$lambda, fit_reg$epsilon, fit_reg$logLik, fit_reg$n_converged,
          logLik_s, AIC_s, BIC_s,
          AIC_mix_reg, BIC_mix_reg,
          m_hat_mode, fit_mode$sigma, fit_mode$lambda, fit_mode$epsilon,
          fit_mode$logLik, AIC_mix_mode, BIC_mix_mode, fit_mode$n_converged
        )
        
        # ---- which model does AIC / BIC prefer for this replicate? ----
        # (only compared when all three models produced a finite value)
        aic_vals <- c(single = AIC_s, mix_reg = AIC_mix_reg, mix_mode = AIC_mix_mode)
        bic_vals <- c(single = BIC_s, mix_reg = BIC_mix_reg, mix_mode = BIC_mix_mode)
        
        aic_winner_vec[r] <- if (all(is.finite(aic_vals))) names(aic_vals)[which.min(aic_vals)] else NA_character_
        bic_winner_vec[r] <- if (all(is.finite(bic_vals))) names(bic_vals)[which.min(bic_vals)] else NA_character_
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
      
      # ---- convergence rates ----
      # conv_rate_*: proportion of the n_reps replicates that produced a
      # usable (finite) fit at all for that parameterization.
      # mean_prop_steps_*: on average, what fraction of the 15 multi-start
      # initializations converged per replicate (a stability diagnostic).
      conv_rate_reg  <- mean(is.finite(rep_results[, "logLik"]))
      conv_rate_mode <- mean(is.finite(rep_results[, "logLik_mix_mode"]))
      mean_prop_steps_reg  <- mean(rep_results[, "n_converged"]      / n_steps, na.rm = TRUE)
      mean_prop_steps_mode <- mean(rep_results[, "n_converged_mode"] / n_steps, na.rm = TRUE)
      
      # ---- AIC / BIC selection rates -- NOT averaged AIC/BIC values ----
      valid_aic   <- !is.na(aic_winner_vec)
      valid_bic   <- !is.na(bic_winner_vec)
      n_valid_aic <- sum(valid_aic)
      n_valid_bic <- sum(valid_bic)
      
      aic_win_single   <- if (n_valid_aic > 0) sum(aic_winner_vec[valid_aic] == "single")   / n_valid_aic else NA
      aic_win_mix_reg  <- if (n_valid_aic > 0) sum(aic_winner_vec[valid_aic] == "mix_reg")  / n_valid_aic else NA
      aic_win_mix_mode <- if (n_valid_aic > 0) sum(aic_winner_vec[valid_aic] == "mix_mode") / n_valid_aic else NA
      
      bic_win_single   <- if (n_valid_bic > 0) sum(bic_winner_vec[valid_bic] == "single")   / n_valid_bic else NA
      bic_win_mix_reg  <- if (n_valid_bic > 0) sum(bic_winner_vec[valid_bic] == "mix_reg")  / n_valid_bic else NA
      bic_win_mix_mode <- if (n_valid_bic > 0) sum(bic_winner_vec[valid_bic] == "mix_mode") / n_valid_bic else NA
      
      # ---- diagnostic: how far apart are the two mixture parameterizations'
      # maximized log-likelihoods, replicate by replicate? Since reg and mode
      # are the SAME likelihood surface, any gap is purely optimizer error.
      # gap > 0 means reg found a higher logLik than mode (reg "won" that rep);
      # gap < 0 means mode found the higher one.
      gap_small_thresh <- 0.01   # nats -- plausible optimizer tolerance noise
      gap_large_thresh <- 1      # nats -- a whole log-likelihood unit; suggests
      # a genuinely worse search, not just noise
      
      logLik_gap <- rep_results[, "logLik"] - rep_results[, "logLik_mix_mode"]
      valid_gap  <- is.finite(logLik_gap)
      n_valid_gap <- sum(valid_gap)
      
      mean_logLik_gap     <- if (n_valid_gap > 0) mean(logLik_gap[valid_gap]) else NA
      min_logLik_gap      <- if (n_valid_gap > 0) min(logLik_gap[valid_gap])  else NA
      max_logLik_gap       <- if (n_valid_gap > 0) max(logLik_gap[valid_gap]) else NA
      mean_abs_logLik_gap  <- if (n_valid_gap > 0) mean(abs(logLik_gap[valid_gap])) else NA
      prop_gap_gt_small    <- if (n_valid_gap > 0) mean(abs(logLik_gap[valid_gap]) > gap_small_thresh) else NA
      prop_gap_gt_large    <- if (n_valid_gap > 0) mean(abs(logLik_gap[valid_gap]) > gap_large_thresh) else NA
      n_reg_better_reps    <- if (n_valid_gap > 0) sum(logLik_gap[valid_gap] >  0) else NA
      n_mode_better_reps   <- if (n_valid_gap > 0) sum(logLik_gap[valid_gap] <  0) else NA
      
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
        
        conv_rate_reg        = conv_rate_reg,
        conv_rate_mode       = conv_rate_mode,
        mean_prop_steps_reg  = mean_prop_steps_reg,
        mean_prop_steps_mode = mean_prop_steps_mode,
        
        n_valid_aic      = n_valid_aic,
        aic_win_single   = aic_win_single,
        aic_win_mix_reg  = aic_win_mix_reg,
        aic_win_mix_mode = aic_win_mix_mode,
        
        n_valid_bic      = n_valid_bic,
        bic_win_single   = bic_win_single,
        bic_win_mix_reg  = bic_win_mix_reg,
        bic_win_mix_mode = bic_win_mix_mode,
        
        n_valid_gap          = n_valid_gap,
        mean_logLik_gap      = mean_logLik_gap,
        min_logLik_gap       = min_logLik_gap,
        max_logLik_gap       = max_logLik_gap,
        mean_abs_logLik_gap  = mean_abs_logLik_gap,
        prop_gap_gt_small    = prop_gap_gt_small,
        prop_gap_gt_large    = prop_gap_gt_large,
        n_reg_better_reps    = n_reg_better_reps,
        n_mode_better_reps   = n_mode_better_reps
      )
      
      cat(sprintf("  -> convergence: reg=%.1f%%  mode=%.1f%%  |  AIC picks mix_mode %.1f%% of reps  |  BIC picks mix_mode %.1f%% of reps\n",
                  100 * conv_rate_reg, 100 * conv_rate_mode,
                  100 * aic_win_mix_mode, 100 * bic_win_mix_mode))
      cat(sprintf("     logLik gap (reg - mode): mean=%.5f  min=%.5f  max=%.5f  |  reg better in %d/%d reps, mode better in %d/%d  |  |gap|>%.2g in %.1f%%, |gap|>%g in %.1f%%\n\n",
                  mean_logLik_gap, min_logLik_gap, max_logLik_gap,
                  n_reg_better_reps, n_valid_gap, n_mode_better_reps, n_valid_gap,
                  gap_small_thresh, 100 * prop_gap_gt_small,
                  gap_large_thresh, 100 * prop_gap_gt_large))
    }
  }
}

summary_df <- do.call(rbind, summary_rows)
print(summary_df)

# =============================================================================
# FINAL SUMMARY -- only the final chosen values per scenario, written to xlsx
# for use in the appendix. Reports recovered estimates, bias, MSE, convergence
# rates, and how often AIC/BIC selected each model (NOT averaged AIC/BIC).
# =============================================================================
final_table <- summary_df[, c("n", "lambda_true", "epsilon_true",
                              "m_mean", "m_bias", "m_mse",
                              "sigma_mean", "sigma_bias", "sigma_mse",
                              "lambda_mean", "lambda_bias", "lambda_mse",
                              "epsilon_mean", "epsilon_bias", "epsilon_mse",
                              "conv_rate_reg", "conv_rate_mode",
                              "aic_win_single", "aic_win_mix_reg", "aic_win_mix_mode",
                              "bic_win_single", "bic_win_mix_reg", "bic_win_mix_mode")]
final_table[, 4:23] <- round(final_table[, 4:23], 4)

# ---- diagnostics table (reg vs. mode optimizer gap) -- separate sheet,
# not part of the appendix table, since it's a numerical-method check rather
# than a substantive result.
diagnostics_table <- summary_df[, c("n", "lambda_true", "epsilon_true",
                                    "n_valid_gap", "mean_logLik_gap", "min_logLik_gap", "max_logLik_gap",
                                    "mean_abs_logLik_gap", "prop_gap_gt_small", "prop_gap_gt_large",
                                    "n_reg_better_reps", "n_mode_better_reps")]
diagnostics_table[, 5:10] <- round(diagnostics_table[, 5:10], 5)

out_path <- "contaminated_lognormal_final_summary.xlsx"
write_xlsx(list(summary = final_table, diagnostics = diagnostics_table), out_path)
cat("\nFinal summary (sheet 'summary') and reg-vs-mode gap diagnostics (sheet 'diagnostics') written to:", out_path, "\n")

# =============================================================================
# PLOTS -- visual counterparts to the tables above, saved as PNGs for the main
# text. The full numeric tables still go into the xlsx above for the appendix.
# =============================================================================
dir.create("plots", showWarnings = FALSE)

lambda_f <- function(df) factor(paste0("lambda=", df$lambda_true))
epsilon_f <- function(df) factor(paste0("epsilon=", df$epsilon_true))
scenario_label <- function(df) paste0("lambda=", df$lambda_true, ", eps=", df$epsilon_true)

# ---- 1. Convergence rate ----
conv_df <- rbind(
  data.frame(n = summary_df$n, lambda_true = summary_df$lambda_true, epsilon_true = summary_df$epsilon_true,
             parameterization = "Regular", conv_rate = summary_df$conv_rate_reg),
  data.frame(n = summary_df$n, lambda_true = summary_df$lambda_true, epsilon_true = summary_df$epsilon_true,
             parameterization = "Mode", conv_rate = summary_df$conv_rate_mode)
)
conv_df$lambda_f  <- lambda_f(conv_df)
conv_df$epsilon_f <- epsilon_f(conv_df)

p_conv <- ggplot(conv_df, aes(x = factor(n), y = conv_rate, color = parameterization, group = parameterization)) +
  geom_line() + geom_point(size = 2) +
  facet_grid(lambda_f ~ epsilon_f) +
  labs(title = "Mixture-model convergence rate by sample size",
       x = "n", y = "Proportion of replicates with a valid fit", color = "Parameterization") +
  ylim(0, 1) +
  theme_minimal()
ggsave("plots/convergence_rate.png", p_conv, width = 8, height = 6, dpi = 150)

# ---- 2. Bias by parameter ----
bias_df <- rbind(
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "m",       bias = summary_df$m_bias),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "sigma",   bias = summary_df$sigma_bias),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "lambda",  bias = summary_df$lambda_bias),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "epsilon", bias = summary_df$epsilon_bias)
)

p_bias <- ggplot(bias_df, aes(x = factor(n), y = bias, color = scenario, group = scenario)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(title = "Bias of recovered parameters by sample size", x = "n", y = "Bias", color = "Scenario") +
  theme_minimal()
ggsave("plots/bias_by_parameter.png", p_bias, width = 9, height = 6, dpi = 150)

# ---- 3. MSE by parameter (log scale -- MSE shrinks by orders of magnitude) ----
mse_df <- rbind(
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "m",       mse = summary_df$m_mse),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "sigma",   mse = summary_df$sigma_mse),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "lambda",  mse = summary_df$lambda_mse),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "epsilon", mse = summary_df$epsilon_mse)
)

p_mse <- ggplot(mse_df, aes(x = factor(n), y = mse, color = scenario, group = scenario)) +
  geom_line() + geom_point() +
  scale_y_log10() +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(title = "MSE of recovered parameters by sample size (log scale)", x = "n", y = "MSE (log10 scale)", color = "Scenario") +
  theme_minimal()
ggsave("plots/mse_by_parameter.png", p_mse, width = 9, height = 6, dpi = 150)

# ---- 4. AIC selection rate ----
aic_win_df <- rbind(
  data.frame(n = summary_df$n, lambda_true = summary_df$lambda_true, epsilon_true = summary_df$epsilon_true,
             model = "Single",         rate = summary_df$aic_win_single),
  data.frame(n = summary_df$n, lambda_true = summary_df$lambda_true, epsilon_true = summary_df$epsilon_true,
             model = "Mixture (reg)",  rate = summary_df$aic_win_mix_reg),
  data.frame(n = summary_df$n, lambda_true = summary_df$lambda_true, epsilon_true = summary_df$epsilon_true,
             model = "Mixture (mode)", rate = summary_df$aic_win_mix_mode)
)
aic_win_df$model     <- factor(aic_win_df$model, levels = c("Single", "Mixture (reg)", "Mixture (mode)"))
aic_win_df$lambda_f  <- lambda_f(aic_win_df)
aic_win_df$epsilon_f <- epsilon_f(aic_win_df)

p_aic <- ggplot(aic_win_df, aes(x = factor(n), y = rate, fill = model)) +
  geom_col(position = "stack") +
  facet_grid(lambda_f ~ epsilon_f) +
  labs(title = "AIC model-selection rate across replicates",
       x = "n", y = "Proportion of replicates selecting each model", fill = "Selected model") +
  theme_minimal()
ggsave("plots/aic_selection_rate.png", p_aic, width = 8, height = 6, dpi = 150)

# ---- 5. BIC selection rate ----
bic_win_df <- rbind(
  data.frame(n = summary_df$n, lambda_true = summary_df$lambda_true, epsilon_true = summary_df$epsilon_true,
             model = "Single",         rate = summary_df$bic_win_single),
  data.frame(n = summary_df$n, lambda_true = summary_df$lambda_true, epsilon_true = summary_df$epsilon_true,
             model = "Mixture (reg)",  rate = summary_df$bic_win_mix_reg),
  data.frame(n = summary_df$n, lambda_true = summary_df$lambda_true, epsilon_true = summary_df$epsilon_true,
             model = "Mixture (mode)", rate = summary_df$bic_win_mix_mode)
)
bic_win_df$model     <- factor(bic_win_df$model, levels = c("Single", "Mixture (reg)", "Mixture (mode)"))
bic_win_df$lambda_f  <- lambda_f(bic_win_df)
bic_win_df$epsilon_f <- epsilon_f(bic_win_df)

p_bic <- ggplot(bic_win_df, aes(x = factor(n), y = rate, fill = model)) +
  geom_col(position = "stack") +
  facet_grid(lambda_f ~ epsilon_f) +
  labs(title = "BIC model-selection rate across replicates",
       x = "n", y = "Proportion of replicates selecting each model", fill = "Selected model") +
  theme_minimal()
ggsave("plots/bic_selection_rate.png", p_bic, width = 8, height = 6, dpi = 150)

cat("\nPlots saved to the 'plots/' directory:\n",
    " - convergence_rate.png\n",
    " - bias_by_parameter.png\n",
    " - mse_by_parameter.png\n",
    " - aic_selection_rate.png\n",
    " - bic_selection_rate.png\n")
