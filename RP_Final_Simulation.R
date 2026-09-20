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
# Model log-likelihood
# =============================================================================

# Contaminated mode parameterized lognormal,
# par = c(log(m), log(sigma), logit(lambda), logit(epsilon)). mu is recovered
# via mu = log(m) + sigma^2 before being plugged in.
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

results      <- list()
summary_rows <- list()

for (n in n_vals) {
  for (lambda in lambda_vals) {
    for (epsilon in epsilon_vals) {

      n1 <- round(n * epsilon)
      n2 <- n - n1

      rep_results <- matrix(NA, nrow = n_reps, ncol = 6)
      colnames(rep_results) <- c(
        "m_hat", "sigma_hat", "lambda_hat", "epsilon_hat", "logLik", "n_converged"
      )

      for (r in seq_len(n_reps)) {
        x1 <- rlnorm(n1, meanlog = log(m) + sigma^2, sdlog = sigma)
        x2 <- rlnorm(n2, meanlog = log(m) + lambda*sigma^2, sdlog = sqrt(lambda)*sigma)
        x  <- c(x1, x2)

        mu0    <- mean(log(x))
        sigma0 <- sd(log(x))
        m0     <- exp(mu0 - sigma0^2)   # mode-space starting value

        # ---- contaminated lognormal, mode parameterization ----
        fit_mode <- fit_mixture_path(x, dlnL_mix_mode, par1_init = log(m0), par2_init = log(sigma0),
                                     n_steps, base_lambda, base_eps, lambda_step, eps_step)
        m_hat_mode <- if (is.na(fit_mode$par1)) NA else exp(fit_mode$par1)

        rep_results[r, ] <- c(
          m_hat_mode, fit_mode$sigma, fit_mode$lambda, fit_mode$epsilon,
          fit_mode$logLik, fit_mode$n_converged
        )
      }

      key <- sprintf("n=%d_lambda=%g_epsilon=%g", n, lambda, epsilon)
      results[[key]] <- rep_results
      cat("Done:", key, "\n")

      saveRDS(results, "results_checkpoint.rds")
      saveRDS(summary_rows, "summary_rows_checkpoint.rds")

      best_rep_idx <- which.max(rep_results[, "logLik"])
      cat(sprintf("  best mixture (mode) fit for this scenario (replicate %d/%d): m=%.3f sigma=%.3f lambda=%.3f epsilon=%.3f logLik=%.2f\n",
                  best_rep_idx, n_reps,
                  rep_results[best_rep_idx, "m_hat"], rep_results[best_rep_idx, "sigma_hat"],
                  rep_results[best_rep_idx, "lambda_hat"], rep_results[best_rep_idx, "epsilon_hat"],
                  rep_results[best_rep_idx, "logLik"]))

      m_hat       <- rep_results[, "m_hat"]
      sigma_hat   <- rep_results[, "sigma_hat"]
      lambda_hat  <- rep_results[, "lambda_hat"]
      epsilon_hat <- rep_results[, "epsilon_hat"]

      # ---- convergence rate ----
      # conv_rate: proportion of the n_reps replicates that produced a usable
      # (finite) fit at all.
      # mean_prop_steps: on average, what fraction of the 15 multi-start
      # initializations converged per replicate (a stability diagnostic).
      conv_rate       <- mean(is.finite(rep_results[, "logLik"]))
      mean_prop_steps <- mean(rep_results[, "n_converged"] / n_steps, na.rm = TRUE)

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

        conv_rate       = conv_rate,
        mean_prop_steps = mean_prop_steps
      )

      cat(sprintf("  -> convergence: %.1f%%\n\n", 100 * conv_rate))
    }
  }
}
summary_df <- do.call(rbind, summary_rows)
print(summary_df)

# =============================================================================
# Post-hoc: report sigma^2 instead of sigma. Recomputed from the raw
# per-replicate sigma_hat already stored in `results` -- NOT a re-run of the
# simulation, and NOT just squaring the existing summary stats (mean(sigma^2)
# != mean(sigma)^2, same issue for bias/MSE).
# =============================================================================
sigma2_true <- sigma^2   # sigma <- 0.5 above, so sigma2_true = 0.25

for (key in names(results)) {
  sigma_hat  <- results[[key]][, "sigma_hat"]
  sigma2_hat <- sigma_hat^2

  summary_rows[[key]]$sigma_mean <- mean(sigma2_hat, na.rm = TRUE)
  summary_rows[[key]]$sigma_sd   <- sd(sigma2_hat, na.rm = TRUE)
  summary_rows[[key]]$sigma_bias <- mean(sigma2_hat, na.rm = TRUE) - sigma2_true
  summary_rows[[key]]$sigma_mse  <- mean((sigma2_hat - sigma2_true)^2, na.rm = TRUE)
}

summary_df <- do.call(rbind, summary_rows)

# =============================================================================
# FINAL SUMMARY -- only the final chosen values per scenario, written to xlsx
# for use in the appendix. Reports recovered estimates, bias, MSE, and
# convergence rate for the mode-parameterized contaminated lognormal.
# =============================================================================
final_table <- summary_df[, c("n", "lambda_true", "epsilon_true",
                              "m_mean", "m_bias", "m_mse",
                              "sigma_mean", "sigma_bias", "sigma_mse",
                              "lambda_mean", "lambda_bias", "lambda_mse",
                              "epsilon_mean", "epsilon_bias", "epsilon_mse",
                              "conv_rate")]
final_table[, 4:16] <- round(final_table[, 4:16], 4)

out_path <- "contaminated_lognormal_mode_final_summary_simulation.xlsx"
write_xlsx(list(summary = final_table), out_path)
cat("\nFinal summary (sheet 'summary') written to:", out_path, "\n")

# =============================================================================
# PLOTS -- reads the final summary directly from the xlsx file, renders
# EVERYTHING in black & white only (no color anywhere). Groups are
# distinguished by linetype + point shape instead of color. Backgrounds
# are white.
#
# NOTE: the "sigma" columns in the xlsx now hold sigma^2 (recomputed
# post-hoc from the raw per-replicate fits -- see the sim script above).
# Labels below say "sigma^2" accordingly; no column-name changes needed.
# =============================================================================

library(readxl)
library(ggplot2)

xlsx_path <- "contaminated_lognormal_mode_final_summary_simulation.xlsx"

summary_df <- read_excel(xlsx_path, sheet = "summary")

dir.create("sim_plots", showWarnings = FALSE)

# ---- white background, black-and-white theme ----
bw_theme <- theme_minimal() + theme(
  panel.background  = element_rect(fill = "white", colour = NA),
  plot.background   = element_rect(fill = "white", colour = NA),
  legend.background = element_rect(fill = "white", colour = NA),
  strip.background  = element_rect(fill = "white", colour = "black"),
  panel.grid.major  = element_line(colour = "grey85"),
  panel.grid.minor  = element_line(colour = "grey92"),
  panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.3)
)

lambda_f <- function(df) factor(paste0("lambda=", df$lambda_true))
epsilon_f <- function(df) factor(paste0("epsilon=", df$epsilon_true))
scenario_label <- function(df) paste0("lambda=", df$lambda_true, ", eps=", df$epsilon_true)

# ---- 1. Convergence rate ----
conv_df <- data.frame(n = summary_df$n, lambda_true = summary_df$lambda_true,
                      epsilon_true = summary_df$epsilon_true, conv_rate = summary_df$conv_rate)
conv_df$lambda_f  <- lambda_f(conv_df)
conv_df$epsilon_f <- epsilon_f(conv_df)

p_conv <- ggplot(conv_df, aes(x = factor(n), y = conv_rate, group = 1)) +
  geom_line(color = "black") + geom_point(color = "black", size = 2.2, shape = 16) +
  facet_grid(lambda_f ~ epsilon_f) +
  labs(x = "n", y = "Proportion of replicates with a valid fit") +
  ylim(0, 1) +
  bw_theme
ggsave("sim_plots/convergence_rate.png", p_conv, width = 8, height = 6, dpi = 150, bg = "white")

# ---- 2. Bias by parameter (4 scenario lines, black & white) ----
bias_df <- rbind(
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "m",         bias = summary_df$m_bias),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "sigma^2",   bias = summary_df$sigma_bias),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "lambda",    bias = summary_df$lambda_bias),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "epsilon",   bias = summary_df$epsilon_bias)
)

scenario_levels <- unique(bias_df$scenario)
linetypes_4 <- c("solid", "dashed", "dotted", "dotdash")
shapes_4    <- c(16, 17, 15, 21)
names(linetypes_4) <- scenario_levels
names(shapes_4)    <- scenario_levels

p_bias <- ggplot(bias_df, aes(x = factor(n), y = bias,
                              linetype = scenario, shape = scenario, group = scenario)) +
  geom_line(color = "black") + geom_point(color = "black", fill = "white", size = 2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.3) +
  scale_linetype_manual(values = linetypes_4) +
  scale_shape_manual(values = shapes_4) +
  facet_wrap(~ parameter, scales = "free_y") +
  labs( x = "n", y = "Bias",
       linetype = "Scenario", shape = "Scenario") +
  bw_theme
ggsave("sim_plots/bias_by_parameter.png", p_bias, width = 9, height = 6, dpi = 150, bg = "white")

# ---- 3. MSE by parameter (log scale, black & white) ----
mse_df <- rbind(
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "m",        mse = summary_df$m_mse),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "sigma^2",  mse = summary_df$sigma_mse),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "lambda",   mse = summary_df$lambda_mse),
  data.frame(n = summary_df$n, scenario = scenario_label(summary_df), parameter = "epsilon",  mse = summary_df$epsilon_mse)
)

p_mse <- ggplot(mse_df, aes(x = factor(n), y = mse,
                            linetype = scenario, shape = scenario, group = scenario)) +
  geom_line(color = "black") + geom_point(color = "black", fill = "white", size = 2) +
  scale_y_log10() +
  scale_linetype_manual(values = linetypes_4) +
  scale_shape_manual(values = shapes_4) +
  facet_wrap(~ parameter, scales = "free_y") +
  labs( x = "n", y = "MSE (log10 scale)",
       linetype = "Scenario", shape = "Scenario") +
  bw_theme
ggsave("sim_plots/mse_by_parameter.png", p_mse, width = 9, height = 6, dpi = 150, bg = "white")

cat("\nBlack-and-white plots (white background) saved to the 'sim_plots/' directory:\n",
    " - convergence_rate.png\n",
    " - bias_by_parameter.png\n",
    " - mse_by_parameter.png\n")
