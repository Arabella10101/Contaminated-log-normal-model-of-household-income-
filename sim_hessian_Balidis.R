#sim using a Hessian check to select params
set.seed(123)

dlnL2 <- function(par, x){
  m     <- exp(par[1]) 
  sigma <- exp(par[2])
  lambda <- exp(par[3]) + 1 #bound lambda > 1
  # Logit transformation to keep epsilon in (0, 1)
  epsilon <- 1 / (1 + exp(-par[4])) 
  
  # Guard against numerical underflow/overflow during optim's line search
  if (!is.finite(m) || !is.finite(sigma) || sigma <= 0 ||
      !is.finite(lambda) || lambda <= 1 || !is.finite(epsilon)){
    return(-1e10)
  }
  
  # Check the actual meanlog/sdlog values that will be passed to dlnorm.
  # lambda*sigma^2 and sqrt(lambda)*sigma can overflow to Inf even when
  # lambda and sigma are each individually finite, so check the combined
  # values directly rather than trusting the inputs alone.
  meanlog1 <- log(m) + sigma^2
  sdlog1   <- sigma
  meanlog2 <- log(m) + lambda * sigma^2
  sdlog2   <- sqrt(lambda) * sigma
  
  if (!is.finite(meanlog1) || !is.finite(sdlog1) || sdlog1 <= 0 ||
      !is.finite(meanlog2) || !is.finite(sdlog2) || sdlog2 <= 0){
    return(-1e10)
  }
  
  dens <- epsilon*dlnorm(x, meanlog=meanlog1, sdlog=sdlog1) +
    (1-epsilon)*dlnorm(x, meanlog=meanlog2, sdlog=sdlog2)
  
  if (any(!is.finite(dens)) || any(dens <= 0)) return(-1e10)
  
  return(sum(log(dens)))
}

#fixed true values
m     <- 2
sigma <- 0.5

#scenarios
n_vals       <- c(100, 1000, 10000)
lambda_vals  <- c(2, 5)
epsilon_vals <- c(0.6, 0.9)
n_reps       <- 100
n_starts     <- 5

# ---------------------------------------------------------------------------
# Data-generating function: fixed, deterministic split by epsilon.
# ---------------------------------------------------------------------------
rcontlnorm <- function(n, m, sigma, lambda, epsilon){
  n1 <- round(epsilon * n)
  n2 <- n - n1
  
  x1 <- rlnorm(n = n1, meanlog = log(m) + sigma^2,        sdlog = sigma)
  x2 <- rlnorm(n = n2, meanlog = log(m) + lambda*sigma^2, sdlog = sqrt(lambda)*sigma)
  
  c(x1, x2)
}

# ---------------------------------------------------------------------------
# Random starting values for optim, on the transformed (unconstrained) scale
# used by dlnL2.
# ---------------------------------------------------------------------------
get_start <- function(x, jitter_sd = 0.5){
  m_guess     <- mean(x)
  sigma_guess <- sd(log(x))
  c(
    log(m_guess)                 + rnorm(1, 0, jitter_sd),  # -> m
    log(max(sigma_guess, 0.05))  + rnorm(1, 0, jitter_sd),  # -> sigma
    log(1)                       + rnorm(1, 0, jitter_sd),  # -> lambda - 1
    rnorm(1, 0, 1)                                          # -> epsilon (logit scale)
  )
}

# ---------------------------------------------------------------------------
# Check that a converged solution is a genuine local maximum, not a saddle
# point or flat region optim stopped at prematurely.
#
# optim() is called with fnscale = -1, so internally it MINIMIZES -dlnL2(par).
# The Hessian it returns is therefore the Hessian of -loglik at the solution.
# At a true maximum of the log-likelihood, -loglik is locally convex there,
# so that Hessian should be positive definite (all eigenvalues > 0).
# ---------------------------------------------------------------------------
is_valid_maximum <- function(hess, tol = 1e-8){
  if (is.null(hess) || any(!is.finite(hess))) return(FALSE)
  eig <- eigen(hess, symmetric = TRUE, only.values = TRUE)$values
  all(is.finite(eig)) && all(eig > tol)
}

# ---------------------------------------------------------------------------
# Multi-start optimizer: call optim() on the full likelihood from several
# random starts, tracking the best result found so far. A start only becomes
# the new "best" if it (a) has a higher log-likelihood AND (b) passes the
# Hessian check above — i.e. it's a genuine local maximum, not a spurious
# stopping point.
# ---------------------------------------------------------------------------
fit_multistart <- function(x, n_starts){
  best_val <- -Inf
  best_par <- rep(NA_real_, 4)
  
  for (s in seq_len(n_starts)){
    par0 <- get_start(x)
    
    fit <- tryCatch(
      optim(par0, dlnL2, x = x, method = "BFGS",
            control = list(fnscale = -1, maxit = 1000), hessian = TRUE),
      error = function(e) NULL
    )
    
    if (!is.null(fit) && fit$convergence == 0 && is.finite(fit$value) &&
        fit$value > best_val && is_valid_maximum(fit$hessian)){
      best_val <- fit$value
      best_par <- fit$par
    }
  }
  
  if (!is.finite(best_val)){
    return(list(par = rep(NA_real_, 4), loglik = NA_real_, convergence = 1L))
  }
  list(par = best_par, loglik = best_val, convergence = 0L)
}

# ---------------------------------------------------------------------------
# Main simulation loop
# ---------------------------------------------------------------------------
scenarios <- expand.grid(n = n_vals, lambda_true = lambda_vals, epsilon_true = epsilon_vals,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

results_list <- vector("list", nrow(scenarios) * n_reps)
row_i <- 1L
total_scenarios <- nrow(scenarios)

for (sc in seq_len(total_scenarios)){
  n_i       <- scenarios$n[sc]
  lambda_i  <- scenarios$lambda_true[sc]
  epsilon_i <- scenarios$epsilon_true[sc]
  
  cat(sprintf("Scenario %d/%d: n=%d, lambda=%.1f, epsilon=%.1f\n",
              sc, total_scenarios, n_i, lambda_i, epsilon_i))
  
  for (r in seq_len(n_reps)){
    
    x <- rcontlnorm(n_i, m = m, sigma = sigma, lambda = lambda_i, epsilon = epsilon_i)
    fit <- fit_multistart(x, n_starts = n_starts)
    
    if (fit$convergence == 0L){
      m_hat       <- exp(fit$par[1])
      sigma_hat   <- exp(fit$par[2])
      lambda_hat  <- exp(fit$par[3]) + 1
      epsilon_hat <- 1 / (1 + exp(-fit$par[4]))
    } else {
      m_hat <- sigma_hat <- lambda_hat <- epsilon_hat <- NA_real_
    }
    
    results_list[[row_i]] <- data.frame(
      n            = n_i,
      lambda_true  = lambda_i,
      epsilon_true = epsilon_i,
      rep          = r,
      m_hat        = m_hat,
      sigma_hat    = sigma_hat,
      lambda_hat   = lambda_hat,
      epsilon_hat  = epsilon_hat,
      loglik       = fit$loglik,
      convergence  = fit$convergence
    )
    row_i <- row_i + 1L
  }
}

results <- do.call(rbind, results_list)

# ---------------------------------------------------------------------------
# Summarize: bias, sd, and RMSE for each parameter, by scenario
# ---------------------------------------------------------------------------
summarize_scenario <- function(df){
  ok <- df[df$convergence == 0, ]
  
  best_row <- ok[which.max(ok$loglik), ]   # the single replicate with the highest log-likelihood
  
  data.frame(
    n_converged = nrow(ok),
    n_total     = nrow(df),
    
    m_best      = best_row$m_hat,
    sigma_best  = best_row$sigma_hat,
    lambda_best = best_row$lambda_hat,
    epsilon_best= best_row$epsilon_hat,
    loglik_best = best_row$loglik,
    
    m_mean  = mean(ok$m_hat),
    m_bias  = mean(ok$m_hat) - m,
    m_sd    = sd(ok$m_hat),
    m_mse   = mean((ok$m_hat - m)^2),
    
    sigma_mean = mean(ok$sigma_hat),
    sigma_bias = mean(ok$sigma_hat) - sigma,
    sigma_sd   = sd(ok$sigma_hat),
    sigma_mse  = mean((ok$sigma_hat - sigma)^2),
    
    lambda_mean = mean(ok$lambda_hat),
    lambda_bias = mean(ok$lambda_hat) - df$lambda_true[1],
    lambda_sd   = sd(ok$lambda_hat),
    lambda_mse  = mean((ok$lambda_hat - df$lambda_true[1])^2),
    
    epsilon_mean = mean(ok$epsilon_hat),
    epsilon_bias = mean(ok$epsilon_hat) - df$epsilon_true[1],
    epsilon_sd   = sd(ok$epsilon_hat),
    epsilon_mse  = mean((ok$epsilon_hat - df$epsilon_true[1])^2)
  )
}

split_keys <- interaction(results$n, results$lambda_true, results$epsilon_true, drop = TRUE)
summary_list <- lapply(split(results, split_keys), summarize_scenario)

summary_df <- do.call(rbind, lapply(names(summary_list), function(k){
  sub <- results[split_keys == k, ][1, c("n", "lambda_true", "epsilon_true")]
  cbind(sub, summary_list[[k]])
}))
rownames(summary_df) <- NULL
summary_df <- summary_df[order(summary_df$n, summary_df$lambda_true, summary_df$epsilon_true), ]

print(summary_df)

saveRDS(results,    "contaminated_lognormal_results_raw.rds")
saveRDS(summary_df, "contaminated_lognormal_results_summary.rds")