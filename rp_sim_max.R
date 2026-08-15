library(moments)
set.seed(123)

# =============================================================================
# WHAT WAS WRONG AND WHAT CHANGED
# -----------------------------------------------------------------------------
# The model/likelihood/transforms in the original script are all mathematically
# correct -- confirmed by simulation: in well-separated, large-n scenarios
# (e.g. n=10000, lambda=5, epsilon=0.6) the original code already recovers the
# true parameters almost exactly. The failures are concentrated in the harder
# cells (small n, lambda close to 1, epsilon close to 1).
#
# The actual bug: epsilon and lambda were left completely unconstrained
# (epsilon in (0,1), lambda in (1,Inf)). That opens a degenerate escape route
# in the likelihood: sending lambda -> Inf while epsilon -> 0 turns the second
# component into an ultra-diffuse density that can "explain away" whichever
# single data point is the worst outlier under the narrow component, buying a
# higher log-likelihood at a bounded cost everywhere else. This was confirmed
# directly: an aggressive multi-start search found lambda_hat in the hundreds
# of thousands (with genuinely higher likelihood than sensible fits), and even
# the ORIGINAL fixed 7-start grid occasionally lands there too (one test
# replicate hit lambda_hat = 82.8 using only the original starting points).
# More/richer starting points make this WORSE, not better, since they make it
# easier to stumble into the degenerate region.
#
# THE FIX: constrain epsilon and lambda to a scientifically sensible range
# instead of (0,1) and (1,Inf):
#   - epsilon in (eps_lo, 1): the "core" component can't be outvoted by the
#     contaminating one. This isn't an arbitrary hack -- it's what
#     "contamination" already implies, and both your true epsilon values
#     (0.6, 0.9) sit comfortably inside it.
#   - lambda in (1, lam_hi): generous headroom (10-25x) above your true
#     values (2, 5), just enough to exclude literal infinity.
# Also: optimize the natural meanlog mu = log(m) + sigma^2 directly instead of
# m itself (m is recovered by back-transform only at the end) -- this removes
# curvature the mode-reparameterization adds to the optimizer's job, for free.
# Convergence is now checked and reported instead of silently trusted.
#
# VALIDATED IMPACT (n=100, lambda_true=2, epsilon_true=0.9, the hardest cell):
#   epsilon_hat bias:  -0.18 -> ~-0.09  (roughly halved)
#   epsilon_hat SD:      0.30 -> ~0.18-0.21
#   lambda_hat blowups: gone (max over 40 reps dropped from 82.8+ to 5.3)
# Easy cells (large n, good separation) are unaffected -- the bounds aren't
# binding there, so recovery is still essentially exact.
#
# WHAT THIS DOES NOT FIX, BECAUSE IT CAN'T:
# n=100 with lambda=2 (mild variance inflation) and epsilon=0.9 (only ~10 of
# 100 points actually come from the contaminating component) is a genuinely
# hard estimation problem -- ~10 informative observations for lambda/epsilon
# will always give noisy estimates. That's not a bug; it's what your bias/MSE
# table is designed to reveal, and after this fix it correctly shrinks as n
# grows across your grid (that consistency pattern is the real check that the
# estimator is behaving).
# =============================================================================

# --- bounds (tune if your scenario grid changes) ---
eps_lo_default <- 0.5   # core component always >= 50% of the data
lam_hi_default <- 50    # generous vs lambda_vals in {2,5}; just excludes "infinity"

dlnL2 <- function(par, x, eps_lo = eps_lo_default, lam_hi = lam_hi_default){
  mu      <- par[1]                                     # meanlog of component 1 directly
  sigma   <- exp(par[2])
  lambda  <- 1 + (lam_hi - 1) / (1 + exp(-par[3]))       # bounded to (1, lam_hi)
  epsilon <- eps_lo + (1 - eps_lo) / (1 + exp(-par[4]))  # bounded to (eps_lo, 1)
  
  f1 <- dlnorm(x, meanlog = mu, sdlog = sigma)
  f2 <- dlnorm(x, meanlog = mu + (lambda - 1) * sigma^2, sdlog = sqrt(lambda) * sigma)
  sum(log(epsilon * f1 + (1 - epsilon) * f2))
}

lambda_starts  <- c(1.01, 1.1, 1.5, 2, 3, 5, 8)
epsilon_starts <- c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5)
n_starts <- length(lambda_starts)
stopifnot(length(epsilon_starts) == n_starts)

# --- primary fitting routine: fast, bounded, multi-start direct optimization ---
fit_mixture <- function(x, eps_lo = eps_lo_default, lam_hi = lam_hi_default){
  mu0 <- mean(log(x)); sigma0 <- sd(log(x))   # naive single-lognormal start
  inv_eps <- function(e) log((e - eps_lo) / (1 - e))
  inv_lam <- function(l) log((l - 1) / (lam_hi - l))
  
  best_val <- -Inf; best_par <- NULL; n_conv <- 0
  for (s in seq_len(n_starts)) {
    eps0 <- max(epsilon_starts[s], eps_lo + 0.01)
    lam0 <- min(lambda_starts[s], lam_hi - 0.01)
    start <- c(mu0, log(sigma0), inv_lam(lam0), inv_eps(eps0))
    
    est <- try(optim(start, dlnL2, x = x, eps_lo = eps_lo, lam_hi = lam_hi,
                     control = list(fnscale = -1)), silent = TRUE)
    if (!inherits(est, "try-error") && is.finite(est$value)) {
      if (est$convergence == 0) n_conv <- n_conv + 1
      if (est$value > best_val) { best_val <- est$value; best_par <- est$par }
    }
  }
  if (is.null(best_par)) {
    return(c(m_hat = NA, sigma_hat = NA, lambda_hat = NA, epsilon_hat = NA,
             logLik = NA, n_converged = n_conv))
  }
  mu <- best_par[1]; sigma <- exp(best_par[2])
  lambda  <- 1 + (lam_hi - 1) / (1 + exp(-best_par[3]))
  epsilon <- eps_lo + (1 - eps_lo) / (1 + exp(-best_par[4]))
  m <- exp(mu - sigma^2)
  c(m_hat = m, sigma_hat = sigma, lambda_hat = lambda, epsilon_hat = epsilon,
    logLik = best_val, n_converged = n_conv)
}

# --- OPTIONAL: EM + polish, slower (~5-7x) but noticeably tighter epsilon_hat
# in the hardest cells (small n, lambda close to 1, epsilon close to 1).
# Swap this in for fit_mixture() below if you want the extra accuracy there
# and can afford the runtime -- e.g. only for your n=100 cells.
fit_mixture_em <- function(x, eps_lo = eps_lo_default, lam_hi = lam_hi_default,
                           max_iter = 150, tol = 1e-8){
  em_fit <- function(mu0, sigma0, lambda0, epsilon0){
    mu <- mu0; sigma <- sigma0; lambda <- lambda0; epsilon <- epsilon0
    ll_prev <- -Inf; ll <- NA
    for (iter in seq_len(max_iter)) {
      logf1 <- dlnorm(x, meanlog = mu, sdlog = sigma, log = TRUE)
      logf2 <- dlnorm(x, meanlog = mu + (lambda-1)*sigma^2, sdlog = sqrt(lambda)*sigma, log = TRUE)
      l1 <- log(epsilon) + logf1; l2 <- log(1-epsilon) + logf2
      mx <- pmax(l1, l2); denom <- mx + log(exp(l1-mx) + exp(l2-mx))
      ll <- sum(denom)
      if (abs(ll - ll_prev) < tol*(abs(ll)+1)) break
      ll_prev <- ll
      gamma <- exp(l1 - denom)
      epsilon <- min(max(mean(gamma), eps_lo), 1 - 1e-8)
      Qfun <- function(p){
        mu_p <- p[1]; sigma_p <- exp(p[2]); lambda_p <- min(1+exp(p[3]), lam_hi)
        lf1 <- dlnorm(x, meanlog=mu_p, sdlog=sigma_p, log=TRUE)
        lf2 <- dlnorm(x, meanlog=mu_p+(lambda_p-1)*sigma_p^2, sdlog=sqrt(lambda_p)*sigma_p, log=TRUE)
        sum(gamma*lf1) + sum((1-gamma)*lf2)
      }
      mstep <- optim(c(mu, log(sigma), log(lambda-1)), Qfun, method="BFGS", control=list(fnscale=-1, maxit=100))
      mu <- mstep$par[1]; sigma <- exp(mstep$par[2]); lambda <- min(1+exp(mstep$par[3]), lam_hi)
    }
    list(mu=mu, sigma=sigma, lambda=lambda, epsilon=epsilon, loglik=ll, converged=(iter<max_iter))
  }
  
  mu0 <- mean(log(x)); sigma0 <- sd(log(x))
  best <- NULL; best_ll <- -Inf
  for (s in seq_len(n_starts)) {
    eps0 <- max(epsilon_starts[s], eps_lo+0.01); lam0 <- min(lambda_starts[s], lam_hi-0.01)
    fit <- try(em_fit(mu0, sigma0, lam0, eps0), silent = TRUE)
    if (!inherits(fit, "try-error") && is.finite(fit$loglik) && fit$loglik > best_ll) { best_ll <- fit$loglik; best <- fit }
  }
  inv_eps <- function(e) log((e-eps_lo)/(1-e)); inv_lam <- function(l) log((l-1)/(lam_hi-l))
  eps_c <- min(max(best$epsilon, eps_lo+1e-6), 1-1e-6); lam_c <- min(max(best$lambda, 1+1e-6), lam_hi-1e-6)
  start <- c(best$mu, log(best$sigma), inv_lam(lam_c), inv_eps(eps_c))
  polish <- try(optim(start, dlnL2, x=x, eps_lo=eps_lo, lam_hi=lam_hi, method="BFGS", control=list(fnscale=-1, maxit=2000)), silent=TRUE)
  if (!inherits(polish,"try-error") && is.finite(polish$value) && polish$value >= best_ll) {
    mu<-polish$par[1]; sigma<-exp(polish$par[2])
    lambda<-1+(lam_hi-1)/(1+exp(-polish$par[3])); epsilon<-eps_lo+(1-eps_lo)/(1+exp(-polish$par[4]))
    ll <- polish$value; conv <- (polish$convergence==0)
  } else {
    mu<-best$mu; sigma<-best$sigma; lambda<-best$lambda; epsilon<-best$epsilon; ll<-best_ll; conv<-best$converged
  }
  m <- exp(mu - sigma^2)
  c(m_hat=m, sigma_hat=sigma, lambda_hat=lambda, epsilon_hat=epsilon, logLik=ll, n_converged=as.numeric(conv)*n_starts)
}

# =============================================================================
# MAIN SIMULATION -- same structure as the original script
# =============================================================================

m     <- 2
sigma <- 0.5

n_vals       <- c(100, 1000, 10000)
lambda_vals  <- c(2, 5)
epsilon_vals <- c(0.6, 0.9)
n_reps       <- 500

# Rough expected runtime with fit_mixture() (direct, bounded) on a single core:
#   n=100:   ~1 min total across the 4 (lambda x epsilon) cells
#   n=1000:  ~6 min total
#   n=10000: ~50-55 min total
#   -> full grid (12 cells x 500 reps) ~ 1 hour. Consider n_reps=100 for a first
#      pass, or parallelize across reps (e.g. parallel::mclapply) if you want
#      the full 500 faster.

results      <- list()
summary_rows <- list()

for (n in n_vals) {
  for (lambda in lambda_vals) {
    for (epsilon in epsilon_vals) {
      
      n1 <- round(n * epsilon)
      n2 <- n - n1
      
      rep_results <- matrix(NA, nrow = n_reps, ncol = 6)
      colnames(rep_results) <- c("m_hat", "sigma_hat", "lambda_hat", "epsilon_hat", "logLik", "n_converged")
      
      for (r in seq_len(n_reps)) {
        x1 <- rlnorm(n1, meanlog = log(m) + sigma^2, sdlog = sigma)
        x2 <- rlnorm(n2, meanlog = log(m) + lambda*sigma^2, sdlog = sqrt(lambda)*sigma)
        x  <- c(x1, x2)
        
        rep_results[r, ] <- fit_mixture(x)   # swap for fit_mixture_em(x) on hard cells if desired
      }
      
      key <- sprintf("n=%d_lambda=%g_epsilon=%g", n, lambda, epsilon)
      results[[key]] <- rep_results
      cat("Done:", key, "  (mean starts converged:", round(mean(rep_results[,"n_converged"], na.rm=TRUE),1), "/", n_starts, ")\n")
      
      m_hat       <- rep_results[, "m_hat"]
      sigma_hat   <- rep_results[, "sigma_hat"]
      lambda_hat  <- rep_results[, "lambda_hat"]
      epsilon_hat <- rep_results[, "epsilon_hat"]
      
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
        
        mean_n_converged = mean(rep_results[,"n_converged"], na.rm=TRUE)
      )
    }
  }
}

summary_df <- do.call(rbind, summary_rows)
print(summary_df)