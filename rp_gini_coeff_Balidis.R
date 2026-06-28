## ----------------------------------------------------------
## Gini coefficient for the contaminated mode-parametrized
## log-normal distribution, eqs (2)-(4)
## ----------------------------------------------------------

library(pracma)   # provides erf()

# Parameters (m, sigma, lambda, epsilon)
m       <- 33496.2178306
sigma   <- 0.8235000
lambda_ <- 1.9491272
epsilon <- 0.5043277

sigma2 <- sigma^2

## H(x; m, sigma^2): CDF of the mode-parametrized log-normal
H <- function(x, m, s2) {
  0.5 + 0.5 * erf((log(x) - log(m) - s2) / sqrt(2 * s2))
}

## F(x; eps, m, sigma^2, lambda): contaminated CDF, eq (3)
F_contam <- function(x, eps, m, s2, lam) {
  eps * H(x, m, s2) + (1 - eps) * H(x, m, lam * s2)
}

## E(X; eps, m, sigma^2, lambda): contaminated mean, eq (4)
EX_contam <- function(eps, m, s2, lam) {
  m * (eps * exp(1.5 * s2) + (1 - eps) * exp(1.5 * lam * s2))
}

## Gini coefficient, eq (2)
## Uses a log-scale substitution x = exp(t) for numerical stability,
## since with large m the survival function (1-F)^2 decays very slowly
## on the raw x-scale and integrate() over (0, Inf) becomes unreliable.
gini <- function(eps, m, s2, lam, tol = 1e-12) {
  
  mu <- EX_contam(eps, m, s2, lam) #computes the denom
  
  survival_sq <- function(x) (1 - F_contam(x, eps, m, s2, lam))^2 #defines the integral
  
  # find a finite upper cutoff X_hi where the survival function is negligible
  #Integrating literally to Inf is numerically unreliable here because m 
  #is huge (33,496), so the tail decays extremely slowly on the raw scale. Instead, uniroot solves for 
  #the point X_hi where the survival function 1-F(x) drops to tol (1e-12) — effectively zero. Beyond X_hi, 
  #the contribution to the integral is negligible, so it's safe to treat X_hi as "infinity" for integration 
  #purposes.
  target <- function(x) (1 - F_contam(x, eps, m, s2, lam)) - tol
  X_hi <- uniroot(target, lower = m, upper = m * 1e8)$root
  
  # integrand on the log scale: x = exp(t), dx = x dt
  #Instead of integrating over x directly, substitute x = e^t 
  #(so t = log(x)). Since dx = x·dt, the integrand becomes survival_sq(x) * x. This works better 
  #numerically because the distribution lives on a log scale (it's log-normal-based) — spacing the 
  #integration points evenly in log(x) matches where the function actually varies, rather than wasting 
  #resolution on a linear x-grid that's mostly flat.
  integrand_t <- function(t) {
    x <- exp(t)
    survival_sq(x) * x
  }
  
  #Performs the actual numerical integration in t-space, from t_lo (corresponding to x ≈ 0) 
  #to t_hi (corresponding to x = X_hi, the cutoff found above)
  t_lo <- log(1e-6)     # effectively x -> 0
  t_hi <- log(X_hi)
  
  res <- integrate(integrand_t, lower = t_lo, upper = t_hi,
                   rel.tol = 1e-10, subdivisions = 1000)
  
  #final calc to get Gini coeff
  G <- 1 - res$value / mu
  list(G = G, abs.error = res$abs.error, X_hi = X_hi, mean = mu)
}

result <- gini(epsilon, m, sigma2, lambda_)

cat("Cutoff X_hi used   :", result$X_hi, "\n")
cat("Mean E(X)          :", result$mean, "\n")
cat("Gini coefficient   :", result$G, "\n")
cat("Integration error  :", result$abs.error, "\n")
