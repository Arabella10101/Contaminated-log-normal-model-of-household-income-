#contaminated model simulation
library(moments) 

set.seed(123)

dlnL2 <- function(par, x){
  m     <- exp(par[1]) 
  sigma <- exp(par[2])
  lambda <- exp(par[3]) + 1 #bound lambda > 1
  # Logit transformation to keep epsilon in (0, 1)
  epsilon <- 1 / (1 + exp(-par[4])) 
  
  return(sum(log(epsilon*dlnorm(x, meanlog=log(m)+sigma^2, sdlog=sigma) + (1-epsilon)*dlnorm(x, meanlog=log(m)+lambda*sigma^2, sdlog=sqrt(lambda)*sigma))))
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

results <- list()

for (n in n_vals) {
  for (lambda in lambda_vals) {
    for (epsilon in epsilon_vals) {
      
      n1 <- round(n * epsilon)
      n2 <- n - n1
      
      rep_results <- matrix(NA, nrow = n_reps, ncol = 5)
      colnames(rep_results) <- c("m_hat", "sigma_hat", "lambda_hat", "epsilon_hat", "logLik")
      
      for (r in seq_len(n_reps)) {
        
        x1 <- rlnorm(n1, meanlog = log(m) + sigma^2,        sdlog = sigma)
        x2 <- rlnorm(n2, meanlog = log(m) + lambda*sigma^2, sdlog = sqrt(lambda)*sigma)
        x  <- c(x1, x2)
        
        best_val <- -Inf
        best_par <- NULL
        
        for (s in seq_len(n_starts)) {
          
          #The shiftfact multiplies each starting parameter by a random factor centered at 1 but varying roughly 0.3
          #so each of the optimization attempts begin at a slightly different point
          shiftfact <- if (s == 1) c(1, 1, 1, 1) else exp(rnorm(4, sd = 0.3))
          
          start <- c(
            log(m) * shiftfact[1],
            log(sigma) * shiftfact[2],
            log(max(lambda - 1, 0.01)) * shiftfact[3],
            log(epsilon / (1 - epsilon)) * shiftfact[4]
          )
          
          est <- try(
            optim(par = start, fn = dlnL2, x = x, control = list(fnscale = -1)),
            silent = TRUE
          )
          
          if (!inherits(est, "try-error") && est$value > best_val) {
            best_val <- est$value
            best_par <- est$par
          }
        }
        
        #extract params and back transform
        final_params <- exp(best_par)
        final_params[3] <- exp(best_par[3]) + 1
        final_params[4] <- 1 / (1 + exp(-best_par[4]))
        
        rep_results[r, ] <- c(final_params, best_val)
      }
      
      key <- sprintf("n=%d_lambda=%g_epsilon=%g", n, lambda, epsilon)
      results[[key]] <- rep_results
      cat("Done:", key, "\n")
    }
  }
}

##summary across the 100 replicates per scenario
summary_table <- do.call(rbind, lapply(names(results), function(key) {
  mat <- results[[key]]
  data.frame(
    scenario     = key,
    m_mean       = mean(mat[, "m_hat"]),
    sigma_mean   = mean(mat[, "sigma_hat"]),
    lambda_mean  = mean(mat[, "lambda_hat"]),
    epsilon_mean = mean(mat[, "epsilon_hat"]),
    m_sd         = sd(mat[, "m_hat"]),
    sigma_sd     = sd(mat[, "sigma_hat"]),
    lambda_sd    = sd(mat[, "lambda_hat"]),
    epsilon_sd   = sd(mat[, "epsilon_hat"])
  )
}))

print(summary_table)

#best parameter set (highest log likelihood) per scenario
best_table <- do.call(rbind, lapply(names(results), function(key) {
  mat <- results[[key]]
  best_row <- mat[which.max(mat[, "logLik"]), ]
  data.frame(scenario = key, t(best_row))
}))

print(best_table)




#graphs

#Kurtosis of cmLN

raw_moment_lnorm <- function(k, meanlog, sdlog) {
  exp(k * meanlog + k^2 * sdlog^2 / 2)
}

cmLN_kurtosis <- function(sigma, lambda, epsilon) {
  meanlog1 <- sigma^2
  sdlog1   <- sigma
  meanlog2 <- lambda * sigma^2
  sdlog2   <- sqrt(lambda) * sigma
  
  M <- sapply(1:4, function(k) {
    epsilon * raw_moment_lnorm(k, meanlog1, sdlog1) +
      (1 - epsilon) * raw_moment_lnorm(k, meanlog2, sdlog2)
  })
  
  mu1 <- M[1]
  var <- M[2] - mu1^2
  mu4 <- M[4] - 4*mu1*M[3] + 6*mu1^2*M[2] - 3*mu1^4
  
  mu4 / var^2
}

LN_kurtosis <- function(sigma) {
  exp(4*sigma^2) + 2*exp(3*sigma^2) + 3*exp(2*sigma^2) - 6 + 3
}

sigma_seq <- seq(0.01, 2, length.out = 200)

epsilon_vals <- c(0.4, 0.5, 0.9)
lambda_vals  <- c(2, 5, 10)

#different epsilon
png("C:/Users/arabe/Documents/Research_Project/Pictures/cmLN_kurtosis_plots_a.png",
    width = 700, height = 600, res = 150)
par(mar = c(4, 4, 2, 1))

vals_a <- sapply(epsilon_vals, function(e) sapply(sigma_seq, function(s) cmLN_kurtosis(s, 5, e)))
base_a <- sapply(sigma_seq, LN_kurtosis)
plot(sigma_seq, base_a, type = "l", lwd = 2, ylim = c(0, 100),
     xlab = expression(sigma), ylab = "Kurtosis", main = "")
for (i in seq_along(epsilon_vals)) lines(sigma_seq, vals_a[, i], lty = i + 1)
legend("topright", legend = c("LN", paste0("cmLN, \u03b5=", epsilon_vals)), lty = 1:4, bty = "n", cex = 0.8)

dev.off()

#different lambda
png("C:/Users/arabe/Documents/Research_Project/Pictures/cmLN_kurtosis_plots_b.png",
    width = 700, height = 600, res = 150)
par(mar = c(4, 4, 2, 1))

vals_b <- sapply(lambda_vals, function(l) sapply(sigma_seq, function(s) cmLN_kurtosis(s, l, 0.6)))
base_b <- sapply(sigma_seq, LN_kurtosis)
plot(sigma_seq, base_b, type = "l", lwd = 2, ylim = c(0, 100),
     xlab = expression(sigma), ylab = "Kurtosis", main = "")
for (i in seq_along(lambda_vals)) lines(sigma_seq, vals_b[, i], lty = i + 1)
legend("topright", legend = c("LN", paste0("cmLN, \u03bb=", lambda_vals)), lty = 1:4, bty = "n", cex = 0.8)

dev.off()









#Mode-parameterized lognormal density

dmLN <- function(x, m, sigma) {
  dlnorm(x, meanlog = log(m) + sigma^2, sdlog = sigma)
}

x_seq <- seq(0.01, 20, length.out = 500)

#different modes, sigma fixed
sigma_fixed <- 0.5
m_vals <- c(1, 3, 6)

png("C:/Users/arabe/Documents/Research_Project/Pictures/mLN_density_modes.png",
    width = 700, height = 600, res = 150)
par(mar = c(4, 4, 2, 1))

plot(x_seq, dmLN(x_seq, m_vals[1], sigma_fixed), type = "l", lwd = 2, lty = 1,
     xlab = "x", ylab = "Density", main = "")
for (i in 2:length(m_vals)) {
  lines(x_seq, dmLN(x_seq, m_vals[i], sigma_fixed), lty = i)
}
legend("topright", legend = paste0("m=", m_vals), lty = 1:length(m_vals), bty = "n", cex = 0.8)

dev.off()

## ---- Panel 2: different sigmas, m fixed ----
m_fixed <- 3
sigma_vals <- c(0.25, 0.5, 1)

png("C:/Users/arabe/Documents/Research_Project/Pictures/mLN_density_sigmas.png",
    width = 700, height = 600, res = 150)
par(mar = c(4, 4, 2, 1))

plot(x_seq, dmLN(x_seq, m_fixed, sigma_vals[1]), type = "l", lwd = 2, lty = 1,
     xlab = "x", ylab = "Density", main = "")
for (i in 2:length(sigma_vals)) {
  lines(x_seq, dmLN(x_seq, m_fixed, sigma_vals[i]), lty = i)
}
legend("topright", legend = paste0("\u03c3=", sigma_vals), lty = 1:length(sigma_vals), bty = "n", cex = 0.8)

dev.off()





