library(ggplot2)
library(dplyr)

out_dir <- "C:\\Users\\arabe\\Documents\\Research_Project\\Pictures"

# --- Plot 1: Varying Mode ---
x <- seq(0.001, 20, length.out = 2000)
sigma2 <- 0.25
modes <- c(1, 3, 6)
lty_vec <- c("solid", "dashed", "dotted")
lwd_vec <- c(2.5, 1.5, 1.5)

dens <- lapply(modes, function(m) {
  mu <- log(m) + sigma2
  dlnorm(x, meanlog = mu, sdlog = sqrt(sigma2))
})

png(file.path(out_dir, "lognormal_varying_mode.png"),
    width = 2000, height = 1500, res = 300)

plot(x, dens[[1]], type = "l", lwd = lwd_vec[1], lty = lty_vec[1],
     xlim = c(0, 20), ylim = c(0, 0.75),
     xlab = "x", ylab = "f(x)",
     main = expression("Log-Normal: Varying Mode " * (sigma^2 == 0.25)),
     bty = "l")
for (i in 2:3) {
  lines(x, dens[[i]], lwd = lwd_vec[i], lty = lty_vec[i])
}
legend("topright", legend = paste0("m=", modes),
       lty = lty_vec, lwd = lwd_vec, bty = "n")

dev.off()

# --- Plot 2: Varying sigma^2 ---
x <- seq(0.001, 20, length.out = 2000)
m <- 2
sigma2s <- c(0.25, 0.5, 1, 2)
lty_vec <- c("solid", "dashed", "dotted", "dotdash")
lwd_vec <- c(2.2, 1.5, 1.5, 1.5)

dens <- lapply(sigma2s, function(s2) {
  mu <- log(m) + s2
  dlnorm(x, meanlog = mu, sdlog = sqrt(s2))
})

png(file.path(out_dir, "lognormal_varying_sigma2.png"),
    width = 2000, height = 1500, res = 300)

plot(x, dens[[1]], type = "l", lwd = lwd_vec[1], lty = lty_vec[1],
     xlim = c(0, 20), ylim = c(0, 0.75),
     xlab = "x", ylab = "f(x)",
     main = expression("Log-Normal: Varying " * sigma^2 * " (m = 2)"),
     bty = "l")
for (i in 2:4) {
  lines(x, dens[[i]], lwd = lwd_vec[i], lty = lty_vec[i])
}
legend("topright", legend = paste0("σ² = ", sigma2s),
       lty = lty_vec, lwd = lwd_vec, bty = "n")

dev.off()











# Log-Normal Mixture: Varying Lambda
# Mode (m) and sigma^2 held fixed; lambda controls how much more
# spread out the second ("contaminating") component is.
#
# Mixture weights come from your n1/n2 sample sizes (60000/40000 -> 0.6/0.4)

m     <- 3          # fixed mode
sigma <- sqrt(0.25)  # fixed sigma^2 = 0.25
p1    <- 0.6         # weight on component 1 (n1 / (n1+n2))
p2    <- 0.4         # weight on component 2 (n2 / (n1+n2))

# Analytic mixture density (mode-anchored lognormal mixture)
f <- function(x, lambda) {
  mu1 <- log(m) + sigma^2
  sd1 <- sigma
  mu2 <- log(m) + lambda * sigma^2
  sd2 <- sqrt(lambda) * sigma
  p1 * dlnorm(x, meanlog = mu1, sdlog = sd1) +
    p2 * dlnorm(x, meanlog = mu2, sdlog = sd2)
}

x <- seq(0.001, 20, length.out = 1000)

lambdas <- c(2, 5, 10)
ltys    <- c(1, 2, 3)     # solid, dashed, dotted
lwds    <- c(2, 1, 1)     # thicker line for lambda=1, matching original style

y_list <- lapply(lambdas, function(l) f(x, l))
ymax   <- max(unlist(y_list))

png(file.path(out_dir, "cmLN_varying_lambda.png"),
    width = 2000, height = 1500, res = 300)

plot(x, y_list[[1]], type = "l", lwd = lwds[1], lty = ltys[1],
     ylim = c(0, ymax * 1.1),
     xlab = "x", ylab = "f(x)",
     main = expression(paste("Contaminated Log Normal: Varying ", lambda,
                             "  (m = 3, ", sigma^2, " = 0.25)")))
lines(x, y_list[[2]], lwd = lwds[2], lty = ltys[2])
lines(x, y_list[[3]], lwd = lwds[3], lty = ltys[3])

legend("topright",
       legend = c(expression(lambda == 1), expression(lambda == 3), expression(lambda == 6)),
       lty = ltys, lwd = lwds, bty = "n")
