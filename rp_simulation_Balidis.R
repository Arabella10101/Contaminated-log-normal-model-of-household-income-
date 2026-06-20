#simulated log normal variations



#single model with mu
set.seed(42)

n_samples <-  100000
mu <- 1.6     
sigma <- 0.5  

# Simulate log-normal data
x <- rlnorm(n = n_samples, meanlog = mu, sdlog = sigma)


#likelihood function of log normal with parameters restricted
dlnL <- function(par, x){
  mu <- exp(par[1]) #restrict parameters
  sigma <- exp(par[2])
  return(sum(dlnorm(x, meanlog = mu, sdlog = sigma, log = T)))
}

#initialize values
start <- log(c(mu,sigma))

#using optim to get mles
est <- optim(par = start, fn = dlnL, x=x, control = list(fnscale=-1))

#extract estimated parameters
par <- est$par
mu_est <- exp(par[1])
print(mu_est)
sigma_est <- exp(par[2])
print(sigma_est)
loglike <- est$value
print(loglike)


#contaminated model reparamterized
m <- 2
sigma <- 0.5
lambda <- 5
epsilon <- 0.6

#How many come from each distribution
n1 <- 60000
n2 <- 40000

x1 <- rlnorm(n = n1, meanlog = log(m) + sigma, sdlog = sigma)
x2 <- rlnorm(n = n2, meanlog = log(m) + lambda * sigma, sdlog = lambda * sigma)

#Concat
x <- c(x1, x2)

dlnL2 <- function(par, x){
  m     <- exp(par[1]) 
  sigma <- exp(par[2])
  lambda <- exp(par[3]) + 1 #bound lambda > 1
  # Logit transformation to keep epsilon in (0, 1)
  epsilon <- 1 / (1 + exp(-par[4])) 
  
  return(sum(log(epsilon*dlnorm(x, meanlog=log(m)+sigma, sdlog=sigma) + (1-epsilon)*dlnorm(x, meanlog=log(m)+lambda*sigma, sdlog=lambda*sigma))))
}

# initialization 
# Use logit for epsilon: log(eps / (1-eps))
start <- c(log(m), log(sigma), log(lambda-1), log(epsilon / (1 - epsilon)))

est <- optim(par = start, fn = dlnL2, x = x, control = list(fnscale = -1))

final_params <- exp(est$par)
final_params[3] <- exp(est$par[3])+1
final_params[4] <- 1 / (1 + exp(-est$par[4])) # Apply inverse logit for epsilon
print(final_params)
print(est$convergence) #it converges

#so it doesnt work better when i just use mu 

#Define true parameters on the log-scale
true_pars <- c(log(m), log(sigma), log(lambda-1), log(epsilon / (1 - epsilon)))

# 2. Compare the Log-Likelihoods
cat("Log-Likelihood with TRUE parameters:", dlnL2(true_pars, x), "\n")
cat("Log-Likelihood with OPTIMIZED parameters:", est$value, "\n")

#via troubleshooting why not giving close results to initial parameters 
#opt > true then optim found parameters that better map the sample then org









#iteration to check estimates consistent across diff samples
n_iterations <- 50
results <- matrix(NA, nrow = n_iterations, ncol = 4)
colnames(results) <- c("m", "sigma", "lambda", "epsilon")

for(i in 1:n_iterations) {
  # Generate a new random sample for each iteration
  # Decide how many come from each distribution
  n1 <- 60000
  n2 <- 40000
  
  # Generate directly
  x1 <- rlnorm(n = n1, meanlog = log(m) + sigma, sdlog = sigma)
  x2 <- rlnorm(n = n2, meanlog = log(m) + lambda * sigma, sdlog = lambda * sigma)
  
  # Combine them into one single mixture vector of length 100,000
  x_sim <- c(x1, x2)
  
  # Run optim with original start parameters
  est <- optim(par = start, fn = dlnL2, x = x_sim, control = list(fnscale = -1))
  
  # Back-transform and store results
  params <- exp(est$par)
  params[3] <- exp(est$par[3])+1
  params[4] <- 1 / (1 + exp(-est$par[4]))
  results[i, ] <- params
}

#Analyze the averages
colMeans(results) #parameters optim maps
apply(results, 2, sd) #model is stable


#using previously estimated parameters as new start
results2 <- matrix(NA, nrow = n_iterations, ncol = 4)
colnames(results2) <- c("m", "sigma", "lambda", "epsilon")

#Set the initial starting point
current_start <- c(log(m), log(sigma), log(lambda-1), log(epsilon / (1 - epsilon)))

for(i in 1:n_iterations) {
  # Run optim on the same sample
  est <- optim(par = current_start, fn = dlnL2, x = x, control = list(fnscale = -1))
  
  # Update current_start
  current_start <- est$par
  
  #store results
  params <- exp(est$par)
  params[3] <- exp(est$par[3])+1
  params[4] <- 1 / (1 + exp(-est$par[4]))
  results2[i, ] <- params
}

print(results2)
#consistently gives same result



















