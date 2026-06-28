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

dlnL2 <- function(par, x){
  m     <- exp(par[1]) 
  sigma <- exp(par[2])
  lambda <- exp(par[3])+1
  epsilon <- 1 / (1 + exp(-par[4])) 
  
  return(sum(log(epsilon*dlnorm(x, meanlog=log(m)+sigma^2, sdlog=sigma) + (1-epsilon)*dlnorm(x, meanlog=log(m)+lambda*sigma^2, sdlog=sqrt(lambda)*sigma))))
}


# Estimate m: The peak of the data correct?
get_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
start_m <- get_mode(household_income$INCOME)

#correct?
start_sigma <- sd(log(household_income$INCOME[household_income$INCOME < median(household_income$INCOME)]))


#prospective lambda n epsilon to check functions work
start_lambda <- 2  
start_eps <- 0.9     

init_pars <- c(log(start_m), log(start_sigma), log(start_lambda-1), log(start_eps / (1 - start_eps)))

est <- optim(par = init_pars, fn = dlnL2, x = household_income$INCOME, control = list(fnscale = -1))

final_params <- exp(est$par)
final_params[3] <- exp(est$par[3])+1
final_params[4] <- 1 / (1 + exp(-est$par[4])) 
print(final_params)
print(est$convergence)


# Define the number of steps and the "nudge" magnitude
n_steps <- 50
lambda_step <- 0.1   # Increase lambda by 0.1 each time
eps_step   <- -0.01  # Decrease epsilon by 0.01 each time

results <- matrix(NA, nrow = n_steps, ncol = 7)
colnames(results) <- c("m", "sigma", "lambda", "epsilon", "curr lambda", "curr epsilon", "logLik")

# Starting baseline
base_lambda <- 1
base_eps   <- 0.99

for(i in 1:n_steps) {
  # Calculate current starting nudge
  curr_L <- base_lambda + (i * lambda_step)
  curr_E <- base_eps   + (i * eps_step)
  
  # Ensure epsilon stays within bounds (0, 1) to avoid math errors
  curr_E <- max(0.01, min(0.99, curr_E))
  
  init_pars <- c(log(start_m), log(start_sigma), log(curr_L - 1), log(curr_E / (1 - curr_E)))
  
  est <- optim(par = init_pars, fn = dlnL2, x = household_income$INCOME, control = list(fnscale = -1))
  
  params <- exp(est$par)
  params[3] <- exp(est$par[3]) + 1
  params[4] <- 1 / (1 + exp(-est$par[4]))
  
  results[i, ] <- c(params, curr_L, curr_E, est$value)
}

# Pick the row with the highest log-likelihood
best_idx <- which.max(results[, "logLik"])
best_result <- results[best_idx, ]

best_result


#using previously estimated parameters as new start to check convergence
n_iterations <- 50
results2 <- matrix(NA, nrow = n_iterations, ncol = 4)
colnames(results2) <- c("m", "sigma", "lambda", "epsilon")

#based on previous result 
current_start <- est$par 

for(i in 1:n_iterations) {
  est <- optim(par = current_start, fn = dlnL2, x = household_income$INCOME, control = list(fnscale = -1))
  
  current_start <- est$par
  
  params <- exp(est$par)
  params[3] <- exp(est$par[3])+1
  params[4] <- 1 / (1 + exp(-est$par[4]))
  results2[i, ] <- params
}

print(results2)











