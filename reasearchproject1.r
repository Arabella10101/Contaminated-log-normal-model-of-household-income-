library(readr)
library(dplyr)
library(moments)
library(ggplot2)
library(scales)
library(tidyr)


options(scipen = 999) #turn off sci notation 

households_data <- read_csv("C:/Users/arabe/OneDrive/Documents/Research_Project/Fact_IES2023_Households.csv")

household_income <- households_data %>% 
  select(INCOME, INCOME_PCP)

# Basic summary (Min, Max, Mean, Quartiles)
summary(household_income)

# Standard Deviation and Skewness
sd_income <- sd(household_income$INCOME, na.rm = TRUE)
skew_income <- skewness(household_income$INCOME, na.rm = TRUE)

print(paste("Standard Deviation:", sd_income))
print(paste("Skewness:", skew_income))


# Histogram of Income
ggplot(household_income, aes(x = INCOME)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 50) +
  theme_minimal() +
  labs(title = "Distribution of Household Income", x = "Annual Income", y = "Frequency")

# Log-transformed Density Plot (helps visualize skewed data)
ggplot(household_income, aes(x = INCOME)) +
  geom_density(fill = "orange", alpha = 0.5) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Log-Transformed Income Density", x = "Income (Log Scale)", y = "Density")

# Calculating specific percentiles
percentiles <- quantile(household_income$INCOME, probs = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE)
print(percentiles)

# Boxplot to identify outliers
ggplot(household_income, aes(y = INCOME)) +
  geom_boxplot(fill = "lightgreen") +
  theme_minimal() +
  labs(title = "Boxplot of Income (Identifying Outliers)")















#Simulation stuff

# 1. Define your parameters
target_mode <- 45000  # The peak of your income distribution (e.g., R45,000)
sigma <- 0.8          # The "spread" or inequality (usually 0.5 to 1.2 for income)
n_obs <- 10000        # Number of households to simulate

# 2. Calculate the required mu for that mode
mu <- log(target_mode) + sigma^2

# 3. Generate the simulated data
set.seed(123) # For reproducibility
simulated_income <- rlnorm(n_obs, meanlog = mu, sdlog = sigma)

# 4. Create the Dataframe
sim_df <- data.frame(Income = simulated_income)

# 5. Plot to verify the peak (Mode)
ggplot(sim_df, aes(x = Income)) + 
  geom_density(fill = "steelblue", alpha = 0.6) + 
  geom_vline(xintercept = target_mode, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(labels = label_comma()) + # Remove limits here
  coord_cartesian(xlim = c(0, 300000)) +       # Zoom here instead
  theme_minimal() +
  labs(title = paste("Simulated Log-Normal (Mode =", target_mode, ")"),
       x = "Annual Income (R)",
       y = "Density")






#sim n theoretical

# 1. Setup Parameters (Using the Mode-based parameters from before)
target_mode <- 45000
sigma <- 0.8
mu <- log(target_mode) + sigma^2
n_obs <- 10000

# 2. Generate the Simulated Data
set.seed(123)
sim_data <- data.frame(Income = rlnorm(n_obs, meanlog = mu, sdlog = sigma))

# 1. Transform your data first
sim_data$log_income <- log10(sim_data$Income)

# 2. Calculate the mean and sd of the LOGGED data
mu_log10 <- mean(sim_data$log_income)
sd_log10 <- sd(sim_data$log_income)

# 3. Plot using the logged data
# Use the ORIGINAL simulated income (not the log10 one)
ggplot(sim_data, aes(x = Income)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 100, fill = "steelblue", color = "white", alpha = 0.5) +
  
  # Increase n to 1000 for a smooth, mathematical curve
  stat_function(fun = dlnorm, 
                args = list(meanlog = mu, sdlog = sigma), 
                color = "darkred", size = 1.2, n = 1000) +
  
  coord_cartesian(xlim = c(0, 500000)) +
  scale_x_continuous(labels = label_comma()) +
  theme_minimal() +
  labs(title = "Simulated and theoretical log normal",
       x = "Annual Income (R)",
       y = "Density")



#multiple dist

# 1. Setup Parameters
target_mode <- 45000
sigmas <- c(0.5, 0.8, 1.1) # Low, Medium, and High inequality
n_obs <- 10000

# 2. Run Simulations
sim_list <- lapply(sigmas, function(s) {
  mu <- log(target_mode) + s^2
  data.frame(Income = rlnorm(n_obs, meanlog = mu, sdlog = s), Sigma = as.factor(s))
})

multi_sim_df <- bind_rows(sim_list)

# 3. Plot the Comparison
ggplot(multi_sim_df, aes(x = Income, fill = Sigma, color = Sigma)) +
  # Using geom_density creates the shaded areas based on the simulated data
  geom_density(alpha = 0.3, size = 1) +
  
  # Focus on the relevant income range
  coord_cartesian(xlim = c(0, 600000)) +
  scale_x_continuous(labels = label_comma()) +
  theme_minimal() +
  labs(title = "Simulated Log Normal Income Distributions by Sigma",
       subtitle = paste("All samples generated with a Mode of", label_comma()(target_mode)),
       x = "Annual Income (R)",
       y = "Density")












#graph from simulated params - from rpsim
#idk how useful this is but its made so its here now

# 1. Define your parameters
# True Parameters
m_true <- 2.0
sigma_true <- 0.5
lambda_true <- 5.0
epsilon_true <- 0.6


# 2. Back-transform to original scale
# exp() for m, sigma, lamda; inv-logit for epsilon
m_est     <- exp(est$par[1])
sigma_est <- exp(est$par[2])
lambda_est <- exp(est$par[3])
eps_est   <- 1 / (1 + exp(-est$par[4]))

# Find the 99th percentile of your data to set a smart limit
x_limit <- quantile(x, 0.99)

#plot with the limit
hist(x, breaks = 100, freq = FALSE, col = "grey80", border = "white",
     main = "True vs. Estimated Mixture Density",
     xlab = "x", ylab = "Density", 
     xlim = c(0, x_limit), # <--- This limits the x-axis to the visible range
     ylim = c(0, 0.15))

#add the curves
curve(mixture_pdf(x, m_true, sigma_true, lambda_true, epsilon_true), 
      add = TRUE, col = "blue", lwd = 2, xlim = c(0, x_limit))

curve(mixture_pdf(x, m_est, sigma_est, lambda_est, epsilon_est), 
      add = TRUE, col = "red", lwd = 2, lty = 2, xlim = c(0, x_limit))

#legend
legend("topright", legend = c("True Params", "Estimated Params"), 
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))

