#install.packages(c("circular", "ggplot2", "gridExtra"))
# Create a directory for saving all outputs
output_dir <- "circular_volatility_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

###########################################
library(circular)

# Set seed for reproducibility
set.seed(123)

# Number of time periods and assets
T <- 100
N <- 3

# Simulate angular returns from von Mises distribution for three assets
mu <- c(pi/4, pi/2, 3*pi/4)  # mean directions
kappa <- c(2, 3, 4)          # concentration parameters
# True values of parameters used in the simulation
true_mu <- c(pi/4, pi/2, 3*pi/4)  # True mean directions
true_kappa <- c(2, 3, 4)          # True concentration parameters


# Generate angular returns for each asset
theta <- matrix(0, nrow = T, ncol = N)
for (i in 1:N) {
  theta[, i] <- rvonmises(T, mu[i], kappa[i])
}

# Convert angles to circular format for plotting
theta_circ <- as.circular(theta)

###########################################

# Function to estimate MLE for von Mises parameters
estimate_vm_params <- function(theta) {
  n <- length(theta)
  C <- mean(cos(theta))
  S <- mean(sin(theta))
  Rbar <- sqrt(C^2 + S^2)
  mu_hat <- atan2(S, C)  # Estimated mean direction
  kappa_hat <- ifelse(Rbar < 0.53, 2 * Rbar + Rbar^3 + 5 * Rbar^5 / 6, 
                      ifelse(Rbar < 0.85, -0.4 + 1.39 * Rbar + 0.43 / (1 - Rbar), 
                             1 / (3 * Rbar - 4 * Rbar^2 + Rbar^3)))  # Estimation of kappa
  
  return(list(mu = mu_hat, kappa = kappa_hat))
}

# Estimate parameters for each asset
vm_params <- apply(theta, 2, estimate_vm_params)
vm_params

# Print the estimated parameters
for (i in 1:N) {
  cat(sprintf("Asset %d - Mean direction (mu): %.4f, Concentration (kappa): %.4f\n", i, vm_params[[i]]$mu, vm_params[[i]]$kappa))
}


#######################################

# Define CVM model with angular cross-correlation
cvm_model <- function(theta, sigma2_prev, alpha0, alpha1, beta1) {
  N <- ncol(theta)
  sigma2 <- rep(0, N)
  
  for (i in 1:N) {
    angular_sum <- 0
    for (j in 1:N) {
      if (i != j) {
        angular_sum <- angular_sum + cos(theta[, i] - theta[, j])
      }
    }
    sigma2[i] <- alpha0 + alpha1 * angular_sum + beta1 * sigma2_prev[i]
  }
  
  return(sigma2)
}

# Initialize previous volatility values (can be zero or some baseline)
sigma2_prev <- rep(0, N)

# Define parameters for the CVM
alpha0 <- 0.01
alpha1 <- 0.05
beta1 <- 0.9

# Compute volatility for each asset using the CVM model
sigma2_est <- cvm_model(theta, sigma2_prev, alpha0, alpha1, beta1)
sigma2_est


##############################

# Rayleigh test for uniformity (Goodness-of-fit test for circular data)
library(circular)

for (i in 1:N) {
  test_result <- rayleigh.test(theta[, i])
  print(test_result)
  cat(sprintf("Goodness-of-fit Rayleigh test for Asset %d: p-value = %.4f\n", i, test_result$p.value))
}

#############################################

library(ggplot2)
library(gridExtra)

# Convert angular returns to degrees for better visualization
theta_deg <- theta * 180 / pi
theta_df <- data.frame(time = 1:T, Asset1 = theta_deg[, 1], Asset2 = theta_deg[, 2], Asset3 = theta_deg[, 3])

# Plot angular returns for each asset
p1 <- ggplot(theta_df, aes(x = time)) +
  geom_line(aes(y = Asset1), color = "blue") + labs(title = "Angular Returns of Asset 1", y = "Degrees")
p2 <- ggplot(theta_df, aes(x = time)) +
  geom_line(aes(y = Asset2), color = "green") + labs(title = "Angular Returns of Asset 2", y = "Degrees")
p3 <- ggplot(theta_df, aes(x = time)) +
  geom_line(aes(y = Asset3), color = "red") + labs(title = "Angular Returns of Asset 3", y = "Degrees")

# Save and print the plots
png(filename = file.path(output_dir, "angular_returns.png"), width = 1000, height = 800)
grid.arrange(p1, p2, p3, ncol = 1)
dev.off()

####################################################

vol_df <- data.frame(Asset = 1:N, Volatility = sigma2_est)

# Plot the estimated volatilities
vol_plot <- ggplot(vol_df, aes(x = factor(Asset), y = Volatility)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Estimated Volatility for Each Asset", x = "Asset", y = "Volatility")

# Save and print the volatility plot
png(filename = file.path(output_dir, "volatility_estimates.png"))
print(vol_plot)
dev.off()

####################################################

# Save estimated parameters to a CSV
param_df <- data.frame(Asset = 1:N,
                       MeanDirection = sapply(vm_params, function(x) x$mu),
                       Concentration = sapply(vm_params, function(x) x$kappa))
write.csv(param_df, file = file.path(output_dir, "estimated_parameters.csv"), row.names = FALSE)

# Save estimated volatilities to a CSV
write.csv(vol_df, file = file.path(output_dir, "volatility_estimates.csv"), row.names = FALSE)

##################################################


#######################################################


# Error between estimated and true parameters
mu_error <- sapply(vm_params, function(x) x$mu) - true_mu
kappa_error <- sapply(vm_params, function(x) x$kappa) - true_kappa

# Create a dataframe for parameter error
error_df <- data.frame(Asset = 1:N,
                       Mu_Error = mu_error,
                       Kappa_Error = kappa_error)

# Plot the estimation errors for Mu and Kappa
library(ggplot2)

mu_error_plot <- ggplot(error_df, aes(x = factor(Asset), y = Mu_Error)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(title = "Error in Estimated Mean Direction (Mu)", x = "Asset", y = "Mu Error")

kappa_error_plot <- ggplot(error_df, aes(x = factor(Asset), y = Kappa_Error)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Error in Estimated Concentration (Kappa)", x = "Asset", y = "Kappa Error")

# Save and display the plots
png(filename = file.path(output_dir, "parameter_estimation_error.png"), width = 1000, height = 600)
grid.arrange(mu_error_plot, kappa_error_plot, ncol = 2)
dev.off()


##################

# Simulate true volatilities for comparison (e.g., from GARCH(1,1) model)
true_sigma2 <- matrix(0.02, nrow = T, ncol = N)
for (t in 2:T) {
  true_sigma2[t, ] <- alpha0 + alpha1 * apply(cos(theta[t, ] - theta[t - 1, ]), 2, sum) + beta1 * true_sigma2[t - 1, ]
}

# Predicted vs Actual Volatility Plot
predicted_vs_actual_df <- data.frame(Time = 1:T,
                                     Predicted_Asset1 = sigma2_est[1],
                                     Predicted_Asset2 = sigma2_est[2],
                                     Predicted_Asset3 = sigma2_est[3],
                                     Actual_Asset1 = true_sigma2[, 1],
                                     Actual_Asset2 = true_sigma2[, 2],
                                     Actual_Asset3 = true_sigma2[, 3])

# Plot predicted vs actual volatility for each asset
p1_vol <- ggplot(predicted_vs_actual_df, aes(x = Time)) +
  geom_line(aes(y = Actual_Asset1), color = "red", linetype = "dashed") +
  geom_line(aes(y = Predicted_Asset1), color = "blue") +
  labs(title = "Predicted vs Actual Volatility (Asset 1)", y = "Volatility")

p2_vol <- ggplot(predicted_vs_actual_df, aes(x = Time)) +
  geom_line(aes(y = Actual_Asset2), color = "red", linetype = "dashed") +
  geom_line(aes(y = Predicted_Asset2), color = "blue") +
  labs(title = "Predicted vs Actual Volatility (Asset 2)", y = "Volatility")

p3_vol <- ggplot(predicted_vs_actual_df, aes(x = Time)) +
  geom_line(aes(y = Actual_Asset3), color = "red", linetype = "dashed") +
  geom_line(aes(y = Predicted_Asset3), color = "blue") +
  labs(title = "Predicted vs Actual Volatility (Asset 3)", y = "Volatility")

# Save and display the predicted vs actual volatility plots
png(filename = file.path(output_dir, "predicted_vs_actual_volatility.png"), width = 1000, height = 800)
grid.arrange(p1_vol, p2_vol, p3_vol, ncol = 1)
dev.off()


#######################
# Circular histograms (rose plots) for angular returns
library(circular)

# Create circular histograms for each asset
png(filename = file.path(output_dir, "circular_histograms.png"), width = 1000, height = 600)
par(mfrow = c(1, 3))  # Set up a 1x3 layout for the rose plots

rose.diag(as.circular(theta[, 1]), bins = 18, col = "blue", main = "Circular Histogram (Asset 1)")
rose.diag(as.circular(theta[, 2]), bins = 18, col = "green", main = "Circular Histogram (Asset 2)")
rose.diag(as.circular(theta[, 3]), bins = 18, col = "red", main = "Circular Histogram (Asset 3)")

dev.off()


# Goodness-of-fit results table
gof_results <- data.frame(Asset = 1:N,
                          P_value = sapply(1:N, function(i) rayleigh.test(theta[, i])$p.value))

# Print and save the goodness-of-fit results
write.csv(gof_results, file = file.path(output_dir, "goodness_of_fit_results.csv"), row.names = FALSE)
gof_results

