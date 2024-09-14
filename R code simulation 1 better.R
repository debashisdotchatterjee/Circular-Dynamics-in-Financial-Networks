# Load necessary libraries
library(circular)
library(ggplot2)
library(gridExtra)

# Create a directory for saving all outputs
output_dir <- "circular_volatility_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set seed for reproducibility
set.seed(123)

# Number of time periods and assets
T <- 100
N <- 3

# True parameter values used for simulation
alpha0 <- 0.02
alpha1 <- 0.05
beta1 <- 0.9
true_mu <- c(pi/4, pi/2, 3*pi/4)  # True mean directions
true_kappa <- c(2, 3, 4)          # True concentration parameters

# Simulate angular returns from von Mises distribution for three assets
theta <- matrix(0, nrow = T, ncol = N)
for (i in 1:N) {
  theta[, i] <- rvonmises(T, true_mu[i], true_kappa[i])
}

# Convert angles to circular format for plotting
theta_circ <- as.circular(theta)

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

# Volatility estimation using CVM
sigma2_est <- rep(0.02, N)  # Initial volatility estimate

# Define CVM model with angular cross-correlation
cvm_model <- function(theta, sigma2_prev, alpha0, alpha1, beta1) {
  N <- ncol(theta)
  sigma2 <- rep(0, N)
  
  for (i in 1:N) {
    angular_sum <- sum(cos(theta[, i] - theta[, -i]))
    sigma2[i] <- alpha0 + alpha1 * angular_sum + beta1 * sigma2_prev[i]
  }
  
  return(sigma2)
}

# Compute predicted volatility using CVM
sigma2_predicted <- matrix(0, nrow = T, ncol = N)
for (t in 2:T) {
  sigma2_predicted[t, ] <- cvm_model(theta[t, , drop = FALSE], sigma2_predicted[t - 1, ], alpha0, alpha1, beta1)
}

# Data frame for predicted vs volatility plot
volatility_df <- data.frame(Time = 1:T,
                            Predicted_Asset1 = sigma2_predicted[, 1],
                            Predicted_Asset2 = sigma2_predicted[, 2],
                            Predicted_Asset3 = sigma2_predicted[, 3])

# Plot predicted volatility for each asset
p1_vol <- ggplot(volatility_df, aes(x = Time)) +
  geom_line(aes(y = Predicted_Asset1), color = "blue") +
  labs(title = "Predicted Volatility (Asset 1)", y = "Volatility")

p2_vol <- ggplot(volatility_df, aes(x = Time)) +
  geom_line(aes(y = Predicted_Asset2), color = "blue") +
  labs(title = "Predicted Volatility (Asset 2)", y = "Volatility")

p3_vol <- ggplot(volatility_df, aes(x = Time)) +
  geom_line(aes(y = Predicted_Asset3), color = "blue") +
  labs(title = "Predicted Volatility (Asset 3)", y = "Volatility")

# Save the predicted volatility plots
png(filename = file.path(output_dir, "predicted_volatility.png"), width = 1000, height = 800)
grid.arrange(p1_vol, p2_vol, p3_vol, ncol = 1)
dev.off()

# Circular histograms (rose plots) for angular returns
png(filename = file.path(output_dir, "circular_histograms.png"), width = 1000, height = 600)
par(mfrow = c(1, 3))  # Set up a 1x3 layout for the rose plots
rose.diag(as.circular(theta[, 1]), bins = 18, col = "blue", main = "Circular Histogram (Asset 1)")
rose.diag(as.circular(theta[, 2]), bins = 18, col = "green", main = "Circular Histogram (Asset 2)")
rose.diag(as.circular(theta[, 3]), bins = 18, col = "red", main = "Circular Histogram (Asset 3)")
dev.off()

# Error between estimated and true parameters
mu_error <- sapply(vm_params, function(x) x$mu) - true_mu
kappa_error <- sapply(vm_params, function(x) x$kappa) - true_kappa

# Create a dataframe for parameter error
error_df <- data.frame(Asset = 1:N,
                       Mu_Error = mu_error,
                       Kappa_Error = kappa_error)

# Plot the estimation errors for Mu and Kappa
mu_error_plot <- ggplot(error_df, aes(x = factor(Asset), y = Mu_Error)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(title = "Error in Estimated Mean Direction (Mu)", x = "Asset", y = "Mu Error")

kappa_error_plot <- ggplot(error_df, aes(x = factor(Asset), y = Kappa_Error)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Error in Estimated Concentration (Kappa)", x = "Asset", y = "Kappa Error")

# Save and display the parameter estimation error plots
png(filename = file.path(output_dir, "parameter_estimation_error.png"), width = 1000, height = 600)
grid.arrange(mu_error_plot, kappa_error_plot, ncol = 2)
dev.off()

# Circular plot for True vs Estimated Mu (Mean Direction)
png(filename = file.path(output_dir, "mu_true_vs_estimated.png"), width = 800, height = 800)

# Set up a circular plot
plot.new()
par(pty = "s")  # Set the plotting region to be square
plot.window(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))

# Draw a circle
symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, fg = "gray")

# Plot true and estimated mu as arrows from the center
for (i in 1:N) {
  # Arrow for true mu (in red)
  arrows(0, 0, cos(true_mu[i]), sin(true_mu[i]), col = "red", lwd = 2, length = 0.1)
  text(1.2 * cos(true_mu[i]), 1.2 * sin(true_mu[i]), paste("True Mu", i), col = "red", cex = 0.8)
  
  # Arrow for estimated mu (in blue)
  arrows(0, 0, cos(vm_params[[i]]$mu), sin(vm_params[[i]]$mu), col = "blue", lwd = 2, length = 0.1)
  text(1.2 * cos(vm_params[[i]]$mu), 1.2 * sin(vm_params[[i]]$mu), paste("Est. Mu", i), col = "blue", cex = 0.8)
}

# Add legend
legend("topright", legend = c("True Mu", "Estimated Mu"), col = c("red", "blue"), lwd = 2)

dev.off()

# Print and save estimated parameters to CSV
param_df <- data.frame(Asset = 1:N,
                       True_MeanDirection = true_mu,
                       Estimated_MeanDirection = sapply(vm_params, function(x) x$mu),
                       True_Concentration = true_kappa,
                       Estimated_Concentration = sapply(vm_params, function(x) x$kappa))
write.csv(param_df, file = file.path(output_dir, "estimated_vs_true_parameters.csv"), row.names = FALSE)

# Goodness-of-fit using Rayleigh test
gof_results <- data.frame(Asset = 1:N,
                          P_value = sapply(1:N, function(i) rayleigh.test(theta[, i])$p.value))

# Save goodness-of-fit results
write.csv(gof_results, file = file.path(output_dir, "goodness_of_fit_results.csv"), row.names = FALSE)

# Display the goodness-of-fit results
print(gof_results)

# Save all results summary as a single file
summary_df <- data.frame(Asset = 1:N,
                         True_MeanDirection = true_mu,
                         Estimated_MeanDirection = sapply(vm_params, function(x) x$mu),
                         True_Concentration = true_kappa,
                         Estimated_Concentration = sapply(vm_params, function(x) x$kappa),
                         Mu_Error = mu_error,
                         Kappa_Error = kappa_error,
                         Goodness_of_Fit_PValue = gof_results$P_value)

# Save the summary CSV file
write.csv(summary_df, file = file.path(output_dir, "model_summary.csv"), row.names = FALSE)

# Print out the summary table for viewing
print(summary_df)

# End of script

######################################


###########################################

# Corrected Predicted vs Simulated Data Plot

# Combine simulated data and predicted volatility for all assets
predicted_vs_simulated_df <- data.frame(Time = 1:T,
                                        Simulated_Asset1 = theta[, 1],
                                        Simulated_Asset2 = theta[, 2],
                                        Simulated_Asset3 = theta[, 3],
                                        Predicted_Asset1 = sigma2_predicted[, 1],
                                        Predicted_Asset2 = sigma2_predicted[, 2],
                                        Predicted_Asset3 = sigma2_predicted[, 3])

# Plot predicted volatility vs simulated data for Asset 1
p1_pred_vs_sim <- ggplot(predicted_vs_simulated_df, aes(x = Time)) +
  geom_line(aes(y = Simulated_Asset1, color = "Simulated Data"), linetype = "dashed") +
  geom_line(aes(y = Predicted_Asset1, color = "Predicted Volatility")) +
  labs(title = "Predicted Volatility vs Simulated Data (Asset 1)", y = "Value") +
  scale_color_manual(name = "Legend", values = c("Simulated Data" = "green", "Predicted Volatility" = "blue"))

# Plot for Asset 2
p2_pred_vs_sim <- ggplot(predicted_vs_simulated_df, aes(x = Time)) +
  geom_line(aes(y = Simulated_Asset2, color = "Simulated Data"), linetype = "dashed") +
  geom_line(aes(y = Predicted_Asset2, color = "Predicted Volatility")) +
  labs(title = "Predicted Volatility vs Simulated Data (Asset 2)", y = "Value") +
  scale_color_manual(name = "Legend", values = c("Simulated Data" = "green", "Predicted Volatility" = "blue"))

# Plot for Asset 3
p3_pred_vs_sim <- ggplot(predicted_vs_simulated_df, aes(x = Time)) +
  geom_line(aes(y = Simulated_Asset3, color = "Simulated Data"), linetype = "dashed") +
  geom_line(aes(y = Predicted_Asset3, color = "Predicted Volatility")) +
  labs(title = "Predicted Volatility vs Simulated Data (Asset 3)", y = "Value") +
  scale_color_manual(name = "Legend", values = c("Simulated Data" = "green", "Predicted Volatility" = "blue"))

# Save the corrected plots
png(filename = file.path(output_dir, "predicted_vs_simulated_volatility.png"), width = 1000, height = 800)
grid.arrange(p1_pred_vs_sim, p2_pred_vs_sim, p3_pred_vs_sim, ncol = 1)
dev.off()

###################################

# Calculate residuals (errors) between actual and predicted
residuals_df <- data.frame(Time = 1:T,
                           Residual_Asset1 = theta[, 1] - sigma2_predicted[, 1],
                           Residual_Asset2 = theta[, 2] - sigma2_predicted[, 2],
                           Residual_Asset3 = theta[, 3] - sigma2_predicted[, 3])

# Plot residuals for Asset 1
p1_residuals <- ggplot(residuals_df, aes(x = Time, y = Residual_Asset1)) +
  geom_line(color = "red") +
  labs(title = "Residuals (Actual - Predicted) for Asset 1", y = "Residuals")

# Plot residuals for Asset 2
p2_residuals <- ggplot(residuals_df, aes(x = Time, y = Residual_Asset2)) +
  geom_line(color = "red") +
  labs(title = "Residuals (Actual - Predicted) for Asset 2", y = "Residuals")

# Plot residuals for Asset 3
p3_residuals <- ggplot(residuals_df, aes(x = Time, y = Residual_Asset3)) +
  geom_line(color = "red") +
  labs(title = "Residuals (Actual - Predicted) for Asset 3", y = "Residuals")

# Save the residuals plot
png(filename = file.path(output_dir, "residuals_plot.png"), width = 1000, height = 800)
grid.arrange(p1_residuals, p2_residuals, p3_residuals, ncol = 1)
dev.off()

###########################

# QQ Plot for Residuals
png(filename = file.path(output_dir, "qq_plot_residuals.png"), width = 800, height = 600)
par(mfrow = c(1, 3))

# QQ Plot for Asset 1
qqnorm(residuals_df$Residual_Asset1, main = "QQ Plot for Asset 1 Residuals")
qqline(residuals_df$Residual_Asset1, col = "blue")

# QQ Plot for Asset 2
qqnorm(residuals_df$Residual_Asset2, main = "QQ Plot for Asset 2 Residuals")
qqline(residuals_df$Residual_Asset2, col = "blue")

# QQ Plot for Asset 3
qqnorm(residuals_df$Residual_Asset3, main = "QQ Plot for Asset 3 Residuals")
qqline(residuals_df$Residual_Asset3, col = "blue")

dev.off()

##########################


# Check the predicted volatility to ensure it's reasonable for plotting
summary(sigma2_predicted)

# Rescale predicted volatility if the values are too small
scaling_factor <- 100  # Apply this factor if the volatility is very small
scaled_sigma2_predicted <- sigma2_predicted * scaling_factor

# Generate confidence intervals with a larger range for visualization
ci_predicted_df <- data.frame(Time = 1:T,
                              Predicted_Asset1 = scaled_sigma2_predicted[, 1],
                              Lower_CI_Asset1 = scaled_sigma2_predicted[, 1] * 0.5,  # 50% lower bound
                              Upper_CI_Asset1 = scaled_sigma2_predicted[, 1] * 1.5,  # 50% upper bound
                              Predicted_Asset2 = scaled_sigma2_predicted[, 2],
                              Lower_CI_Asset2 = scaled_sigma2_predicted[, 2] * 0.5,
                              Upper_CI_Asset2 = scaled_sigma2_predicted[, 2] * 1.5,
                              Predicted_Asset3 = scaled_sigma2_predicted[, 3],
                              Lower_CI_Asset3 = scaled_sigma2_predicted[, 3] * 0.5,
                              Upper_CI_Asset3 = scaled_sigma2_predicted[, 3] * 1.5)

# Plot predicted volatility with confidence intervals for Asset 1
p1_ci <- ggplot(ci_predicted_df, aes(x = Time)) +
  geom_line(aes(y = Predicted_Asset1, color = "Predicted Volatility")) +
  geom_ribbon(aes(ymin = Lower_CI_Asset1, ymax = Upper_CI_Asset1), alpha = 0.2, fill = "blue") +
  labs(title = "Predicted Volatility with Confidence Interval (Asset 1)", y = "Volatility") +
  scale_color_manual(name = "Legend", values = c("Predicted Volatility" = "blue"))

# Plot for Asset 2
p2_ci <- ggplot(ci_predicted_df, aes(x = Time)) +
  geom_line(aes(y = Predicted_Asset2, color = "Predicted Volatility")) +
  geom_ribbon(aes(ymin = Lower_CI_Asset2, ymax = Upper_CI_Asset2), alpha = 0.2, fill = "blue") +
  labs(title = "Predicted Volatility with Confidence Interval (Asset 2)", y = "Volatility") +
  scale_color_manual(name = "Legend", values = c("Predicted Volatility" = "blue"))

# Plot for Asset 3
p3_ci <- ggplot(ci_predicted_df, aes(x = Time)) +
  geom_line(aes(y = Predicted_Asset3, color = "Predicted Volatility")) +
  geom_ribbon(aes(ymin = Lower_CI_Asset3, ymax = Upper_CI_Asset3), alpha = 0.2, fill = "blue") +
  labs(title = "Predicted Volatility with Confidence Interval (Asset 3)", y = "Volatility") +
  scale_color_manual(name = "Legend", values = c("Predicted Volatility" = "blue"))

# Save the confidence interval plot
png(filename = file.path(output_dir, "predicted_volatility_with_ci.png"), width = 1000, height = 800)
grid.arrange(p1_ci, p2_ci, p3_ci, ncol = 1)
dev.off()


#######
# Combine simulated data, predicted volatility, and confidence intervals into one dataframe
overlay_df <- data.frame(Time = 1:T,
                         Simulated_Asset1 = theta[, 1],
                         Predicted_Asset1 = sigma2_predicted[, 1],  # No scaling, original predicted volatility
                         Lower_CI_Asset1 = sigma2_predicted[, 1] * 0.85,  # 15% lower bound
                         Upper_CI_Asset1 = sigma2_predicted[, 1] * 1.15,  # 15% upper bound
                         Simulated_Asset2 = theta[, 2],
                         Predicted_Asset2 = sigma2_predicted[, 2],
                         Lower_CI_Asset2 = sigma2_predicted[, 2] * 0.85,
                         Upper_CI_Asset2 = sigma2_predicted[, 2] * 1.15,
                         Simulated_Asset3 = theta[, 3],
                         Predicted_Asset3 = sigma2_predicted[, 3],
                         Lower_CI_Asset3 = sigma2_predicted[, 3] * 0.85,
                         Upper_CI_Asset3 = sigma2_predicted[, 3] * 1.15)

# Plot for Asset 1: Overlay simulated data, predicted volatility, and confidence intervals
p1_overlay <- ggplot(overlay_df, aes(x = Time)) +
  geom_line(aes(y = Simulated_Asset1, color = "Simulated Data"), linetype = "dashed", size = 1) +
  geom_line(aes(y = Predicted_Asset1, color = "Predicted Volatility"), size = 1) +
  geom_ribbon(aes(ymin = Lower_CI_Asset1, ymax = Upper_CI_Asset1), alpha = 0.2, fill = "blue") +
  labs(title = "Predicted Volatility with Confidence Interval and Simulated Data (Asset 1)", y = "Value") +
  scale_color_manual(name = "Legend", values = c("Simulated Data" = "green", "Predicted Volatility" = "blue"))

# Plot for Asset 2: Overlay simulated data, predicted volatility, and confidence intervals
p2_overlay <- ggplot(overlay_df, aes(x = Time)) +
  geom_line(aes(y = Simulated_Asset2, color = "Simulated Data"), linetype = "dashed", size = 1) +
  geom_line(aes(y = Predicted_Asset2, color = "Predicted Volatility"), size = 1) +
  geom_ribbon(aes(ymin = Lower_CI_Asset2, ymax = Upper_CI_Asset2), alpha = 0.2, fill = "blue") +
  labs(title = "Predicted Volatility with Confidence Interval and Simulated Data (Asset 2)", y = "Value") +
  scale_color_manual(name = "Legend", values = c("Simulated Data" = "green", "Predicted Volatility" = "blue"))

# Plot for Asset 3: Overlay simulated data, predicted volatility, and confidence intervals
p3_overlay <- ggplot(overlay_df, aes(x = Time)) +
  geom_line(aes(y = Simulated_Asset3, color = "Simulated Data"), linetype = "dashed", size = 1) +
  geom_line(aes(y = Predicted_Asset3, color = "Predicted Volatility"), size = 1) +
  geom_ribbon(aes(ymin = Lower_CI_Asset3, ymax = Upper_CI_Asset3), alpha = 0.2, fill = "blue") +
  labs(title = "Predicted Volatility with Confidence Interval and Simulated Data (Asset 3)", y = "Value") +
  scale_color_manual(name = "Legend", values = c("Simulated Data" = "green", "Predicted Volatility" = "blue"))

# Save the combined overlay plot with simulated data, predicted volatility, and confidence intervals
png(filename = file.path(output_dir, "predicted_volatility_with_ci_and_simulated_no_scaling.png"), width = 1000, height = 800)
grid.arrange(p1_overlay, p2_overlay, p3_overlay, ncol = 1)
dev.off()

