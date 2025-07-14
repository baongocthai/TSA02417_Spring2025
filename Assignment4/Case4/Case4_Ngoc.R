# ---- Question 1.1: Simulate 5 times ----
set.seed(123)

# Parameters
n <- 100
a <- 0.9
b <- 1
sigma1 <- 1
X0 <- 5
num_simulations <- 5

# Matrix to store simulations
simulations <- matrix(NA, nrow = n, ncol = num_simulations)

# Simulate 5 trajectories
for (i in 1:num_simulations) {
  X <- numeric(n)
  X[1] <- X0
  for (t in 2:n) {
    X[t] <- a * X[t - 1] + b + rnorm(1, mean = 0, sd = sigma1)
  }
  simulations[, i] <- X
}

# Plot all trajectories
matplot(simulations, type = "l", lty = 1, col = rainbow(num_simulations),
        main = "5 Simulated Trajectories of X_t",
        xlab = "Time", ylab = "X_t")
legend("topleft", legend = paste("Sim", 1:num_simulations),
       col = rainbow(num_simulations), lty = 1)


# ---- Question 1.2: Latent state & observations ----
# Parameters
n <- 100
a <- 0.9
b <- 1
sigma1 <- 1
sigma2 <- 1
X0 <- 5

# Initialize vectors
X <- numeric(n)
Y <- numeric(n)
X[1] <- X0
Y[1] <- X[1] + rnorm(1, mean = 0, sd = sigma2)

# Simulate the process
for (t in 2:n) {
  X[t] <- a * X[t - 1] + b + rnorm(1, mean = 0, sd = sigma1)
  Y[t] <- X[t] + rnorm(1, mean = 0, sd = sigma2)
}

# Plot the latent state and observations
plot(Y, type = "l", col = "blue", lty = 1, ylim = range(c(X, Y)),
     main = "Latent State vs Noisy Observations",
     xlab = "Time", ylab = "Value")
lines(X, col = "red", lty = 2)
legend("topright", legend = c("Observation Y_t", "True State X_t"),
       col = c("blue", "red"), lty = c(1, 2))

# check if the difference between the 2 are white noise
residual = X - Y
par(mfrow = c(1, 2)) # 1 row, 2 columns
acf(residual, main='ACF of residuals', lag.max = 30)
pacf(residual, main='PACF of residuals', lag.max = 30)

# ---- Question 1.3: Kalman filter ----
# Simulate data (from Part 1.2)
set.seed(123)
n <- 100
a <- 0.9
b <- 1
sigma1 <- 1
sigma2 <- 1
X0 <- 5
P0 <- 1
R = sigma2^2

X <- numeric(n)
Y <- numeric(n)
X[1] <- X0
Y[1] <- X[1] + rnorm(1, 0, sigma2)

for (t in 2:n) {
  X[t] <- a * X[t - 1] + b + rnorm(1, 0, sigma1)
  Y[t] <- X[t] + rnorm(1, 0, sigma2)
}

# Apply Kalman filter
theta <- c(a, b, sigma1)
kf_result <- myKalmanFilter(Y, theta, R, X0, P0)

# Plot results
par(mfrow = c(1, 1))
time <- 1:n
upper <- kf_result$x_pred + 1.96 * sqrt(kf_result$P_pred)
lower <- kf_result$x_pred - 1.96 * sqrt(kf_result$P_pred)

# Plot with shaded 95% confidence interval
plot(time, Y, type = "l", col = "blue", ylim = range(c(Y, X, upper, lower)),
     ylab = "Value", xlab = "Time", main = "Kalman Filter Prediction with 95% CI")

# Add shaded area for 95% confidence interval
polygon(c(time, rev(time)), c(upper, rev(lower)),
        col = rgb(1, 0, 0, alpha = 0.3), border = NA)

# Add lines
lines(time, X, col = "black", lty = 2, lwd = 3)
lines(time, kf_result$x_pred, col = "red")

# Add legend
legend("topleft", legend = c("Observation Y_t", "True State X_t", "Predicted State", "95% CI"),
       col = c("blue", "black", "red", rgb(1, 0, 0, alpha = 0.3)),
       lty = c(1, 2, 1, NA), pch = c(NA, NA, NA, 15), pt.cex = 2, bty = "n")

# ---- Question 1.4: Kalman filter Gaussian noise ----
# Required functions: myLogLikFun and myKalmanFilter must be defined before this
simulate_data <- function(a, b, sigma1, sigma2, n) {
  x <- numeric(n)
  y <- numeric(n)
  x[1] <- rnorm(1, 0, sqrt(10))  # Initial state
  y[1] <- x[1] + rnorm(1, 0, sigma2)
  
  for (t in 2:n) {
    x[t] <- a * x[t - 1] + b + rnorm(1, 0, sigma1)
    y[t] <- x[t] + rnorm(1, 0, sigma2)
  }
  
  return(y)
}

# Parameters
a_true <- 1
b_true <- 0.9
sigma1 <- 1
R <- 1  # sigma2^2
n <- 100
num_simulations <- 100

# Store estimates
estimates <- matrix(NA, nrow = num_simulations, ncol = 3)
colnames(estimates) <- c("a", "b", "sigma1")

set.seed(123)  # For reproducibility

for (i in 1:num_simulations) {
  y_sim <- simulate_data(a_true, b_true, sigma1, sqrt(R), n)
  
  # Initial guess for parameters
  init_theta <- c(0.5, 0.5, 0.5)
  
  # Bound sigma1 to be positive using method = "L-BFGS-B"
  result <- tryCatch({
    optim(
      par = init_theta,
      fn = myLogLikFun,
      y = y_sim,
      R = R,
      method = "L-BFGS-B",
      lower = c(-Inf, -Inf, 1e-4),
      upper = c(Inf, Inf, Inf)
    )
  }, error = function(e) {
    warning("Optimization failed: ", conditionMessage(e))
    list(par = rep(NA, 3))
  })
  print(result$par)
  
  estimates[i, ] <- result$par
  cat("Simulation", i, "done\n")
}

# Save or inspect the results
estimates_df <- as.data.frame(estimates)
summary(estimates_df)

# Optional: visualize estimates
boxplot(estimates_df, main = "Parameter Estimates (a=1, b=0.9, sigma1=1)")

# ---- Question 1.5.1: Kalman filter Student's t-distribution noise  ----
# Define degrees of freedom
nu_vals <- c(100, 5, 2, 1)
n <- 100
simulations <- 100
a <- 1
b <- 0.9
sigma1 <- 1
sigma2 <- 1

# Density plot setup
x_vals <- seq(-5, 5, length.out = 1000)
densities <- data.frame(x = x_vals)

# Add standard normal density
densities$Normal <- dnorm(x_vals, mean = 0, sd = 1)

# Add t-distribution densities
for (nu in nu_vals) {
  densities[[paste0("t_", nu)]] <- dt(x_vals, df = nu)
}

# Convert to long format for ggplot
dens_long <- pivot_longer(densities, cols = -x, names_to = "Distribution", values_to = "Density")

# Plot densities
ggplot(dens_long, aes(x = x, y = Density, color = Distribution)) +
  geom_line(size = 1) +
  labs(
    title = "Comparison of t-distribution Densities vs Standard Normal",
    x = "Value", y = "Density"
  ) +
  theme_minimal()

# ---- Question 1.5.2: Kalman filter Student's t-distribution noise (optimized) ----
simulate_data_t <- function(a, b, df1, sigma1, sigma2, n) {
  x <- numeric(n)
  y <- numeric(n)
  x[1] <- sigma1 * rt(n,df1)  # Initial state
  y[1] <- x[1] + rnorm(1, 0, sigma2)
  
  for (t in 2:n) {
    x[t] <- a * x[t - 1] + b + sigma1 * rt(n,df1)
    y[t] <- x[t] + rnorm(1, 0, sigma2)
  }
  
  return(y)
}

# Parameters
a_true <- 1
b_true <- 0.9
df <- 100
sigma1 <- 1
R <- 1  # sigma2^2
n <- 100
num_simulations <- 100

# Store estimates
estimates <- matrix(NA, nrow = num_simulations, ncol = 3)
colnames(estimates) <- c("a", "b", "sigma1")

set.seed(123)  # For reproducibility

for (i in 1:num_simulations) {
  y_sim <- simulate_data_t(a_true, b_true, df, sigma1, sqrt(R), n)
  
  # Initial guess for parameters
  init_theta <- c(0.5, 0.5, 0.5)
  
  # Bound sigma1 to be positive using method = "L-BFGS-B"
  result <- optim(
    par = init_theta,
    fn = myLogLikFun,
    y = y_sim,
    R = R,
    method = "L-BFGS-B",
    lower = c(-Inf, -Inf, 1e-4),  # avoid sigma1 = 0
    upper = c(Inf, Inf, Inf)
  )
  
  estimates[i, ] <- result$par
  cat("Simulation", i, "done\n")
}

# Save or inspect the results
estimates_df <- as.data.frame(estimates)
summary(estimates_df)

# Optional: visualize estimates
boxplot(estimates_df, main = "Parameter Estimates (a=1, b=0.9, sigma1=1, df=1)")
