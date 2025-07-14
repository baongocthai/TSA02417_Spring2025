myLogLikFun <- function(theta, y, R, x_prior = 0, P_prior = 10) {
  # Parameters
  a <- theta[1]     # State transition coefficient
  b <- theta[2]     # Bias
  sigma1 <- theta[3]  # Standard deviation of process noise
  
  # Run Kalman filter (user-defined function for 1D)
  kf_result <- myKalmanFilter(y, theta = c(a, b, sigma1), R = R,
                              x_prior = x_prior, P_prior = P_prior)
  
  # Extract innovation and innovation variance
  err <- kf_result$innovation         # vector of innovations (residuals)
  S <- kf_result$innovation_var       # vector of innovation variances
  
  # Compute total log-likelihood
  # logL_vec <- -0.5 * (log(2 * pi * S) + (err^2) / S)
  # logL <- sum(logL_vec)
  logL <- -0.5 * sum(log(2 * pi * S) + (err^2) / S)
  
  return(-logL)  # Return negative log-likelihood for optimization
}
