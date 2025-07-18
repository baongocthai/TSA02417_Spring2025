myKalmanFilter <- function(
    y,             # Vector of observations y_t
    theta,         # Model parameters for X_{t+1} = a*X_t + b + e_t
    R,             # Measurement noise variance
    x_prior = 0,   # Initial prior mean for X_0
    P_prior = 10   # Initial prior variance for X_0
) {
  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  N <- length(y)
  x_pred  <- numeric(N)  # Predicted means
  P_pred  <- numeric(N)  # Predicted variances
  x_filt  <- numeric(N)  # Filtered means
  P_filt  <- numeric(N)  # Filtered variances
  innovation     <- numeric(N)  # Pre-fit residuals: y[t] - x_pred[t]
  innovation_var <- numeric(N)  # Innovation covariance: P_pred[t] + R
  
  # Placeholder for Kalman gain (calculated inside loop)
  K_t <- numeric(N)
  
  for (t in seq_len(N)) {
    # the prediction step
    if (t == 1) {
      x_pred[t] <- x_prior  # the mean prediction using the prior
      P_pred[t] <- P_prior  # the variance prediction using the prior
    } else {
      x_pred[t] <- a * x_filt[t - 1] + b # the mean prediction using the previous filtered estimate
      P_pred[t] <- a^2 * P_filt[t - 1] + sigma1^2 # the variance prediction using the previous filtered estimate
    }
    
    # Update step
    innovation[t] <- y[t] - x_pred[t] # the prediction error
    innovation_var[t] <- P_pred[t] + R # the prediction error variance
    K_t[t] <- P_pred[t] / innovation_var[t] # the Kalman gain
    
    x_filt[t] <- x_pred[t] + K_t[t] * innovation[t] # the filtered estimate
    P_filt[t] <- (1 - K_t[t]) * P_pred[t] # the filtered estimate variance
  }
  
  return(list(
    x_pred = x_pred,
    P_pred = P_pred,
    x_filt = x_filt,
    P_filt = P_filt,
    innovation = innovation,
    innovation_var = innovation_var
  ))
}
