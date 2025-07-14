#####################################################
################## Assignment 2  ####################
#####################################################

# Part 2
library(forecast)  

simulate_seasonal_arma <- function(p = 0, d = 0, q = 0, 
                                   P = 0, D = 0, Q = 0, 
                                   s = 12, n = 300, 
                                   ar_coefs = NULL, ma_coefs = NULL, 
                                   sar_coefs = NULL, sma_coefs = NULL, 
                                   seed = 123) {
  
  set.seed(seed)  
  
  # Define ARMA model coefficients if provided
  ar <- if (!is.null(ar_coefs)) ar_coefs else rep(0, p)
  ma <- if (!is.null(ma_coefs)) ma_coefs else rep(0, q)
  sar <- if (!is.null(sar_coefs)) sar_coefs else rep(0, P)
  sma <- if (!is.null(sma_coefs)) sma_coefs else rep(0, Q)
  
  # Simulate the ARMA process
  model <- list(order = c(p, d, q), 
                seasonal = list(order = c(P, D, Q), period = s))
  
  data_sim <- arima.sim(model = list(ar = c(ar, sar), ma = c(ma, sma)), n = n)
  
  # Plot time series
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE))
  plot(data_sim, main = paste("Simulated Seasonal ARMA(", p, ",", d, ",", q, ") × (", 
                              P, ",", D, ",", Q, ")_", s, sep=""), 
       col = "blue", type = "l", ylab = "Value", xlab = "Time")
  
  # ACF plot
  acf(data_sim, main = "Autocorrelation Function (ACF)")
  
  # PACF plot
  pacf(data_sim, main = "Partial Autocorrelation Function (PACF)")
  
  layout(1)  
  par(mfrow = c(1, 1))
  return(data_sim)
}


# A (1, 0, 0) × (0, 0, 0)12 model with the parameter ϕ1 = 0.6.
simulate_seasonal_arma(p=1, ar_coefs = c(0.6))

# A (0, 0, 0) × (1, 0, 0)12 model with the parameter Φ1 = −0.9.
simulate_seasonal_arma(P=1, sar_coefs = c(-0.9))

