# Load necessary library
library(forecast)

# Function to simulate and plot seasonal ARIMA models
simulate_seasonal_arima <- function(model, title) {
  set.seed(123)
  sim_data <- arima.sim(model = model, n = 200)
  
  # Set layout for plots
  layout(rbind(1,2:3))
  par(mar=c(3,3,1,1), mgp=c(2, 0.7,0))
  
  # Plot time series
  plot(sim_data, main = paste("Simulated", title), ylab = "Value", col = "blue")
  
  # Set layout for ACF and PACF plots
  acf(sim_data, lag.max=50, main = paste("ACF of", title))
  pacf(sim_data, lag.max=50, main = paste("PACF of", title))
}

# 2.1: (1,0,0) × (0,0,0)[12] with phi1 = 0.6
simulate_seasonal_arima(list(order = c(1,0,0), seasonal = list(order = c(0,0,0), period = 12), ar = 0.6), "Q2.1: (1,0,0) × (0,0,0)[12]")

# 2.2: (0,0,0) × (1,0,0)[12] with Phi1 (sar) = -0.9
simulate_seasonal_arima(list(order = c(0,0,0), seasonal = list(order = c(1,0,0), period = 12), sar = -0.9), "Q2.2: (0,0,0) × (1,0,0)[12]")

# 2.3: (1,0,0) × (0,0,1)[12] with phi1 = 0.9, Theta1 (sma) = -0.7
simulate_seasonal_arima(list(order = c(1,0,0), seasonal = list(order = c(0,0,1), period = 12), ar = 0.9, sma = -0.7), "Q2.3: (1,0,0) × (0,0,1)[12]")

# 2.4: (1,0,0) × (1,0,0)[12] with phi1 = -0.6, Phi1 (sar) = -0.8
simulate_seasonal_arima(list(order = c(1,0,0), seasonal = list(order = c(1,0,0), period = 12), ar = -0.6, sar = -0.8), "Q2.4: (1,0,0) × (1,0,0)[12]")

# 2.5: (0,0,1) × (0,0,1)[12] with theta1 = 0.4, Theta1 (sma) = -0.8
simulate_seasonal_arima(list(order = c(0,0,1), seasonal = list(order = c(0,0,1), period = 12), ma = 0.4, sma = -0.8), "Q2.5: (0,0,1) × (0,0,1)[12]")

# 2.6: (0,0,1) × (1,0,0)[12] with theta1 = -0.4, Phi1 (sar) = 0.7
simulate_seasonal_arima(list(order = c(0,0,1), seasonal = list(order = c(1,0,0), period = 12), ma = -0.4, sar = 0.7), "Q2.6: (0,0,1) × (1,0,0)[12]")
