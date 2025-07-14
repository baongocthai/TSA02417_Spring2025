## ---------------------------------------------------------------------------------------
# Function for plot
plotit <- function(x, title){
  layout(rbind(1,2:3))
  par(mar=c(3,3,1,1), mgp=c(2, 0.7,0))
  plot(x, ylab="X", main = title)
  acf(x, lag.max=50, lwd=2)
  pacf(x, lag.max=50, lwd=2)
}

# Load necessary library
library(forecast)
# Set seed for reproducibility
set.seed(123)
## ---------------------------------------------------------------------------------------
# Define ARIMA parameters 2.1
ar_order <- 1  # AR order
i_order <- 0   # I order (differencing)
ma_order <- 0  # MA order
seasonal_order <- c(0, 0, 0)  # Seasonal ARIMA order (P, D, Q)
seasonal_period <- 12  # Seasonal period (e.g., monthly data)
phi = 0.6 #phi
n = 300

# Generate ARIMA model with seasonal component
arima_model <- arima.sim(model = list(order = c(ar_order, i_order, ma_order), 
                                      seasonal = list(order = seasonal_order, period = seasonal_period), 
                                      ar = c(phi)), n = n)

plotit(arima_model, "A (1, 0, 0) × (0, 0, 0)12 model with the parameter ϕ1 = 0.6")

## ---------------------------------------------------------------------------------------
# Define ARIMA parameters 2.2
set.seed(123)
ar_order <- 0  # AR order
i_order <- 0   # I order (differencing)
ma_order <- 0  # MA order
seasonal_order <- c(1, 0, 0)  # Seasonal ARIMA order (P, D, Q)
seasonal_period <- 12  # Seasonal period (e.g., monthly data)
phi = -0.9 #phi
n = 300

# Generate ARIMA model with seasonal component
arima_model <- arima.sim(model = list(order = c(ar_order, i_order, ma_order), seasonal = list(order = seasonal_order, period = 12), 
                                     sar = phi), n=n)

plotit(arima_model, "A (0, 0, 0) × (1, 0, 0)12 model with the parameter Φ1 = −0.9")



