# Load necessary packages
rm(list = ls())
library(forecast)
library(ggplot2)
library(tseries)

# 2.1 Read data and prepare log-transformed series
data <- read.csv("assignment3_2025/datasolar.csv")
Yt <- data$power
mu <- 5.72
Xt <- log(Yt) - mu

# Plot the original series
plot.ts(Yt, main="Monthly Solar Power Generation (MWh)", ylab="Generation (MWh)", xlab="Time")

# Fit the seasonal AR(1)(12) model manually
# (1 + φ1*B)(1 + Φ1*B^12)(Xt) = εt
phi1 <- -0.38
Phi1 <- -0.94
n <- length(Xt)
resid <- rep(NA, n)

for (t in 14:n) {
  resid[t] <- Xt[t] + phi1 * Xt[t - 1] + Phi1 * Xt[t - 12] + phi1 * Phi1 * Xt[t - 13]
}

# Plot residuals
resid_clean <- na.omit(resid)

# Set up the plotting area to have 3 rows and 1 column
par(mfrow = c(3, 1))

# Plot the residuals time series
ts.plot(resid_clean, main = "Residuals from AR(1)(12) Model", ylab = "Residuals")

# Plot the ACF of residuals
acf(resid_clean, main = "ACF of Residuals")

# Plot the PACF of residuals
pacf(resid_clean, main = "PACF of Residuals")

# Normality check
qqnorm(resid_clean)
qqline(resid_clean)
print(shapiro.test(resid_clean))

# 2.2 Forecast the next 12 months
h <- 12
Xt_forecast <- numeric(h)
Xt_forecast[1] <- -phi1 * Xt[n] - Phi1 * Xt[n - 11] - phi1 * Phi1 * Xt[n - 12]

# resid[t] <- Xt[t] + phi1 * Xt[t - 1] + Phi1 * Xt[t - 12] + phi1 * Phi1 * Xt[t - 13]

for (k in 2:h) {
  x_lag1 <- ifelse(k == 1, Xt[n], Xt_forecast[k - 1])
  x_lag12 <- ifelse(k <= 12, Xt[n - 12 + k], Xt_forecast[k - 12])
  x_lag13 <- ifelse(k <= 13, Xt[n - 13 + k], Xt_forecast[k - 13])
  
  Xt_forecast[k] <- -phi1 * x_lag1 - Phi1 * x_lag12 - phi1 * Phi1 * x_lag13
}

Yt_forecast <- exp(Xt_forecast + mu)

par(mfrow = c(1, 1))
# 2.3 95% prediction intervals (using AR(1) only)
sigma2 <- 0.222
se <- sqrt(sigma2 * cumsum((1 + phi1^(1:h))^2))
lower <- exp(Xt_forecast + mu - 1.96 * se)
upper <- exp(Xt_forecast + mu + 1.96 * se)

# Plot forecast
# ts_all <- ts(c(Yt, Yt_forecast), start=1, frequency=12)
# time_all <- time(ts_all)
# plot(time_all, ts_all, type="l", ylim=c(min(Yt), max(upper)), main="Forecast with 95% Prediction Interval",
#      xlab="Time", ylab="Power (MWh)", col="black")
# lines(time_all[(length(Yt) + 1):(length(Yt) + h)], lower, col="red", lty=2)
# lines(time_all[(length(Yt) + 1):(length(Yt) + h)], upper, col="red", lty=2)

# Assuming Yt_forecast, lower, and upper are vectors with the same length
time <- seq_along(Yt_forecast)

# Plot the forecast
plot(time, Yt_forecast, type = "l", col = "blue", ylim = range(c(lower, upper, Yt_forecast)),
     xlab = "Time", ylab = "Value", main = "Forecast with 95% Confidence Interval")

# Add the shaded area for the confidence interval
polygon(c(time, rev(time)), c(lower, rev(upper)), col = rgb(0, 0, 1, 0.2), border = NA)

# Add the forecast line again to be on top of the shaded area
lines(time, Yt_forecast, col = "blue")


# Tabulower# Tabupper.tri()# Table of forecast values
forecast_table <- data.frame(
  Month = 1:h,
  Forecast = round(Yt_forecast, 2),
  Lower_95 = round(lower, 2),
  Upper_95 = round(upper, 2)
)

print(forecast_table)
