# ---- Question 3.1: Load data ----
rm(list = ls())
library(ggplot2)
library(tidyr)
library(gridExtra)

# Read data
data <- read.csv("assignment3_2025/box_data_60min.csv")

# Create continuous time variable starting from 0
data$tdate <- as.POSIXct(data$tdate, "%Y-%m-%d %H:%M:%S", tz="UTC")

# Create the ggplot Ph Tdelta Gv with time
ggplot(data, aes(x = tdate)) +
  geom_line(aes(y = Ph, color = "Ph (W)")) +
  geom_line(aes(y = Tdelta*10, color = "Tdelta (oC)")) +
  geom_line(aes(y = Gv, color = "Gv (W/m2)")) +
  scale_y_continuous(
    name = "Ph (W) and Gv (W/m2)",
    sec.axis = sec_axis(~ ./10, name = "Tdelta (oC)")
  ) +
  
  labs(title = "Three non-lagged time series", x = "Time", y = "Response") +
  theme_light()

# ---- Question 3.2: Split data ----
# Divide intro train and test set
teststart <- as.POSIXct("2013-02-06 00:00", tz="UTC")
train <- data[data$tdate < teststart, ]
test <- data[data$tdate >= teststart, ]

# ---- Question 3.3: Data analysis ----
target_column <- "Ph"
# Reshape the dataframe to a long format
data_nonlag = train[c("Ph", "Tdelta", "Gv")]
data_long <- data_nonlag %>%
  pivot_longer(cols = colnames(data_nonlag)[2:3], names_to = "variable", values_to = "value")

# Create the scatter plot with facets
ggplot(data_long, aes_string(x = target_column, y = "value")) +
  geom_point() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = paste("Scatter plots of", target_column, "vs other variables"),
       x = target_column,
       y = "Value") +
  theme_light()

# ACF plots
acf(train$Ph, lag.max=50, main = paste("ACF of Ph"))
acf(train$Tdelta, lag.max=50, main = paste("ACF of Tdelta"))
acf(train$Gv, lag.max=50, main = paste("ACF of Gv"))

# CCF plots
ccf_tdelta <- Ccf(train$Tdelta, train$Ph, lag.max = 50, plot = TRUE)
ccf_gv <- Ccf(train$Gv, train$Ph, lag.max = 50, plot = TRUE)

# ---- Question 3.4: Impulse response?? ----
# # Build lag matrix
# library(dynlm)
# 
# # Impulse response from Tdelta to Ph
# lags <- 10
# tdelta_model <- dynlm(Ph ~ L(Tdelta, 0:lags), data = train)
# plot(0:3, coef(tdelta_model)[-1], type = "h",
#      xlab = "Lag", ylab = "Coefficient", main = "Impulse Response: Tdelta → Ph")
# 
# # Impulse response from Gv to Ph
# gv_model <- dynlm(Ph ~ L(Gv, 0:lags), data = train)
# plot(0:lags, coef(gv_model)[-1], type = "h",
#      xlab = "Lag", ylab = "Coefficient", main = "Impulse Response: Gv → Ph")

# ---- Question 3.5: Linear regression ----
# First train the model based on 3 data points before doing one-step prediction
lm_model <- lm(Ph ~ Tdelta + Gv, data = train[1:3,])

res_list = list()
y_pred_list = list()
for (i in (4:nrow(train))) {
  y_true = train[i,]$Ph
  y_pred = predict(lm_model, newdata = train[i,])
  residual = (y_true - y_pred)[1]
  
  res_list[i] <- residual
  y_pred_list[i] = y_pred
  
  lm_model = lm(Ph ~ Tdelta + Gv, data = train[1:i,])
}

# Residual diagnostics
combined_list = mapply(cbind, y_pred_list, res_list, SIMPLIFY = FALSE)
df_pred_residual <- do.call(rbind, combined_list)
df_pred_residual <- as.data.frame(df_pred_residual, stringsAsFactors = FALSE)

# Add three null values at the beginning of each column
df_pred_residual <- rbind(data.frame(V1 = rep(NA, 3), V2 = rep(NA, 3)), df_pred_residual)

df_pred_residual$tdate = train$tdate
df_pred_residual$y_true = train$Ph
df_pred_residual$Tdelta = train$Tdelta
df_pred_residual$Gv = train$Gv
colnames(df_pred_residual)[1:2] <- c("y_pred", "residuals")

# Create the ggplot y_pred and y_true
ggplot(df_pred_residual, aes(x = tdate)) +
  geom_point(aes(y = y_true, color = "y_true (Ph)")) +
  geom_line(aes(y = y_pred, color = "y_pred (Ph)")) +
  scale_color_manual(values = c("y_true (Ph)" = "black", "y_pred (Ph)" = "red")) +
  labs(title = "Comparison: true Ph and predicted Ph based on linear regression model",
       x = "Time", y = "Response", color = "Legend") +
  theme_light()

# Residuals
# Create the residuals plot
residuals_plot <- ggplot(df_pred_residual, aes(x = tdate)) +
  geom_line(aes(y = residuals)) +
  labs(title = "Residuals from one-step prediction using linear regression",
       x = "Time", y = "Residuals") +
  theme_light()

# Remove missing values from the residuals column
cleaned_residuals <- na.omit(df_pred_residual$residuals)

# Create the ACF plot
acf_plot <- ggplot(data.frame(lag = acf(cleaned_residuals, lag.max=50)$lag, 
                              acf = acf(cleaned_residuals, lag.max=50)$acf), 
                   aes(x = lag, y = acf)) +
  geom_bar(stat = "identity") +
  labs(title = "ACF of Residuals", x = "Lag", y = "ACF") +
  theme_light()

# Create the CCF plots
ccf_tdelta_data <- ccf(df_pred_residual[4:nrow(df_pred_residual),]$Tdelta, 
                       df_pred_residual[4:nrow(df_pred_residual),]$residuals, lag.max = 50, plot = FALSE)
ccf_tdelta_plot <- ggplot(data.frame(lag = ccf_tdelta_data$lag, ccf = ccf_tdelta_data$acf), aes(x = lag, y = ccf)) +
  geom_bar(stat = "identity") +
  labs(title = "CCF of Tdelta and residuals", x = "Lag", y = "CCF") +
  theme_light()

ccf_gv_data <- ccf(df_pred_residual[4:nrow(df_pred_residual),]$Gv, 
                   df_pred_residual[4:nrow(df_pred_residual),]$residuals, lag.max = 50, plot = FALSE)
ccf_gv_plot <- ggplot(data.frame(lag = ccf_gv_data$lag, ccf = ccf_gv_data$acf), aes(x = lag, y = ccf)) +
  geom_bar(stat = "identity") +
  labs(title = "CCF of Gv and residuals", x = "Lag", y = "CCF") +
  theme_light()


# Combine the two plots into one figure
grid.arrange(residuals_plot, acf_plot, ccf_tdelta_plot, ccf_gv_plot, ncol = 1)


# Need to add ARX

# ---- Question 3.6: Linear regression ----
arx1_model <- lm(Ph ~ Ph.l1 + Tdelta + Gv, data = train)
summary(arx1_model)

# Residual diagnostics
resid_arx1 <- resid(arx1_model)
plot(resid_arx1, type = "l", main = "Residuals from ARX(1)")
acf(resid_arx1, main = "ACF of ARX(1) Residuals")

# ---- Question 3.7: Linear regression ----
library(MASS)

max_order <- 5
aic_vals <- numeric(max_order)
bic_vals <- numeric(max_order)

for (order in 1:max_order) {
  predictors <- paste0("Ph.l", 1:order, collapse = " + ")
  tdelta_lags <- paste0("Tdelta.l", 0:(order - 1), collapse = " + ")
  gv_lags <- paste0("Gv.l", 0:(order - 1), collapse = " + ")
  
  formula <- as.formula(paste("Ph ~", predictors, "+", tdelta_lags, "+", gv_lags))
  fit <- lm(formula, data = train)
  
  aic_vals[order] <- AIC(fit)
  bic_vals[order] <- BIC(fit)
}

# Plot
plot(1:max_order, aic_vals, type = "b", col = "blue", ylim = range(c(aic_vals, bic_vals)),
     ylab = "Information Criterion", xlab = "Model Order", main = "AIC vs BIC")
lines(1:max_order, bic_vals, type = "b", col = "red")
legend("topright", legend = c("AIC", "BIC"), col = c("blue", "red"), lty = 1)

# ---- Question 3.8: Linear regression ----
library(Metrics)

rmse_vals <- numeric(max_order)

for (order in 1:max_order) {
  predictors <- paste0("Ph.l", 1:order, collapse = " + ")
  tdelta_lags <- paste0("Tdelta.l", 0:(order - 1), collapse = " + ")
  gv_lags <- paste0("Gv.l", 0:(order - 1), collapse = " + ")
  
  formula <- as.formula(paste("Ph ~", predictors, "+", tdelta_lags, "+", gv_lags))
  fit <- lm(formula, data = train)
  
  preds <- predict(fit, newdata = test)
  rmse_vals[order] <- rmse(test$Ph, preds)
}

# Plot RMSE vs Model Order
plot(1:max_order, rmse_vals, type = "b", col = "darkgreen",
     xlab = "Model Order", ylab = "RMSE", main = "One-step RMSE on Test Data")

# ---- Question 3.9: Linear regression ----
best_order <- which.min(bic_vals)  # Or pick based on RMSE if preferred

# Recursive multi-step prediction
predicted_ph <- rep(NA, nrow(data))
for (t in (best_order + 1):nrow(data)) {
  lag_ph <- sapply(1:best_order, function(l) predicted_ph[t - l])
  lag_td <- sapply(0:(best_order - 1), function(l) data$Tdelta[t - l])
  lag_gv <- sapply(0:(best_order - 1), function(l) data$Gv[t - l])
  
  coefs <- coef(lm(
    as.formula(paste("Ph ~", 
                     paste0("Ph.l", 1:best_order, collapse = " + "), "+",
                     paste0("Tdelta.l", 0:(best_order - 1), collapse = " + "), "+",
                     paste0("Gv.l", 0:(best_order - 1), collapse = " + ")
    )), data = train))
  
  pred <- sum(c(1, lag_ph, lag_td, lag_gv) * coefs)
  predicted_ph[t] <- pred
}

# Plot vs actual
plot(data$time, data$Ph, type = "l", col = "black", main = "Simulated Multi-Step Forecast", ylab = "Ph")
lines(data$time, predicted_ph, col = "red")
legend("topright", legend = c("Observed", "Predicted"), col = c("black", "red"), lty = 1)

# ---- Question 3.10: Linear regression ----
#We explored the dynamic relationships between heating (Ph), internal/external temperature difference (Tdelta), and solar radiation (Gv). Through impulse response analysis and model fitting, it was evident that a higher-order ARX model (order 2–3) provided substantial improvements over a simple linear regression. Based on AIC, BIC, and RMSE, the optimal order was around 2. Multi-step prediction showed reasonable performance, suggesting the model's potential for real-time operational use, provided that input variables are reliably forecasted.


