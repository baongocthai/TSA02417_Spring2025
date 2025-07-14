# ---- Question 3.1: Load data ----
rm(list = ls())
library(ggplot2)
library(tidyr)
library(gridExtra)
library(forecast)

# Read data
data <- read.csv("assignment3_2025/box_data_60min.csv")

# Create continuous time variable starting from 0
data$tdate <- as.POSIXct(data$tdate, "%Y-%m-%d %H:%M:%S", tz="UTC")

# Create the ggplot Ph Tdelta Gv with time
ggplot(data, aes(x = tdate)) +
  geom_line(aes(y = Ph, color = "Ph (W)"), size = 1) +
  geom_line(aes(y = Tdelta*10, color = "Tdelta (oC)"), size = 1) +
  geom_line(aes(y = Gv, color = "Gv (W/m2)"), size = 1) +
  scale_y_continuous(
    name = "Ph (W) and Gv (W/m2)",
    sec.axis = sec_axis(~ ./10, name = "Tdelta (oC)")
  ) +
  
  labs(title = "Three non-lagged time series", x = "Time", y = "Response") +
  theme_light()

# ---- Question 3.2: Split data ----
# Divide intro train and test set
teststart <- as.POSIXct("2013-02-06 00:00", tz="UTC")
train <- data[data$tdate <= teststart, ]
test <- data[data$tdate > teststart, ]

# ---- Question 3.3: Data analysis ----
target_column <- "Ph"
# Reshape the dataframe to a long format
data_nonlag = train[c("Ph", "Tdelta", "Gv")]
data_long <- data_nonlag %>%
  pivot_longer(cols = colnames(data_nonlag)[2:3], names_to = "variable", values_to = "value")

# Create the scatter plot with facets
# First plot: Scatter plots with ggplot2
scatter_plot <- ggplot(data_long, aes_string(x = target_column, y = "value")) +
  geom_point() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(
    title = paste("Scatter plots of", target_column, "vs other variables"),
    x = target_column,
    y = "Value"
  ) +
  theme_light()

# Second plot: ACF plots using base R
acf_ph <- acf(train$Ph, lag.max = 50, plot = FALSE)
acf_tdelta <- acf(train$Tdelta, lag.max = 50, plot = FALSE)
acf_gv <- acf(train$Gv, lag.max = 50, plot = FALSE)

# Convert ACF plots to ggplot-like objects using autoplot from forecast package
acf_ph_plot <- autoplot(acf_ph) + labs(title = "ACF of Ph")
acf_tdelta_plot <- autoplot(acf_tdelta) + labs(title = "ACF of Tdelta")
acf_gv_plot <- autoplot(acf_gv) + labs(title = "ACF of Gv")

# Third plot: CCF plots using base R
ccf_tdelta <- Ccf(train$Tdelta, train$Ph, lag.max = 50, plot = FALSE)
ccf_gv <- Ccf(train$Gv, train$Ph, lag.max = 50, plot = FALSE)

# Convert CCF plots to ggplot-like objects using autoplot from forecast package
ccf_tdelta_plot <- autoplot(ccf_tdelta) + labs(title = "CCF of Tdelta vs Ph")
ccf_gv_plot <- autoplot(ccf_gv) + labs(title = "CCF of Gv vs Ph")

# Combine all plots into one figure
combined_plot <- grid.arrange(
  scatter_plot, acf_ph_plot, acf_tdelta_plot, acf_gv_plot,
  ccf_tdelta_plot, ccf_gv_plot,
  ncol = 2, nrow = 3
)

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

# Overall model: no one-step predictions
lm_model_overall <- lm(Ph ~ Tdelta + Gv, data = train)
y_pred_overall = predict(lm_model_overall, newdata = train)

df_pred_residual$y_pred_overall = y_pred_overall

# Plot for y_pred_overall
pred_overall = ggplot(df_pred_residual, aes(x = tdate)) +
  geom_point(aes(y = y_true, color = "y_true (Ph)")) +
  geom_line(aes(y = y_pred_overall, color = "y_pred (Ph)")) +
  scale_color_manual(values = c("y_true (Ph)" = "black", "y_pred (Ph)" = "red")) +
  labs(title = "True and predicted Ph: linear regression model",
       x = "Time", y = "Response", color = "Legend") +
  theme_light()

# Create the ggplot y_pred and y_true
pred = ggplot(df_pred_residual, aes(x = tdate)) +
  geom_point(aes(y = y_true, color = "y_true (Ph)")) +
  geom_line(aes(y = y_pred, color = "y_pred (Ph)")) +
  scale_color_manual(values = c("y_true (Ph)" = "black", "y_pred (Ph)" = "red")) +
  labs(title = "True and predicted Ph: one-step prediction, linear regression model",
       x = "Time", y = "Response", color = "Legend") +
  theme_light()

# Residuals
# Create the residuals plot
residuals_plot <- ggplot(df_pred_residual, aes(x = tdate)) +
  geom_line(aes(y = residuals)) +
  labs(title = "Residuals from one-step prediction",
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

# Combine the plots
combined_plot <- grid.arrange(
  pred_overall,                                # First plot (spans the entire row)
  pred,
  arrangeGrob(
    residuals_plot, acf_plot,          # Subsequent 4 plots (2x2 arrangement)
    ccf_tdelta_plot, ccf_gv_plot,
    ncol = 2                           # Set 2 columns for the 2x2 layout
  ),
  ncol = 1                             # Main layout: 1 column with 2 rows
)

# Need to add ARX

# ---- Question 3.6: ARX(1) ----
arx1_model <- lm(Ph ~ Ph.l1 + Tdelta + Gv, data = train[1:3,])

res_list = list()
y_pred_list = list()
for (i in (4:nrow(train))) {
  y_true = train[i,]$Ph
  y_pred = predict(arx1_model, newdata = train[i,])
  residual = (y_true - y_pred)[1]
  
  res_list[i] <- residual
  y_pred_list[i] = y_pred
  
  arx1_model <- lm(Ph ~ Ph.l1 + Tdelta + Gv, data = train[1:i,])
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

# Remove missing values from the residuals column
cleaned_residuals <- na.omit(df_pred_residual$residuals)

# Overall model: no one-step predictions
lm_model_overall <- lm(Ph ~ Ph.l1 + Tdelta + Gv, data = train)
y_pred_overall = predict(lm_model_overall, newdata = train)

df_pred_residual$y_pred_overall = y_pred_overall

# Plot for y_pred_overall
pred_overall = ggplot(df_pred_residual, aes(x = tdate)) +
  geom_point(aes(y = y_true, color = "y_true (Ph)")) +
  geom_line(aes(y = y_pred_overall, color = "y_pred (Ph)")) +
  scale_color_manual(values = c("y_true (Ph)" = "black", "y_pred (Ph)" = "red")) +
  labs(title = "True and predicted Ph: ARX model (first order)",
       x = "Time", y = "Response", color = "Legend") +
  theme_light()

# Create the ggplot y_pred and y_true
pred <- ggplot(df_pred_residual, aes(x = tdate)) +
  geom_point(aes(y = y_true, color = "y_true (Ph)")) +
  geom_line(aes(y = y_pred, color = "y_pred (Ph)")) +
  scale_color_manual(values = c("y_true (Ph)" = "black", "y_pred (Ph)" = "red")) +
  labs(title = "True and predicted Ph: one-step prediction, ARX model (first order)",
       x = "Time", y = "Response", color = "Legend") +
  theme_light()

# Residuals
# Create the residuals plot
residuals_plot <- ggplot(df_pred_residual, aes(x = tdate)) +
  geom_line(aes(y = residuals)) +
  labs(title = "Residuals from one-step prediction",
       x = "Time", y = "Residuals") +
  theme_light()

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

# Combine the plots
combined_plot <- grid.arrange(
  pred_overall,
  pred,                                # First plot (spans the entire row)
  arrangeGrob(
    residuals_plot, acf_plot,          # Subsequent 4 plots (2x2 arrangement)
    ccf_tdelta_plot, ccf_gv_plot,
    ncol = 2                           # Set 2 columns for the 2x2 layout
  ),
  ncol = 1                             # Main layout: 1 column with 2 rows
)

# ---- Question 3.7: ARX(increase order), one-step predictions (dont really make sense) ----
max_order <- 10
aic_vals <- numeric(max_order)
bic_vals <- numeric(max_order)
res_order = list()
y_pred_order = list()

## Function to do iterative fitting for any model order
rolling_lm <- function(formula, train){
  arx_model <- lm(formula, data = train[1:3,])

  res_list = list()
  y_pred_list = list()
  for (i in (4:nrow(train))) {
    y_true = train[i,]$Ph
    y_pred = predict(arx_model, newdata = train[i,])
    residual = (y_true - y_pred)[1]
    
    res_list[i] <- residual
    y_pred_list[i] = y_pred
    
    arx_model <- lm(formula, data = train[1:i,])
  }
  num_params = length(coef(arx_model))
  # Return residuals and predictions as a list
  return(list(residuals = res_list, predictions = y_pred_list, num_params = num_params))
}

## Function to calculate AIC & BIC from residuals
calculate_AIC_BIC <- function(residuals,num_params) {
  # Extract RSS (Residual Sum of Squares)
  # residuals <- model$residuals
  RSS <- sum(residuals^2)
  
  # Extract number of parameters (including intercept)
  # num_params <- length(coef(model))
  
  # Extract number of observations
  num_obs <- length(residuals)
  
  # Calculate AIC and BIC
  AIC <- 2 * num_params + num_obs * log(RSS / num_obs)
  BIC <- num_params * log(num_obs) + num_obs * log(RSS / num_obs)
  
  return(list(AIC = AIC, BIC = BIC))
}

# calculate AIC and BIC for models with different orders
for (order in 1:max_order) {
  # Create model with different orders
  predictors <- paste0("Ph.l", 1:order, collapse = " + ")
  tdelta_lags <- paste0("Tdelta.l", 0:(order - 1), collapse = " + ")
  gv_lags <- paste0("Gv.l", 0:(order - 1), collapse = " + ")
  
  # fit model one step prediction
  formula <- as.formula(paste("Ph ~", predictors, "+", tdelta_lags, "+", gv_lags))
  fit <- rolling_lm(formula, train)
  
  # calculate AIC, BIC
  y_pred = unlist(fit$predictions)
  residuals = unlist(fit$residuals)
  num_params = fit$num_params
  AIC_BIC = calculate_AIC_BIC(residuals,num_params)
  
  # save residuals and predictions for each order
  res_order[order] <- residuals
  y_pred_order[order] = y_pred
  
  aic_vals[order] <- AIC_BIC$AIC
  bic_vals[order] <- AIC_BIC$BIC
  print(order)
}

# Plot
plot(1:max_order, aic_vals, type = "b", col = "blue", ylim = range(c(aic_vals, bic_vals)),
     ylab = "Information Criterion", xlab = "Model Order", main = "AIC vs BIC")
lines(1:max_order, bic_vals, type = "b", col = "red")
legend("topleft", legend = c("AIC", "BIC"), col = c("blue", "red"), lty = 1)

# ---- Question 3.7: ARX (increase order, no one-step predictions) ----
max_order <- 10
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
legend("bottomleft", legend = c("AIC", "BIC"), col = c("blue", "red"), lty = 1)

# ---- Question 3.8: Test period one step predictions ----
library(Metrics)
max_order <- 10
rmse_test <- numeric(max_order)
y_true = test$Ph
y_pred_order = data.frame(matrix(nrow = nrow(test), ncol = 0))
res_order = data.frame(matrix(nrow = nrow(test), ncol = 0))

for (order in 1:max_order) {
  predictors <- paste0("Ph.l", 1:order, collapse = " + ")
  tdelta_lags <- paste0("Tdelta.l", 0:(order - 1), collapse = " + ")
  gv_lags <- paste0("Gv.l", 0:(order - 1), collapse = " + ")
  
  formula <- as.formula(paste("Ph ~", predictors, "+", tdelta_lags, "+", gv_lags))
  fit <- lm(formula, data = train)
  
  pred_list = list()
  #One step prediction
  for (i in (1:nrow(test))) {
    preds <- predict(fit, newdata = test[i,])
    pred_list[i] = preds
    
    # retrain model with new data
    data_new = rbind(train, test[1:i,])
    fit <- lm(formula, data = data_new)
  }
  pred_list = unlist(pred_list)
  y_pred_order[[paste0("order", order)]] <- pred_list
  
  residuals = y_true - pred_list
  res_order[[paste0("order", order)]] <- residuals
  
  rmse_vals[order] <- rmse(y_true, pred_list)
  print(order)
}

# Plot RMSE vs Model Order
plot(1:max_order, rmse_vals, type = "b", col = "darkgreen",
     xlab = "Model Order", ylab = "RMSE", main = "One-step RMSE on Test Data")

# plot predictions vs. true Ph
y_pred_order$tdate = test$tdate
y_pred_order$Observed = test$Ph

library(reshape2)

# Reshape the dataframe to long format for plotting
df_long <- melt(y_pred_order, id.vars = "tdate", variable.name = "Series", value.name = "Value")
# Reshape the dataframe to long format for plotting
res_order$tdate = test$tdate
df_long_res <- melt(res_order, id.vars = "tdate", variable.name = "Series", value.name = "Value")

# Plot all predicted series and add the observed series as a dotted black line
pred_plot = ggplot(df_long, aes(x = tdate, y = Value, color = Series)) +
  geom_line(data = df_long[df_long$Series != "Observed", ], aes(group = Series), alpha = 1) +  # Plot predicted series
  geom_point(data = df_long[df_long$Series == "Observed", ],
            color = "black", size = 1.5) +  # Overlay observed series
  labs(
    title = "Time Series Plot: Predicted vs Observed",
    x = "Time",
    y = "Value",
    color = "Series"
  ) +
  theme_light()

# Plot all residuals
residual_plot = ggplot(df_long_res, aes(x = tdate, y = Value, color = Series)) +
  geom_line(data = df_long_res, aes(group = Series), alpha = 1) +  # Plot predicted series
  labs(
    title = "Time Series Plot: Predicted vs Observed",
    x = "Time",
    y = "Value",
    color = "Series"
  ) +
  theme_light()

combined_plot = grid.arrange(pred_plot,residual_plot)

# ---- Question 3.9: Best order multi-step predictions ----
best_order <- 3  # Or pick based on RMSE if preferred

# Recursive multi-step prediction
predicted_ph <- rep(0, nrow(data))
for (t in (best_order + 1):nrow(data)) {
  # lag predictions of Ph
  lag_ph <- sapply(1:best_order, function(l) predicted_ph[t - l])
  
  # lag observed inputs
  lag_td <- sapply(0:(best_order - 1), function(l) data$Tdelta[t - l])
  lag_gv <- sapply(0:(best_order - 1), function(l) data$Gv[t - l])
  
  coefs <- coef(lm(
    as.formula(paste("Ph ~", 
                     paste0("Ph.l", 1:best_order, collapse = " + "), "+",
                     paste0("Tdelta.l", 0:(best_order - 1), collapse = " + "), "+",
                     paste0("Gv.l", 0:(best_order - 1), collapse = " + ")
    )), data = data))
  
  pred <- sum(c(1, lag_ph, lag_td, lag_gv) * coefs)
  predicted_ph[t] <- pred
}

# Plot vs actual
plot(data$tdate, data$Ph, type = "l", col = "black", main = "Simulated Multi-Step Forecast", ylab = "Ph")
lines(data$tdate, predicted_ph, col = "red")
legend("topleft", legend = c("Observed", "Predicted"), col = c("black", "red"), lty = 1)

# ggplot
ggplot(data, aes(x = tdate)) +
  geom_point(aes(y = Ph, color = "Observed")) +
  geom_line(aes(y = predicted_ph, color = "Simulated")) +
  labs(title = "Simulated Multi-Step Forecast", x = "Date", y = "Ph") +
  theme_minimal() +
  scale_color_manual(values = c("Observed" = "black", "Simulated" = "red")) +
  theme(legend.position = "topleft") +
  guides(color = guide_legend(override.aes = list(linetype = c(NA, 1), shape = c(16, NA))))


