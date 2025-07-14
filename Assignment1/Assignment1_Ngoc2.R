
##############################################
### Examples for Time Series Aanalysis 02417, Lecture 3 (2025)
# 1 - Linear regression with matrix notation 
# 2 - Weighted Least Squares
# 3 - From Global to Local model (OLS with fewer datapoints)
# 4 - WLS with "lambda weights" (exponential weights)
# 5 - RLS with forgetting
##############################################


# load some packages:
library(fpp2)
library(dplyr)
library(tidyverse)

# ---- Question 1.1: Load data ----
setwd("C:/Users/ngoc2/OneDrive - Danmarks Tekniske Universitet/Optagelser/2_Time Series Analysis/Assignment/Assignment1")
### Read training data
D <- read.csv("DST_BIL54.csv")

# Create continuous time variable starting from 0
D$time <- as.POSIXct(paste0(D$time,"-01"), "%Y-%m-%d", tz="UTC")

## Year to month for each of them
D$year <- 1900 + as.POSIXlt(D$time)$year + as.POSIXlt(D$time)$mon / 12 
D$year_norm = D$year - min(D$year)

## Make a copy of the data total
D$total_mil <- D$total

## Make the output variable a floating point (i.e.\ decimal number)
D$total <- as.numeric(D$total) / 1E6

## Divide intro train and test set
teststart <- as.POSIXct("2024-01-01", tz="UTC")
Dtrain <- D[D$time < teststart, ]
Dtest <- D[D$time >= teststart, ]

# Now we split the data into test and train:
train <- Dtrain
test <- Dtest

# ---- Question 1.2: Plot ----
### Plot train and test data
ggplot(train, aes(x=year, y=total)) +
  geom_point(col="black") + 
  geom_point(data=test, col="blue")+
  labs(title = "Observations",
       x = "Year",
       y = "Number of motor driven vehicles in Denmark (Million)",
       color = "Legend") +
  theme_minimal()

# ---- Question 2.2: OLS model: Linear regression with matrix notation----
# we will fit a linear model: Y_i = beta_0 + beta_1 * time_i + epsilon_i
# so we have two parameters: beta_0 and beta_1
p <- 2

# we also save the number of observations (26 obs. in the training data):
n <- length(train$year)

# X is the "design matrix"
X <- cbind(1, train$year_norm)
print(X)

# y is vector with observations:
y <- cbind(train$total)
print(y)

# to estimate parameters we solve the "normal equations":
theta_OLS <- solve( t(X) %*% X ) %*% t(X) %*% y
print(theta_OLS)
# these are the parameter estimates!
theta_0 <- theta_OLS[1]
theta_1 <- theta_OLS[2]

# now we compute y_hat values (so thet we can plot the regression line) 
yhat_ols <- X %*% theta_OLS

# plot:
ggplot(train, aes(x=year, y=total)) +
  geom_point(col="black") + 
  geom_line(aes(y=yhat_ols), col="blue", size=.5) + 
  theme_classic()

# we will now calculate the standard errors on the parameters beta_0 and beta_1:

# first compute residuals:
e_ols <- y - yhat_ols

# calculate sum of squared residuals:
RSS_ols <- t(e_ols) %*% e_ols

# calculate sigma^2:
sigma2_ols <- as.numeric(RSS_ols/(n - p))

# calculate variance-covariance matrix of _parameters_:
V_ols <- sigma2_ols * solve(t(X) %*% X)
print(V_ols)

# the variances of the parameters are the values in the diagonal:
diag(V_ols)
# and the standard errors are given by:
sqrt(diag(V_ols))

se_theta_0 <- (sqrt(diag(V_ols)))[1] # standard error of the intercept-parameter
se_theta_1 <- (sqrt(diag(V_ols)))[2] # standard error of the slope-parameter

# now we have both point estimates and standard errors:
# intercept:
theta_0
se_theta_0
# slope:
theta_1
se_theta_1

# ---- Question 2.3: OLS model prediction ----
# Predictions for future values ("forecast")
# now we use the model for predictions on future timepoints
# we use the timepoints from the testdata:
Xtest <- cbind(1, test$year_norm)
print(Xtest)

# compute predictions (we compute all 10 predictions at once):
y_pred <- Xtest %*% theta_OLS
print(y_pred)

# compute prediction variance-covariance matrix:
Vmatrix_pred <- sigma2_ols * (1+ (Xtest %*% solve( t(X) %*% X )) %*% t(Xtest))
# the variances of individual predictions are in the diagonal of the matrix above
print(diag(Vmatrix_pred))

# compute "prediction intervals" 
y_pred_lwr <- y_pred - qt(0.975, df=n-p) * sqrt(diag(Vmatrix_pred))
y_pred_upr <- y_pred + qt(0.975, df=n-p) * sqrt(diag(Vmatrix_pred))

# plot forecast:
ggplot(train, aes(x = year, y = total)) +
  geom_point(aes(color = "Train set"), size = 1) + 
  geom_point(data = test, aes(x = year, y = total, color = "Test set")) +
  geom_line(aes(y = yhat_ols, color = "OLS model"), size = 0.5) +
  geom_point(data = test, aes(x = year, y = y_pred, color = "Forecast"), size = 1) +
  geom_ribbon(data = test, aes(x = year, ymin = y_pred_lwr, ymax = y_pred_upr), inherit.aes = FALSE, alpha = 0.2, fill = "red") +
  labs(
    # title = "Forecast Plot",
    x = "Year",
    y = "Number of motor driven vehicles in Denmark (Million)",
    color = "Legend"
  ) +
  scale_color_manual(
    values = c("Train set" = "black", "Test set" = "darkorange", "OLS model" = "blue", "Forecast" = "red")
  ) +
  theme_classic()

# qq plot of residuals:
qqnorm(e_ols)
qqline(e_ols)

# plot residuals versus x (year):
ggplot(train, aes(x = year)) +
  geom_point(aes(y = e_ols), col = "blue") +
  geom_line(aes(y = e_ols), col = "blue") +
  scale_color_manual(
    values = c("Residual" = "blue")
  ) +
  labs(
    x = "Year",  # Label for the x-axis
    y = "Residuals"  # Label for the y-axis
  ) +
  theme_classic()

# plot some white noise:
set.seed(876573)
white_noise = rnorm(n=nrow(Dtrain), mean = 0, sd = sqrt(sigma2_ols))
qqnorm(white_noise)
qqline(white_noise)

ggplot(train, aes(x=year)) +
  geom_point(aes(y=white_noise), col="blue") +
  geom_line(aes(y=white_noise), col="blue") +
  scale_color_manual(
    values = c("Residual" = "blue")
  ) +
  labs(
    x = "Year",  # Label for the x-axis
    y = "Residuals: White noise"  # Label for the y-axis
  ) +
  theme_classic()

# ---- Question 2.4: WLS model: Weighted Least Squares with matrix notation ----
lambda <- 0.9
weights <- lambda^(nrow(train):1 - 1) # weights become smaller as you go back in time
W <- diag(weights)

# estimate parameters with WLS
theta_WLS <- solve( t(X) %*% W %*% X) %*% t(X) %*% W %*% y
print(theta_WLS)

yhat_wls <- X %*% theta_WLS
e_wls <- y - yhat_wls
RSS_wls <- t(e_wls) %*% W %*% e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - p))

V_wls <- sigma2_wls *  solve( t(X) %*% W %*% X)
print(sqrt(diag(V_wls))) 

# predictions with WLS estimates:
y_pred_wls <- Xtest %*% theta_WLS
Vmatrix_pred <- sigma2_wls * (1 + (Xtest %*% solve( t(X) %*% W %*% X)) %*% t(Xtest) )
y_pred_lwr_wls <- y_pred_wls - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
y_pred_upr_wls <- y_pred_wls + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))

ggplot(train, aes(x=year, y=total)) +
  geom_point(col="black") + 
  geom_line(aes(y=yhat_ols, color="OLS"), size=.5) +
  geom_point(data=test, aes(x=year, y=y_pred, color="OLS"), size=.5) +
  geom_ribbon(data=test, aes(x=year, ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.2, fill="blue") +
  geom_point(data=test, aes(x=year, y=total), col="orange") +
  geom_line(aes(y=yhat_wls, color="WLS"), size=.5) +
  geom_point(data=test, aes(x=year, y=y_pred_wls, color="WLS"), size=.5) +
  geom_ribbon(data=test, aes(x=year, ymin=y_pred_lwr_wls, ymax=y_pred_upr_wls), inherit.aes=FALSE, alpha=0.2, fill="red") +
  labs(title="Model Comparison: OLS vs WLS",
       x="Year",
       y="Total",
       color="Model") +
  scale_color_manual(values=c("OLS"="blue", "WLS"="red")) +
  theme_classic()

# # Was our guess of rho a good guess?  
# # lets calculate the correlation between lag-1 residuals:
# cor(e_wls[1:(n-1)], e_wls[2:n], method="pearson")
# # maybe we could try with rho = 0.7 

# ---- Week03 - From Global to Local model (OLS with fewer datapoints) ----
# 
# # Fit data with less datapoints:
# fit_austa <- function(start, stop) {
#   n <- stop - start + 1
#   p <- 2
#   theta_OLS <- solve(t(X[start:stop,])%*%X[start:stop,])%*%t(X[start:stop,])%*%y[start:stop]
#   
#   yhat <- X%*%theta_OLS
#   Xpred <- cbind(1, seq(1980+stop, 1980+stop+9, 1))
#   ypred <- Xpred%*%theta_OLS
#   pred <- data.frame(Xpred, ypred)
#   
#   RSS_ols <- sum((yhat[start:stop] - y[start:stop])^2)
#   sigma2_ols <- as.numeric(RSS_ols/(n - p))
#   Vmatrix_pred <- sigma2_ols*(1+(Xpred%*%solve(t(X[start:stop,])%*%X[start:stop,]))%*%t(Xpred))
#   pred$ypred_lwr <- pred$ypred - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
#   pred$ypred_upr <- pred$ypred + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
#   
#   myplot <- ggplot(data, aes(x=year, y=int_visit)) +
#     geom_point() + 
#     geom_point(data=data[start:stop,], col="red") + 
#     geom_line(data=data[1:(stop),], aes(y=yhat), col="red", size=.5) + 
#     geom_point(data=pred, aes(x=Xpred[,2], y=ypred), col="red", size=.5) + 
#     geom_ribbon(data=pred, aes(x=Xpred[,2],ymin=ypred_lwr, ymax=ypred_upr), inherit.aes=FALSE, alpha=0.2, fill="red") +
#     xlim(1980, 2020) + ylim(0, 8)
#   
#   print(myplot)
# }
# 
# fit_austa(1,26)
# fit_austa(17,26)
# fit_austa(21,26)
# fit_austa(25,26)
# 
# ---- Question 2.4 (extra): WLS model: WLS with "lambda weights" (exponential weights) ----
lambda = 0.9
weights <- lambda^((n-1):0)
# plot the weights:
# barplot(weights, names=1:26)


SIGMA <- diag(n)
diag(SIGMA) <- 1/weights
# print lower right corner to check:
print(SIGMA[20:26,20:26]) # looks fine :-)


# estimate parameters with WLS
theta_WLS <- solve(t(X)%*%solve(SIGMA)%*%X)%*%(t(X)%*%solve(SIGMA)%*%y)
print(theta_WLS)
yhat_wls <- X%*%theta_WLS

ggplot(train, aes(x=year, y=total)) +
  geom_point(col="black") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5, linetype=2) +
  geom_line(aes(y=yhat_wls), col="blue", size=.5)

# (back to slides)

# predictions from the WLS model:
y_pred_wls <- Xtest%*%theta_WLS

# prediction intervals:
e_wls <- y - yhat_wls
RSS_wls <- t(e_wls)%*%solve(SIGMA)%*%e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - p))
Vmatrix_pred <- sigma2_wls * (1 + (Xtest %*% solve(t(X)%*%solve(SIGMA)%*%X)) %*% t(Xtest) )
y_pred_lwr_wls <- y_pred_wls - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
y_pred_upr_wls <- y_pred_wls + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))

# plot WLS (blue) together with OLS (red) and true test data (black):
ggplot(train, aes(x=year, y=total)) +
  geom_point(col="red") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5, linetype=2) +
  geom_point(data=test, aes(x=year,y=y_pred), col="red", size=.5) +
  geom_ribbon(data=test, aes(x=year,ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.1, fill="red") +
  geom_point(data=test, aes(x=year,y=total), col="black") +
  geom_line(aes(y=yhat_wls), col="blue", size=.5) +
  geom_point(data=test, aes(x=year,y=y_pred_wls), col="blue", size=.5) +
  geom_ribbon(data=test, aes(x=year,ymin=y_pred_lwr_wls, ymax=y_pred_upr_wls), inherit.aes=FALSE, alpha=0.2, fill="blue")

# these prediction intervals are very (too) narrow!

# better prediction intervals:
Tmemory <- sum(1/diag(SIGMA))
print(Tmemory)

sigma2_wls <- as.numeric(RSS_wls/(Tmemory - p))
Vmatrix_pred <- sigma2_wls * (1 + (Xtest %*% solve(t(X)%*%solve(SIGMA)%*%X)) %*% t(Xtest) )
y_pred_lwr_wls <- y_pred_wls - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
y_pred_upr_wls <- y_pred_wls + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))

# plot WLS (blue) together with OLS (red) and true test data (black):
ggplot(train, aes(x=year, y=total)) +
  geom_point(col="red") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5, linetype=2) +
  geom_point(data=test, aes(x=year,y=y_pred), col="red", size=.5) +
  geom_ribbon(data=test, aes(x=year,ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.1, fill="red") +
  geom_point(data=test, aes(x=year,y=total), col="black") +
  geom_line(aes(y=yhat_wls), col="blue", size=.5) +
  geom_point(data=test, aes(x=year,y=y_pred_wls), col="blue", size=.5) +
  geom_ribbon(data=test, aes(x=year,ymin=y_pred_lwr_wls, ymax=y_pred_upr_wls), inherit.aes=FALSE, alpha=0.2, fill="blue")

# ---- Question 4.2: RLS without forgetting: lambda = 1 ----
# Now iterate through data:
# For each step:
# - calculate R and theta
# - calculate one-step prediction

lambda <- 1

# we will use the entire dataset (not only the training data from before):
X <- cbind(1, D$year_norm)
y <- cbind(D$total)

n <- length(X[,1])

# initialise containers for parameter estimates (Theta) and one step predictions:
Theta <- matrix(NA, nrow=n, ncol=p)
OneStepPred <- matrix(NA, nrow=n)

# 1 # very first step:
x1 <- X[1,]

# R_0 = matrix(c(0.1, 0, 0, 0.1), nrow=2, ncol=2, byrow=TRUE)
# R_1 <- lambda*R_0 + x1%*%t(x1) # R is a pxp matrix
R_1 <- x1%*%t(x1) # R is a pxp matrix
h_1 <- x1*y[1]    # h is a px1 vector (but R prints it in a row..)

# to estimate theta we need to invert R:
# Theta_1 = solve(R_1) %*% h_1
# in this very first step R cannot be inverted - too soon to estimate parameters!
# (we cannot estimate p parameters drom only one datapoint)

# 2 # second step - first time to estimate parameters and make prediction
x2 <- X[2,]
R_2 <- lambda*R_1 + x2 %*% t(x2)
h_2 <- lambda*h_1 + x2 * y[2]

solve(R_2)
# R is now invertible (we can estimate p parameters from p observations)

# we estimate theta (for the first time - so not yet using "update" formula):
Theta[2,] <- solve(R_2) %*% h_2

# we predict one step ahead:
OneStepPred[2+1] <- X[2+1,]%*%Theta[2,]

# 3 # third step - first time to use update formula
x3 <- X[3,]
R_3 <- lambda*R_2 + x3 %*% t(x3)
Theta[3,] <- Theta[2,] + solve(R_3) %*% x3 %*% (y[3] - t(x3) %*% Theta[2,])

# we predict one step ahead:
OneStepPred[3+1] <- X[3+1,]%*%Theta[3,]

# next many steps # - update and predict
R <- R_3

for(i in 4:n){
  x <- X[i, ]
  # Update
  R <- lambda*R + x %*% t(x)
  Theta[i, ] <- Theta[i-1, ] + solve(R) %*% x %*% (y[i] - t(x) %*% Theta[i-1, ])
}

# predict
residuals <- numeric(n)
for(i in 4:n-1){
  OneStepPred[i+1] <- X[i+1, ] %*% Theta[i, ]
  residuals[i+1] <- OneStepPred[i+1] - y[i]
}

par(mfrow = c(1, 2))
# Plot estimate of intercept:
plot(Theta[,1])

# Plot estimate of slope:
plot(Theta[,2])

# Plot one step predictions:
ggplot(D, aes(x=year, y=total)) +
  geom_point(col="black", aes(color="Total")) +
  geom_point(aes(y=OneStepPred, color="One Step Prediction"), size=1) + 
  geom_line(aes(y=OneStepPred, color="One Step Prediction")) +
  labs(
    title = "Yearly Data with Predictions, lambda=1",
    x = "Year",
    y = "Total",
    color = "Legend"
  ) +
  scale_color_manual(values = c("Total" = "black", "One Step Prediction" = "blue"))+
  theme_classic()

# Plot one-step ahead residuals
burn_in <- 4
residuals <- residuals[(burn_in + 1):n]

# Time axis for plotting after burn-in
time_axis <- D$year[(burn_in + 1):n]


# Create a data frame for the plot
data <- data.frame(
  time_axis = time_axis,
  residual = residuals)

ggplot(data, aes(x = time_axis)) +
  geom_line(aes(y = residual, color = paste("Lambda =", lambda))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "One-Step Ahead Residuals (RLS without Forgetting)",
    x = "Time",
    y = "Residuals",
    color = "Legend"
  ) + theme_classic()

# ---- Question 4.2: RLS with forgetting: lambda = 0.7 ----
# Now iterate through data:
# For each step:
# - calculate R and theta
# - calculate one-step prediction

lambda <- 0.7

# we will use the entire dataset (not only the training data from before):
X <- cbind(1, D$year_norm)
y <- cbind(D$total)

n <- length(X[,1])

# initialise containers for parameter estimates (Theta) and one step predictions:
Theta <- matrix(NA, nrow=n, ncol=p)
OneStepPred <- matrix(NA, nrow=n)

# 1 # very first step:
x1 <- X[1,]

# R_0 = matrix(c(0.1, 0, 0, 0.1), nrow=2, ncol=2, byrow=TRUE)
# R_1 <- lambda*R_0 + x1%*%t(x1) # R is a pxp matrix
R_1 <- x1%*%t(x1) # R is a pxp matrix
h_1 <- x1*y[1]    # h is a px1 vector (but R prints it in a row..)

# to estimate theta we need to invert R:
# Theta_1 = solve(R_1) %*% h_1
# in this very first step R cannot be inverted - too soon to estimate parameters!
# (we cannot estimate p parameters drom only one datapoint)

# 2 # second step - first time to estimate parameters and make prediction
x2 <- X[2,]
R_2 <- lambda*R_1 + x2 %*% t(x2)
h_2 <- lambda*h_1 + x2 * y[2]

solve(R_2)
# R is now invertible (we can estimate p parameters from p observations)

# we estimate theta (for the first time - so not yet using "update" formula):
Theta[2,] <- solve(R_2) %*% h_2

# we predict one step ahead:
OneStepPred[2+1] <- X[2+1,]%*%Theta[2,]

# 3 # third step - first time to use update formula
x3 <- X[3,]
R_3 <- lambda*R_2 + x3 %*% t(x3)
Theta[3,] <- Theta[2,] + solve(R_3) %*% x3 %*% (y[3] - t(x3) %*% Theta[2,])

# we predict one step ahead:
OneStepPred[3+1] <- X[3+1,]%*%Theta[3,]

# next many steps # - update and predict
R <- R_3

for(i in 4:n){
  x <- X[i, ]
  # Update
  R <- lambda*R + x %*% t(x)
  Theta[i, ] <- Theta[i-1, ] + solve(R) %*% x %*% (y[i] - t(x) %*% Theta[i-1, ])
}

# predict
residuals <- numeric(n)
for(i in 4:n-1){
  OneStepPred[i+1] <- X[i+1, ] %*% Theta[i, ]
  residuals[i+1] <- OneStepPred[i+1] - y[i]
}

par(mfrow = c(1, 2))
# Plot estimate of intercept:
plot(Theta[,1])

# Plot estimate of slope:
plot(Theta[,2])

# Plot one step predictions:
ggplot(D, aes(x=year, y=total)) +
  geom_point(col="black", aes(color="Total")) +
  geom_point(aes(y=OneStepPred, color="One Step Prediction"), size=1) + 
  geom_line(aes(y=OneStepPred, color="One Step Prediction")) +
  labs(
    title = "Yearly Data with Predictions, lambda=0.7",
    x = "Year",
    y = "Total",
    color = "Legend"
  ) +
  scale_color_manual(values = c("Total" = "black", "One Step Prediction" = "blue"))+
  theme_classic()

# Plot one-step ahead residuals
burn_in <- 4
residuals <- residuals[(burn_in + 1):n]

# Time axis for plotting after burn-in
time_axis <- D$year[(burn_in + 1):n]

# Create a data frame for the plot
data <- data.frame(
  time_axis = time_axis,
  residual = residuals)

ggplot(data, aes(x = time_axis)) +
  geom_line(aes(y = residual, color = paste("Lambda =", lambda))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "One-Step Ahead Residuals (RLS with Forgetting)",
    x = "Time",
    y = "Residuals",
    color = "Legend"
  ) + theme_classic()

# ---- Question 4.6: RLS with forgetting: optimization ----

# Function to implement Recursive Least Squares (RLS) with forgetting
rls_forecast <- function(train, x, lambda, k) {
  N <- nrow(train)
  X <- cbind(1, train$year_norm)
  y <- cbind(train$total)

  # 1 # very first step:
  x1 <- X[1,]
  
  # R_0 = matrix(c(0.1, 0, 0, 0.1), nrow=2, ncol=2, byrow=TRUE)
  # R_1 <- lambda*R_0 + x1%*%t(x1) # R is a pxp matrix
  R_1 <- x1%*%t(x1) # R is a pxp matrix
  h_1 <- x1*y[1]    # h is a px1 vector (but R prints it in a row..)
  
  # to estimate theta we need to invert R:
  # Theta_1 = solve(R_1) %*% h_1
  # in this very first step R cannot be inverted - too soon to estimate parameters!
  # (we cannot estimate p parameters drom only one datapoint)
  
  # 2 # second step - first time to estimate parameters and make prediction
  x2 <- X[2,]
  R_2 <- lambda*R_1 + x2 %*% t(x2)
  h_2 <- lambda*h_1 + x2 * y[2]
  
  # solve(R_2)
  # R is now invertible (we can estimate p parameters from p observations)
  
  # we estimate theta (for the first time - so not yet using "update" formula):
  Theta_2 <- solve(R_2) %*% h_2
  
  # Initialise R and theta
  # R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2) # Initial R matrix
  # theta <- matrix(c(0, 0), nrow = 2) # Initial theta vector
  R = R_2
  theta = Theta_2
  h = h_2
  
  # Matrices to store parameter estimates
  theta_history <- matrix(NA, nrow = N, ncol = 2)
  
  # Vector to store one-step-ahead predictions
  y_hat <- numeric(N)
  
  # RLS loop
  for (t in 3:N) {
    # Regressor vector
    # X <- matrix(c(1, x[t-1]), nrow = 2)
    x_ <- X[t,]
    y_ <- y[t]
    
    # Prediction
    y_hat[t] <- t(x_) %*% theta
    
    # Prediction error
    # epsilon <- training_data[t] - y_hat[t]
    
    # Update R
    R <- lambda * R + x_ %*% t(x_)
    
    # Update h
    h <- lambda*h + x_ * y_
    
    # Update theta
    # theta <- theta + solve(R) %*% X %*% epsilon
    theta <- solve(R) %*% h
    
    # Store theta
    theta_history[t,] <- theta
  }
  
  # Make k-step-ahead predictions
  predictions <- numeric(N - k)
  for (a in 1:(N - k)) {
    # X_future <- matrix(c(1, x[a + k]), nrow = 2)  # Regressor vector for k steps ahead
    X_future <- X[a+k,]
    predictions[a] <- t(X_future) %*% theta_history[a, ]  # Use theta estimated at time t
  }
  
  return(list(predictions = predictions, theta_history = theta_history))
}

# Function to calculate RMSE for a given horizon k and lambda
calculate_rmse <- function(training_data, x, lambda, k) {
  # Run RLS and get predictions
  rls_results <- rls_forecast(training_data, x, lambda, k)
  predictions <- rls_results$predictions
  
  # Actual values for the prediction period
  actual_values <- training_data[(k + 1):nrow(training_data),]$total
  
  # Calculate squared errors
  # squared_errors <- (predictions - actual_values)^2
  errors = ifelse(is.na(predictions) | is.na(actual_values), NA, predictions - actual_values)
  squared_errors = errors^2
  
  # Calculate RMSE
  rmse <- sqrt(mean(squared_errors, na.rm = TRUE))
  
  return(rmse)
}

# Set range of lambda values
lambda_values <- seq(0.1, 0.99, by = 0.01)

# Set horizons
horizons <- 1:12

# Matrix to store RMSE values
rmse_matrix <- matrix(NA, nrow = length(lambda_values), ncol = length(horizons),
                      dimnames = list(lambda = lambda_values, horizon = horizons))

# Loop through lambda values and horizons
for (lambda in lambda_values) {
  for (k in horizons) {
    # Calculate RMSE
    rmse <- calculate_rmse(train, X, lambda, k)
    
    # Store RMSE in matrix
    rmse_matrix[as.character(lambda), as.character(k)] <- rmse
  }
}
par(mfrow = c(1, 1))

# Plotting RMSE vs. lambda for each horizon
colors <- rainbow(length(horizons)) #Assign colours to each horizon

plot(lambda_values, rmse_matrix[, 1], type = "l", col = colors[1],
     xlab = "Forgetting Factor (Lambda)", ylab = "RMSE",
     main = "RMSE vs. Lambda for Different Horizons",
     ylim = c(min(rmse_matrix, na.rm = TRUE), max(rmse_matrix, na.rm = TRUE)))

for (i in 2:length(horizons)) {
  lines(lambda_values, rmse_matrix[, i], col = colors[i])
}

legend("topright", legend = paste("Horizon", horizons), col = colors, lty = 1, cex = 0.7)


