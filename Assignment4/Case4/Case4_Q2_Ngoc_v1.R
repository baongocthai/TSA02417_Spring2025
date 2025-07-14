# ---- Question 2.1: Load data ----
# Load necessary libraries
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)

# Load the data
data <- read_csv("assignment4/transformer_data.csv")

# Rename columns for clarity (if needed)
colnames(data) <- c("Time", "Transformer_Temp_Yt", "Outdoor_Temp_Ta", "Solar_Radiation_PhiS", "Load_PhiI")

# Pivot to long format with labels
long_data <- data %>%
  pivot_longer(cols = -Time, names_to = "Variable", values_to = "Value") %>%
  mutate(
    Variable = factor(Variable, levels = c(
      "Transformer_Temp_Yt", "Outdoor_Temp_Ta", "Solar_Radiation_PhiS", "Load_PhiI"
    ),
    labels = c(
      "Transformer Temp Y[t] (°C)",
      "Outdoor Temp T[a,t] (°C)",
      "Solar Radiation Φ[s,t] (W/m²)",
      "Load Φ[I,t] (kW)"
    ))
  )

# Assign a distinct color per variable
variable_colors <- c(
  "Transformer Temp Y[t] (°C)" = "firebrick",
  "Outdoor Temp T[a,t] (°C)" = "dodgerblue",
  "Solar Radiation Φ[s,t] (W/m²)" = "orange",
  "Load Φ[I,t] (kW)" = "forestgreen"
)

# Plot with separate color for each variable
ggplot(long_data, aes(x = Time, y = Value, color = Variable)) +
  geom_line() +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  scale_color_manual(values = variable_colors, guide = "none") +
  labs(
    title = "Transformer Station Time Series Data",
    x = "Time (hours)",
    y = NULL
  ) +
  theme_minimal(base_size = 13)

# ---- Question 2.2: 1D state space model ----
## Load & inspect the data
## data
df <- read_csv("assignment4/transformer_data.csv")

## starting values and bounds (now 8 parameters)
start_par <- c(
  0.80,      # a
  0.10,      # b_Ta
  0.001,     # b_S
  0.10,      # b_I
  1.00,      # C   (start at 1 because Y and X are both in °C)
  log(1.0),  # log_sigma1
  log(1.0),  # log_sigma2
  df$Y[1]    # x0
)

lower <- c(-0.99,  rep(-5, 3),   0.50,  log(1e-5), log(1e-5),  min(df$Y) - 20)
upper <- c( 0.99,  rep( 5, 3),   1.50,  log(1e+2), log(1e+2),  max(df$Y) + 20)

fit <- estimate_dt_1D(start_par, df, lower, upper)
print(fit$par)      # MLEs on the original scale

##  Diagnostics for the 1-D model *with* observation gain  C
library(stats)        # acf(), qqnorm()

par_hat <- fit$par        # <- result from estimate_dt()
logLik  <- kf_logLik_dt_1D(par_hat, df)

## 1-step-ahead innovations and their variances
kal_out_1D <- local({
  a   <- par_hat[1]
  B   <- matrix(par_hat[2:4], 1)   # 1×3
  C   <- par_hat[5]
  s1  <- exp(par_hat[6])
  s2  <- exp(par_hat[7])
  x   <- par_hat[8]
  P   <- 10                        # same vague prior
  inn <- var <- numeric(nrow(df))
  y_preds <- numeric(nrow(df))
  
  for (t in seq_len(nrow(df))) {
    ## ---------- prediction ----------
    x_pred <- a * x + as.numeric(B %*% as.numeric(df[t, c("Ta", "S", "I")]))
    P_pred <- a^2 * P + s1^2
    
    ## ---------- innovation ----------
    y_pred <- C * x_pred
    y_preds[t] <- y_pred
    var[t] <- C^2 * P_pred + s2^2
    inn[t] <- df$Y[t] - y_pred
    
    ## ---------- update --------------
    K       <- (P_pred * C) / var[t]
    x       <- x_pred + K * inn[t]
    P       <- (1 - K * C) * P_pred
  }
  list(innov = inn, innov_var = var, preds = y_preds)
})
residuals = kal_out_1D$innov
y_preds = kal_out_1D$preds
std_resid <- kal_out_1D$innov / sqrt(kal_out_1D$innov_var)

df_pred_1D = as.data.frame(y_preds)
df_pred_1D$time = df$time
df_pred_1D$obs = df$Y
df_pred_1D$residuals = residuals

# par(mfrow = c(2,2))
# plot(std_resid, type = "l", main = "Standardised residuals")
# acf(std_resid, main = "ACF of residuals")
# pacf(std_resid, main = "PACF of residuals")
# qqnorm(std_resid); qqline(std_resid)

k   <- length(par_hat)            # now 8 parameters
n   <- nrow(df)
AIC <- -2 * logLik + 2 * k
BIC <- -2 * logLik + k * log(n)
cat(sprintf("log-Lik = %.2f  |  AIC = %.2f  |  BIC = %.2f\n",
            logLik, AIC, BIC))

# --- Plot observed vs predicted temperature ---
# Load required libraries
library(ggplot2)
library(gridExtra)
library(forecast)  # for ggAcf and ggPacf
library(grid)

# --- Plot 1: Observed vs Predicted ---
p1 <- ggplot(df_pred_1D, aes(x = time)) +
  geom_line(aes(y = obs, color = "Observed"), size = 1) +
  geom_line(aes(y = y_preds, color = "Predicted"), size = 1, linetype = "dashed") +
  labs(title = "Observed vs Predicted Transformer Temperature: 1D state-space model",
       y = "Temperature (°C)", color = "Legend") +
  theme_minimal() +
  theme(
    legend.position = c(0.1, 0.8),
    legend.background = element_rect(fill = alpha('white', 0.8), color = NA),
    legend.box.background = element_rect(color = "grey70")
  )

# --- Plot 2: Residual Time Series ---
p2 <- ggplot(df_pred_1D, aes(x = time, y = residuals)) +
  geom_line(color = "red") +
  labs(title = "Model Residuals Over Time", y = "Residual (°C)", x = "Time") +
  theme_minimal()

# --- Plot 3: QQ Plot with ggplot2 ---
p3 <- ggplot(df_pred_1D, aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line(color = "blue") +
  labs(title = "QQ Plot of Residuals") +
  theme_minimal()

# --- Plot 4: ACF Plot with ggAcf ---
p4 <- ggAcf(df_pred_1D$residuals, main = "ACF of Residuals") +
  theme_minimal()

# --- Plot 5: PACF Plot with ggPacf ---
p5 <- ggPacf(df_pred_1D$residuals, main = "PACF of Residuals") +
  theme_minimal()

# --- Arrange all plots ---
grid.arrange(
  p1,
  arrangeGrob(p2, p3, ncol = 2),
  arrangeGrob(p4, p5, ncol = 2),
  nrow = 3,
  heights = c(2, 1.5, 2)
)

# Display combined figure
# print(final_plot)

# ---- Question 2.3: 2D state space model ----
df <- read_csv("assignment4/transformer_data.csv")

# starting values & bounds 
# start_par <- c(
#   ## A 
#   0.80, 0.05,
#   0.10, 0.70,
#   ## B 
#   0.12, 0.002, 0.30,     # first row (Ta, S, I)
#   0.05, 0.001, 0.10,     # second row
#   ## C 
#   1.0,  0.0,             # mostly first state
#   ## log σ's 
#   log(1), log(1), log(0.05),
#   ## x0 
#   df$Y[1], df$Y[1] - 1,
#   log(10)                # vague prior var
# )
# 
# lower <- c(rep(-0.99, 4),   # A
#            rep(-5, 6),      # B
#            -2, -2,          # C
#            log(1e-5), log(1e-5), log(1e-5),   # σ's
#            rep(-50, 2), log(1e-2))
# 
# upper <- c(rep( 0.99, 4),
#            rep( 5,   6),
#            2,  2,
#            log(1e2), log(1e2), log(1e2),
#            rep(150, 2), log(1e3))

start_par <- c(
  ## A 
  0.80, 0.05,
  0.10, 0.70,
  ## B 
  0.12, 0.002, 0.30,     # first row (Ta, S, I)
  0.05, 0.001, 0.05,     # second row
  ## C 
  1,  0,             # mostly first state
  ## log σ's 
  log(1), log(1), log(0.05),
  ## x0 
  df$Y[1], df$Y[1] - 1,
  log(10)                # vague prior var
)

lower <- c(rep(-1, 4),   # A
           rep(-1, 6),      # B
           -1, -1,          # C
           log(1e-5), log(1e-5), log(1e-5),   # σ's
           rep(min(df$Y) - 20, 2), log(1e-2))

upper <- c(rep( 1, 4),
           rep( 1, 6),
           1,  1,
           log(1e2), log(1e2), log(1e2),
           rep(max(df$Y) + 100, 2), log(1e3))

fit2d <- estimate_dt_2D(start_par, df, lower, upper)

library(stats)

par_hat <- fit2d$par
logLik  <- kf_logLik_dt_2D(par_hat, df)   # maximised log-likelihood

## ----- run filter once to get innovations 
A <- matrix(par_hat[1:4],  2, 2, byrow = TRUE)
B <- matrix(par_hat[5:10], 2, 3, byrow = TRUE)
C <- matrix(par_hat[11:12], 1, 2)

Sigma1 <- diag(exp(par_hat[13:14])^2, 2)
Sigma2 <- exp(par_hat[15])^2

x  <- matrix(par_hat[16:17], 2, 1)
P  <- diag(exp(par_hat[18])^2, 2)

y_preds <- inn <- var <- numeric(nrow(df))

for (t in seq_len(nrow(df))) {
  x_pred <- A %*% x + B %*% t(as.matrix(df[t, c("Ta","S","I")]))
  P_pred <- A %*% P %*% t(A) + Sigma1
  S_t    <- C %*% P_pred %*% t(C) + Sigma2
  inn[t] <- df$Y[t] - C %*% x_pred
  var[t] <- S_t
  K      <- P_pred %*% t(C) / as.numeric(S_t)
  x      <- x_pred + K * inn[t]
  P      <- (diag(2) - K %*% C) %*% P_pred
  y_pred <- C %*% x_pred
  y_preds[t] <- y_pred
}

std_resid <- inn / sqrt(var)

## ----- plots 
par(mfrow = c(2,2))
plot(std_resid, type = "l", main = "Standardised residuals")
acf(std_resid, main = "ACF of residuals")
pacf(std_resid, main = "PACF of residuals")
qqnorm(std_resid); qqline(std_resid)

## ----- information criteria 
k   <- length(par_hat)             # 18 parameters
n   <- nrow(df)
AIC <- -2*logLik + 2*k
BIC <- -2*logLik + k*log(n)
cat(sprintf("log-Lik = %.2f  |  AIC = %.2f  |  BIC = %.2f\n",
            logLik, AIC, BIC))

residuals = inn
y_preds = y_preds

df_pred_2D = as.data.frame(y_preds)
df_pred_2D$time = df$time
df_pred_2D$obs = df$Y
df_pred_2D$residuals = residuals

library(ggplot2)
library(gridExtra)
library(forecast)  # for ggAcf and ggPacf
library(grid)

# --- Plot 1: Observed vs Predicted ---
p1 <- ggplot(df_pred_2D, aes(x = time)) +
  geom_line(aes(y = obs, color = "Observed"), size = 1) +
  geom_line(aes(y = y_preds, color = "Predicted"), size = 1, linetype = "dashed") +
  labs(title = "Observed vs Predicted Transformer Temperature: 2D state-space model",
       y = "Temperature (°C)", color = "Legend") +
  theme_minimal() +
  theme(
    legend.position = c(0.1, 0.8),
    legend.background = element_rect(fill = alpha('white', 0.8), color = NA),
    legend.box.background = element_rect(color = "grey70")
  )

# --- Plot 2: Residual Time Series ---
p2 <- ggplot(df_pred_2D, aes(x = time, y = residuals)) +
  geom_line(color = "red") +
  labs(title = "Model Residuals Over Time", y = "Residual (°C)", x = "Time") +
  theme_minimal()

# --- Plot 3: QQ Plot with ggplot2 ---
p3 <- ggplot(df_pred_2D, aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line(color = "blue") +
  labs(title = "QQ Plot of Residuals") +
  theme_minimal()

# --- Plot 4: ACF Plot with ggAcf ---
p4 <- ggAcf(df_pred_2D$residuals, main = "ACF of Residuals") +
  theme_minimal()

# --- Plot 5: PACF Plot with ggPacf ---
p5 <- ggPacf(df_pred_2D$residuals, main = "PACF of Residuals") +
  theme_minimal()

# --- Arrange all plots ---
grid.arrange(
  p1,
  arrangeGrob(p2, p3, ncol = 2),
  arrangeGrob(p4, p5, ncol = 2),
  nrow = 3,
  heights = c(2, 1.5, 2))


# ---- Question 2.4: Reconstruct ----
for (t in seq_len(nrow(df))) {
  ## prediction
  x_pred <- A %*% x_est + B %*% t(as.matrix(df[t, c("Ta","S","I")]))
  P_pred <- A %*% P_est %*% t(A) + Sigma1
  
  ## update
  S_t   <- C %*% P_pred %*% t(C) + Sigma2
  K_t   <- P_pred %*% t(C) / as.numeric(S_t)
  innov <- df$Y[t] - C %*% x_pred
  x_est <- x_pred + K_t %*% innov
  P_est <- (diag(2) - K_t %*% C) %*% P_pred
  
  x_filt[t, ] <- x_est[,1]
}

## --- gather into long data.frame for ggplot 
plot_df <- data.frame(
  time = df$time,
  Y      = df$Y,
  Ta     = df$Ta,
  Solar  = df$S,
  Load   = df$I,
  State1 = x_filt[,1],
  State2 = x_filt[,2]
)

library(tidyr)
states_long <- pivot_longer(plot_df,
                            cols = c(State1, State2),
                            names_to = "State", values_to = "value")

inputs_long <- pivot_longer(plot_df,
                            cols = c(Ta, Solar, Load),
                            names_to = "Input", values_to = "value")

## --- Figure 1: both states 
ggplot(states_long, aes(time, value, colour = State)) +
  geom_line() +
  labs(y = "Latent state (°C)", x = "Time (h)",
       title = "Reconstructed latent states – 2-D model")

## --- Figure 2: states + inputs 
ggplot() +
  geom_line(data = states_long,
            aes(time, value, colour = State), size = 1) +
  geom_line(data = inputs_long,
            aes(time, scale(value),  linetype = Input), alpha = 1) +
  labs(y = "Scaled value", x = "Time (h)",
       title = "States and (scaled) inputs") +
  scale_linetype_manual(values = c("dotted","dashed","dotdash")) +
  theme(legend.position = "bottom")+
  theme_minimal()


