# Load necessary libraries
rm(list = ls())
library(stats)
library(ggplot2)
library(gridExtra)

# Set seed for reproducibility
set.seed(42)

# parameters
n = 200
phi1 = -0.75  #-0.7 #0.6  #-0.6 #-0.6
phi2 = -0.3   #-0.3 #-0.3 #0.3  #0.5
sd = 1

# Function to simulate AR(2) process
AR2.sim <- function(n, phi1, phi2, sd)
{
  eps <- rnorm(n, mean = 0, sd = sd)
  Xt = c(0,0)
  for (t in 3:n)
  {
    Xt[t] = -phi1*Xt[t-1] - phi2*Xt[t-2] + eps[t]
  }
  return(Xt)
}

# Simulate 5 realizations
realizations <- replicate(5, AR2.sim(n, phi1, phi2, sd))

# Plot the realizations
df <- data.frame(Time = 1:n, Realization1 = realizations[,1], Realization2 = realizations[,2], 
                 Realization3 = realizations[,3], Realization4 = realizations[,4], Realization5 = realizations[,5])

# Calculate and plot empirical ACFs
acf_values <- lapply(1:5, function(i) acf(realizations[,i], plot = FALSE, lag.max = 30)$acf)

df_acf <- data.frame(Lag = 0:30, Realization1 = acf_values[[1]], Realization2 = acf_values[[2]], 
                     Realization3 = acf_values[[3]], Realization4 = acf_values[[4]], Realization5 = acf_values[[5]])


# Theoretical ACF using recursion
theoretical_acf <- numeric(31)
theoretical_acf[1] <- 1
theoretical_acf[2] <- -phi1 / (1 + phi2)

for (k in 3:31) {
  theoretical_acf[k] <- -phi1 * theoretical_acf[k - 1] - phi2 * theoretical_acf[k - 2]
}

df_acf$theoretical = theoretical_acf

# Plot 5 realizations

# First plot: 5 Realizations of AR(2) Process
# Define a common color palette
colors <- c("Realization 1" = "blue", 
            "Realization 2" = "red", 
            "Realization 3" = "green", 
            "Realization 4" = "orange", 
            "Realization 5" = "purple", 
            "Theoretical ACF" = "black")
p1 <- ggplot(df, aes(x = Time)) +
  geom_line(aes(y = Realization1, color = "Realization 1")) +
  geom_line(aes(y = Realization2, color = "Realization 2")) +
  geom_line(aes(y = Realization3, color = "Realization 3")) +
  geom_line(aes(y = Realization4, color = "Realization 4")) +
  geom_line(aes(y = Realization5, color = "Realization 5")) +
  labs(title = "5 Realizations of AR(2) Process", x = "Time", y = "Xt") +
  scale_color_manual(values = colors) +
  theme_light()

# Second plot: Empirical & Theoretical ACF
p2 <- ggplot(df_acf, aes(x = Lag)) +
  geom_line(aes(y = Realization1, color = "Realization 1")) +
  geom_line(aes(y = Realization2, color = "Realization 2")) +
  geom_line(aes(y = Realization3, color = "Realization 3")) +
  geom_line(aes(y = Realization4, color = "Realization 4")) +
  geom_line(aes(y = Realization5, color = "Realization 5")) +
  geom_line(aes(y = theoretical, color = "Theoretical ACF"), linetype = "dashed", size = 0.8) +
  labs(title = "Empirical & Theoretical ACF of 5 Realizations of AR(2) Process", x = "Lag", y = "ACF") +
  scale_color_manual(values = colors) +
  theme_light()
# Combine the two plots
grid.arrange(p1, p2, ncol = 1)


