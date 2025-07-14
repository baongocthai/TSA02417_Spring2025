## --------------------------------------------------------------------------
## Parameter vector (length 8)
##  par = c( a,           # state-transition scalar
##           b_Ta,        # coeff. on outdoor air temperature
##           b_S,         # coeff. on solar radiation
##           b_I,         # coeff. on electrical load
##           C,           # OBSERVATION gain  (NEW!)
##           log_sigma1,  # log process-noise SD
##           log_sigma2,  # log observation-noise SD
##           x0 )         # initial latent state
## --------------------------------------------------------------------------

kf_logLik_dt_1D <- function(par, df) {
  
  ## --- unpack -------------------------------------------------------------
  a      <- par[1]
  B      <- matrix(par[2:4], nrow = 1)          # 1×3
  C      <- par[5]                              # scalar
  sigma1 <- exp(par[6])
  sigma2 <- exp(par[7])
  x_est  <- par[8]
  P_est  <- 10                                  # vague prior variance
  
  ## --- data ---------------------------------------------------------------
  Y <- df$Y
  U <- as.matrix(df[, c("Ta", "S", "I")])
  Tn <- length(Y)
  
  logLik <- 0
  
  for (t in seq_len(Tn)) {
    
    ##  prediction ----------------------------------------------------------
    x_pred <- a * x_est + as.numeric(B %*% U[t, ])
    P_pred <- a^2 * P_est + sigma1^2            # scalar
    
    ##  innovation ----------------------------------------------------------
    y_pred <- C * x_pred
    S_t    <- C^2 * P_pred + sigma2^2           # scalar
    innov  <- Y[t] - y_pred
    
    logLik <- logLik - 0.5 * (log(2*pi*S_t) + innov^2 / S_t)
    
    ##  update --------------------------------------------------------------
    K_t   <- (P_pred * C) / S_t
    x_est <- x_pred + K_t * innov
    P_est <- (1 - K_t * C) * P_pred
  }
  
  return(as.numeric(logLik))
}

## ---------- convenience wrapper for optim() -------------------------------
estimate_dt_1D <- function(start_par, df, lower = NULL, upper = NULL) {
  optim(
    par     = start_par,
    fn      = function(p) -kf_logLik_dt_1D(p, df),
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(maxit = 1000, trace = 1)
  )
}
