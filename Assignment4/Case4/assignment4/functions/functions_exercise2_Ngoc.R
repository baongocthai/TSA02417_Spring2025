## --------------------------------------------------------------------------
## 2-D Kalman filter log-likelihood for transformer data
##
## Parameter vector (length 18)
##  par = c( vec(A   ),          # a11 a12 a21 a22           (4)
##           vec(B   ),          # b11 b12 b13 b21 b22 b23   (6)
##           C1, C2,             # observation row-vector    (2)
##           log_sigma11,        # s.d. of process noise 1
##           log_sigma22,        # s.d. of process noise 2   (diagonal Σ1)
##           log_sigma2,         # s.d. of obs. noise        (scalar Σ2)
##           x01, x02,           # initial state mean
##           log_P0 )            # log prior variance (diag) (1)   → total 18
## --------------------------------------------------------------------------

kf_logLik_dt_2D <- function(par, df) {
  
  ## ----- unpack -----------------------------------------------------------
  A <- matrix(par[1:4],  2, 2, byrow = TRUE)
  B <- matrix(par[5:10], 2, 3, byrow = TRUE)
  C <- matrix(par[11:12], 1, 2)
  
  Sigma1 <- diag(exp(par[13:14])^2, 2)      # 2×2 diagonal
  Sigma2 <- matrix(exp(par[15])^2, 1, 1)    # scalar (1×1)
  
  x_est <- matrix(par[16:17], 2, 1)         # initial mean
  P_est <- diag(exp(par[18])^2, 2)          # initial covariance
  
  ## ----- data -------------------------------------------------------------
  Y  <- as.matrix(df$Y)                       # 1×T
  U  <- as.matrix(df[, c("Ta","S","I")])      # T×3
  Tn <- nrow(df)
  
  logLik <- 0
  I2 <- diag(2)
  
  for (t in seq_len(Tn)) {
    
    ## ---- predict ---------------------------------------------------------
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], 3, 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1
    
    ## ---- innovation ------------------------------------------------------
    y_pred <- C %*% x_pred
    S_t    <- C %*% P_pred %*% t(C) + Sigma2
    innov  <- Y[t] - y_pred
    
    logLik <- logLik - 0.5*(log(2*pi*det(S_t)) + innov^2 / S_t)
    
    ## ---- update ----------------------------------------------------------
    K_t   <- P_pred %*% t(C) / as.numeric(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- (I2 - K_t %*% C) %*% P_pred
  }
  
  as.numeric(logLik)        # maximiser will *maximise* this
}

## ---------- optimiser wrapper --------------------------------------------
estimate_dt_2D <- function(start_par, df, lower = NULL, upper = NULL) {
  optim(
    par     = start_par,
    fn      = function(p) -kf_logLik_dt_2D(p, df),   # minimise –logL
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(maxit = 2000, trace = 1)
  )
}
