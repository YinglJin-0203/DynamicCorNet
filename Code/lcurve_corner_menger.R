#' Title: L-curve corner search by Menger curvature of every three adjacent point
#' 
#' @param sweep_df: dataframe with stress(loss), movement(penalization) and lambda grid
#' @param log_scale: log-transformation of stress and movement
#' @param plot: produce plot to illustrate the search result
#' 
#' 
#' @returns a list of three elements
#' \item{lambda_star}{a scalar. the value of regularization parameter corresponding to L-curve corner}
#' \item{idx}{a integer. the position of lambda_star on the search grid}
#' \item{kappa_M}{a vector. the Menger curvature at each point  of the L-curve}
#' 
lcurve_corner_menger <- function(sweep_df, log_scale = TRUE, plot = TRUE) {
  
  # --- Curve coordinates ---------------------------------------------------
  if (log_scale) {
    x <- log(sweep_df$stress)
    y <- log(sweep_df$movement)
  } else {
    x <- sweep_df$stress
    y <- sweep_df$movement
  }
  
  n   <- length(x)
  lam <- sweep_df$lambda
  
  # --- Local Menger curvature for each interior point ----------------------
  kappa_M <- rep(NA_real_, n)
  h_vec   <- rep(NA_real_, n)
  R_vec   <- rep(NA_real_, n)
  
  for (k in 2:(n - 1)) {
    P1 <- c(x[k-1], y[k-1])
    P2 <- c(x[k],   y[k])
    P3 <- c(x[k+1], y[k+1])
    
    # Edge lengths
    d12 <- sqrt(sum((P2 - P1)^2))    # |P₁P₂|
    d23 <- sqrt(sum((P3 - P2)^2))    # |P₂P₃|
    d13 <- sqrt(sum((P3 - P1)^2))    # |P₁P₃| (local chord)
    
    # Guard against degenerate triangles (collinear or coincident points)
    if (d12 < 1e-14 || d23 < 1e-14 || d13 < 1e-14) {
      kappa_M[k] <- 0
      h_vec[k]   <- 0
      R_vec[k]   <- Inf
      next
    }
    
    # Triangle area via cross product (signed; take absolute value)
    # A = (1/2)|( P₂−P₁ ) × ( P₃−P₁ )|
    v1  <- P2 - P1
    v2  <- P3 - P1
    A   <- abs(v1[1]*v2[2] - v1[2]*v2[1]) / 2
    
    # Perpendicular distance from P₂ to local chord P₁P₃
    # h = 2A / |P₁P₃|
    h   <- 2 * A / d13
    
    # Circumradius R = (|P₁P₂|·|P₂P₃|·|P₁P₃|) / (4A)
    R   <- (d12 * d23 * d13) / (4 * A)
    
    # Menger curvature κ_M = 1/R = 4A / (|P₁P₂|·|P₂P₃|·|P₁P₃|)
    #                             = 2h / (|P₁P₂|·|P₂P₃|)
    kappa_M[k] <- 1 / R
    h_vec[k]   <- h
    R_vec[k]   <- R
  }
  
  # --- Select corner -------------------------------------------------------
  idx         <- which.max(kappa_M)   # NA values are ignored by which.max
  lambda_star <- lam[idx]
  
  # --- Diagnostic plot -----------------------------------------------------
  if (plot) {
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    
    # Panel 1: L-curve with local triangle at selected corner
    plot(x, y,
         type = "b", pch = 19, cex = 0.7, col = "#444444",
         xlab = if (log_scale) "log(Stress norm)" else "Stress norm",
         ylab = if (log_scale) "log(Movement norm)" else "Movement norm",
         main = "L-curve: Local Menger Corner")
    
    # Mark selected corner
    points(x[idx], y[idx],
           pch = 8, cex = 2, lwd = 2, col = "#e41a1c")
    
    legend("topright", bty = "n",
           pch = 8, col = "#e41a1c", pt.lwd = 2, pt.cex = 1.5,
           legend = sprintf("Corner  λ* = %.4f", lambda_star))
    
    # Panel 2: local Menger curvature profile
    plot(lam, kappa_M,
         type = "b", pch = 19, cex = 0.7, col = "#377eb8",
         xlab = "lambda",
         ylab = "Local Menger curvature  κ_M",
         main = "Local Menger Curvature vs λ",
         na.action = na.omit)
    abline(v   = lambda_star,
           col = "#e41a1c", lty = 2, lwd = 2)
    points(lambda_star, kappa_M[idx],
           pch = 8, cex = 2, lwd = 2, col = "#e41a1c")
    legend("topright", bty = "n",
           legend = sprintf("λ* = %.4f", lambda_star),
           col = "#e41a1c", lty = 2, lwd = 2)
    
    par(op)
  }
  
  list(
    lambda_star = lambda_star,
    idx         = idx,
    kappa_M     = kappa_M,
    # h           = h_vec,
    # R           = R_vec
  )
}