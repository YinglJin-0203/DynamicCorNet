#' Title: L-curve corner search of analytically curvature from splines-smoothed L-curve
#' 
#' @param sweep_df: dataframe with stress(loss), movement(penalization) and lambda grid
#' @param log_scale: log-transformation of stress and movement
#' @param plot: produce plot to illustrate the search result
#' @param spar: a scalar or a list of two scalars for splines smoothing parameter of stress and movement
#' 
lcurve_corner_splines <- function(sweep_df, log_scale = TRUE, spar = NULL, plot = TRUE) {
  
  t <- sweep_df$lambda
  
  # --- L-curve coordinates -------------------------------------------------
  if (log_scale) {
    s <- log(sweep_df$stress)
    m <- log(sweep_df$movement)
  } else {
    s <- sweep_df$stress
    m <- sweep_df$movement
  }
  
  # --- Parse spar ----------------------------------------------------------
  # Allow a single scalar or a named list(s = ..., m = ...)
  if (is.list(spar)) {
    spar_s <- spar$s
    spar_m <- spar$m
  } else {
    spar_s <- spar
    spar_m <- spar
  }
  
  # --- Fit independent splines ------------------------
  # s(t): stress as a function of lambda
  # m(t): movement as a function lambda
  spline_s <- smooth.spline(t, s, spar = spar_s)
  spline_m <- smooth.spline(t, m, spar = spar_m)
  
  # --- Evaluate splines and their derivatives on a fine grid ---------------
  t_fine <- seq(min(t), max(t), length.out = 500)
  
  s0_fine <- predict(spline_s, t_fine, deriv = 0)$y   # s(t)
  s1_fine <- predict(spline_s, t_fine, deriv = 1)$y   # s'(t)
  s2_fine <- predict(spline_s, t_fine, deriv = 2)$y   # s''(t)
  
  m0_fine <- predict(spline_m, t_fine, deriv = 0)$y   # m(t)
  m1_fine <- predict(spline_m, t_fine, deriv = 1)$y   # m'(t)
  m2_fine <- predict(spline_m, t_fine, deriv = 2)$y   # m''(t)
  
  # --- Parametric curvature on fine grid -----------------------------------
  # κ(t) = |s' m'' − m' s''| / (s'² + m'²)^(3/2)
  num_fine   <- abs(s1_fine * m2_fine - m1_fine * s2_fine)
  denom_fine <- (s1_fine^2 + m1_fine^2)^1.5
  denom_fine[denom_fine < 1e-14] <- 1e-14       # guard against flat regions
  kappa_fine <- num_fine / denom_fine
  
  # --- Evaluate at original grid points for selection ----------------------
  s1_grid <- predict(spline_s, t, deriv = 1)$y
  s2_grid <- predict(spline_s, t, deriv = 2)$y
  m1_grid <- predict(spline_m, t, deriv = 1)$y
  m2_grid <- predict(spline_m, t, deriv = 2)$y
  
  num_grid   <- abs(s1_grid * m2_grid - m1_grid * s2_grid)
  denom_grid <- (s1_grid^2 + m1_grid^2)^1.5
  denom_grid[denom_grid < 1e-14] <- 1e-14
  kappa_grid <- num_grid / denom_grid
  
  # --- Corner: maximum curvature at original grid points -------------------
  idx         <- which.max(kappa_grid)
  lambda_star <- t[idx]
  
  # --- Diagnostic plot -----------------------------------------------------
  if (plot) {
    op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    
    # Panel 1: L-curve with spline fit and selected corner
    plot(s, m,
         type = "p", pch = 19, cex = 0.7, col = "#444444",
         xlab = if (log_scale) "log(Stress)" else "Stress norm",
         ylab = if (log_scale) "log(Movement)" else "Movement norm",
         main = "L-curve")
    lines(s0_fine, m0_fine, col = "#888888", lwd = 2)   # spline path
    points(s[idx], m[idx],
           pch = 8, cex = 2, lwd = 2, col = "#e41a1c")
    legend("topright", bty = "n",
           pch = c(19, 8), col = c("#444444", "#e41a1c"),
           pt.lwd = c(1, 2), pt.cex = c(0.9, 1.5),
           legend = c("Data", sprintf("Corner  λ*=%.4f", lambda_star)))
    
    # Panel 2: curvature profile on fine grid
    plot(t_fine, kappa_fine,
         type = "l", lwd = 2, col = "#377eb8",
         xlab = "lambda",
         ylab = "Curvature  κ(λ)",
         main = "Spline Curvature Profile")
    points(t, kappa_grid, pch = 19, cex = 0.6, col = "#444444")
    abline(v   = lambda_star, col = "#e41a1c", lty = 2, lwd = 2)
    points(t[idx], kappa_grid[idx],
           pch = 8, cex = 2, lwd = 2, col = "#e41a1c")
    legend("topright", bty = "n",
           legend = c("Spline κ(λ)", "Grid points",
                      sprintf("λ* = %.4f", lambda_star)),
           col    = c("#377eb8", "#444444", "#e41a1c"),
           lty    = c(1, NA, 2), pch = c(NA, 19, NA),
           lwd    = 2, pt.cex = 0.8)
    
    # Panel 3: individual spline fits for s(t) and m(t)
    s0_grid <- predict(spline_s, t, deriv = 0)$y
    m0_grid <- predict(spline_m, t, deriv = 0)$y
    
    ylim_range <- range(c(s, m, s0_fine, m0_fine))
    plot(t, s,
         type = "p", pch = 19, cex = 0.7, col = "#2166ac",
         ylim = ylim_range,
         xlab = "lambda",
         ylab = if (log_scale) "log(norm)" else "norm",
         main = "Spline Fits")
    lines(t_fine, s0_fine, col = "#2166ac", lwd = 2)
    points(t, m, pch = 19, cex = 0.7, col = "#d6604d")
    lines(t_fine, m0_fine, col = "#d6604d", lwd = 2)
    abline(v = lambda_star, col = "#e41a1c", lty = 2, lwd = 2)
    legend("right", bty = "n",
           col = c("#2166ac", "#d6604d"),
           lty = 1, lwd = 2,
           legend = c("Stress spline", "Movement spline"))
    
    par(op)
  }
  
  list(
    lambda_star = lambda_star,
    idx         = idx,
    kappa       = kappa_grid,
    kappa_fine  = kappa_fine,
    t_fine      = t_fine,
    spline_s    = spline_s,
    spline_m    = spline_m
  )
}
