# =============================================================================
# Helper functions for simulating mutlivariate longitudinal data 
# with time-dependent correlation
# Scenarios: (1) Slow/Smooth Dynamics  (2) High Noise / Low Signal
# =============================================================================

# library(MASS)      # mvrnorm
library(Matrix)    # nearPD


#### Data Generation ####

#' Ensure a matrix is a valid correlation matrix
#' M: a square matrix
make_corr <- function(M) {
  M <- (M + t(M)) / 2
  diag(M) <- 1
  pd <- nearPD(M, corr = TRUE, keepDiag = TRUE)
  as.matrix(pd$mat)
}

#' Interpolate between two correlation matrices via convex combination
#' C1, C2: correlation matrices
#' alpha: weight
interp_corr <- function(C1, C2, alpha) make_corr((1 - alpha) * C1 + alpha * C2)

#' Generate a random baseline correlation matrix for p variables
#' p": number of variables
#' spread: standard deviation of generation
random_corr <- function(p, spread = 1) {
  L <- matrix(rnorm(p * p, sd = spread), p, p)
  M <- L %*% t(L)
  cov2cor(M)
}

#' Add symmetric noise to a correlation matrix
#' C: correlation matrix
#' sigma: standard deviation of noise
add_corr_noise <- function(C, sigma) {
  p <- nrow(C)
  E <- matrix(rnorm(p * p, sd = sigma), p, p)
  E <- (E + t(E)) / 2
  make_corr(C + E)
}

#' Generate multivariate observations from a correlation matrix
#'   n_obs : observations per time point
#'   C: correlation matrix
obs_from_corr <- function(C, n_obs) {
  p <- nrow(C)
  X <- MASS::mvrnorm(n_obs, mu = rep(0, p), Sigma = C)
  X   # return sample 
}


# SCENARIO 1: SLOW / SMOOTH DYNAMICS
#   True correlations evolve gradually via sinusoidal interpolation
#   between a small number of anchor matrices.
#   Observation noise is low, so the sample correlation tracks truth well.

gen_smooth <- function(p = 8, T = 60, n_obs = 200, noise_sd = 0.05,
                       n_anchors = 3, seed = 1) {
  set.seed(seed)

  # Build anchor correlation matrices
  anchors <- lapply(seq_len(n_anchors), function(i) random_corr(p))

  # True correlation path: smooth sinusoidal blending between anchors
  true_corrs <- vector("list", T)
  for (t in seq_len(T)) {
    phase  <- (t - 1) / (T - 1) * (n_anchors - 1)   # 0 → n_anchors-1
    seg    <- min(floor(phase), n_anchors - 2)
    alpha  <- phase - seg
    # Smooth (sine-eased) interpolation
    alpha_s <- (1 - cos(pi * alpha)) / 2
    true_corrs[[t]] <- interp_corr(anchors[[seg + 1]], anchors[[seg + 2]], alpha_s)
  }

  # Observed correlations: true + small noise
  obs_corrs <- lapply(true_corrs, function(C) {
    Cn <- add_corr_noise(C, noise_sd)
    obs_from_corr(Cn, n_obs)
  })

  list(
    scenario   = "smooth",
    p          = p,
    T          = T,
    n_obs      = n_obs,
    noise_sd   = noise_sd,
    true_corrs = true_corrs,
    obs_corrs  = obs_corrs
  )
}

#### Result Analysis ####

#### L-CURVE CORNER DETECTION ####


# Moving average 
lcurve_corner_ma <- function(sweep_df) {
  s <- log(sweep_df$stress)
  m <- log(sweep_df$movement)
  lam <- log(sweep_df$lambda)
  
  # Smooth before differentiating (3-point moving average)
  smooth3 <- function(x) {
    n <- length(x)
    c(x[1], (x[1:(n-2)] + x[2:(n-1)] + x[3:n]) / 3, x[n])
  }
  s <- smooth3(s); m <- smooth3(m)
  
  # Numerical derivatives w.r.t. log(lambda)
  ds <- diff(s) / diff(lam); ds <- c(ds[1], ds)
  dm <- diff(m) / diff(lam); dm <- c(dm[1], dm)
  d2s <- diff(ds) / diff(lam); d2s <- c(d2s[1], d2s)
  d2m <- diff(dm) / diff(lam); d2m <- c(d2m[1], d2m)
  
  # Curvature
  kappa <- abs(d2s * dm - ds * d2m) / (ds^2 + dm^2)^1.5
  kappa[!is.finite(kappa)] <- 0
  
  idx <- which.max(kappa)
  list(
    lambda_star = sweep_df$lambda[idx],
    idx         = idx,
    kappa       = kappa
  )
}


# L-CURVE CORNER DETECTION: max perpendicular distance
lcurve_corner_dist <- function(sweep_df) {
  s <- log(sweep_df$stress)
  m <- log(sweep_df$movement)
  
  # Normalize both axes to [0,1] before computing distance
  s_n <- (s - min(s)) / (max(s) - min(s))
  m_n <- (m - min(m)) / (max(m) - min(m))
  
  # Line from first to last point
  p1 <- c(s_n[1], m_n[1])
  p2 <- c(s_n[length(s_n)], m_n[length(m_n)])
  d  <- p2 - p1
  
  # Perpendicular distance from each point to the line
  dist_to_line <- abs(d[2] * s_n - d[1] * m_n + (p2[1]*p1[2] - p2[2]*p1[1])) /
    sqrt(sum(d^2))
  
  idx <- which.max(dist_to_line)
  list(lambda_star = sweep_df$lambda[idx], idx = idx,
       dist_to_line = dist_to_line)
}



# Menger curvature
menger_corner <- function(sweep_df, log_scale = TRUE, plot = TRUE) {
  
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
    h           = h_vec,
    R           = R_vec
  )
}


# splines corner detection
lcorner_splines <- function(sweep_df, log_scale = TRUE, spar = NULL, plot = TRUE) {
  
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

#### Diagnostics ####

# DIAGNOSTIC PLOTS

plot_sweep <- function(sweep_df, scenario_name, lambda_star = NULL) {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  lam_log <- log10(sweep_df$lambda)
  
  # 1. Stress vs log(lambda)
  plot(lam_log, sweep_df$stress, type = "b", pch = 19, col = "#2166ac",
       xlab = "log10(lambda)", ylab = "Stress", main = "Stress")
  if (!is.null(lambda_star))
    abline(v = log10(lambda_star), col = "red", lty = 2, lwd = 2)
  
  # 2. Movement vs log(lambda)
  plot(lam_log, sweep_df$movement, type = "b", pch = 19, col = "#d6604d",
       xlab = "log10(lambda)", ylab = "Movement", main = "Movement")
  if (!is.null(lambda_star))
    abline(v = log10(lambda_star), col = "red", lty = 2, lwd = 2)
  
  # 3. Objective vs log(lambda)
  plot(lam_log, sweep_df$objective, type = "b", pch = 19, col = "#4dac26",
       xlab = "log10(lambda)", ylab = "Objective", main = "Joint Objective")
  if (!is.null(lambda_star))
    abline(v = log10(lambda_star), col = "red", lty = 2, lwd = 2)
  
  # 4. L-curve: log(Stress) vs log(Movement)
  plot(log(sweep_df$stress), log(sweep_df$movement),
       type = "b", pch = 19, col = "#762a83",
       xlab = "log(Stress)", ylab = "log(Movement)", main = "L-Curve")
  if (!is.null(lambda_star)) {
    idx <- which.min(abs(sweep_df$lambda - lambda_star))
    points(log(sweep_df$stress[idx]), log(sweep_df$movement[idx]),
           pch = 8, cex = 2, col = "red", lwd = 2)
  }
  
  mtext(paste("Scenario:", scenario_name), outer = TRUE, line = -1.5, cex = 1.1)
  par(op)
}

# PLOT STRESS PLATEAU DIAGNOSTIC

plot_stress_plateau <- function(sweep_df, plateau, scenario_name) {
  op  <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  lam_log <- log10(sweep_df$lambda)
  
  # Panel 1: normalized stress with plateau upper bound marked
  plot(lam_log, sweep_df$stress_norm,
       type = "b", pch = 19, col = "#2166ac",
       xlab = "log10(lambda)", ylab = "Normalized Stress",
       main = paste("Stress Plateau —", scenario_name))
  abline(v   = log10(plateau$lambda_ub),
         col = "red", lty = 2, lwd = 2)
  legend("topleft",
         legend = sprintf("lambda_UB = %.4f", plateau$lambda_ub),
         col = "red", lty = 2, lwd = 2, bty = "n")
  
  # Panel 2: marginal stress increase with threshold
  plot(lam_log, plateau$delta_s,
       type = "b", pch = 19, col = "#d6604d",
       xlab = "log10(lambda)",
       ylab = "Marginal stress increase / unit log(lambda)",
       main = "Marginal Stress Rate")
  abline(h   = plateau$threshold,
         col = "darkgreen", lty = 2, lwd = 2)
  abline(v   = log10(plateau$lambda_ub),
         col = "red", lty = 2, lwd = 2)
  legend("topright",
         legend = c(sprintf("threshold (eps=%.2f)", 0.05),
                    sprintf("lambda_UB = %.4f", plateau$lambda_ub)),
         col = c("darkgreen", "red"), lty = 2, lwd = 2, bty = "n")
  
  par(op)
}


