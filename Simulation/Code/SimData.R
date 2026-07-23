# =============================================================================
# Dynamic MDS: Simulated Multivariate Longitudinal Data
# Scenarios: (1) Slow/Smooth Dynamics  (2) High Noise / Low Signal
# =============================================================================


# -----------------------------------------------------------------------------
#### SCENARIO 1: SLOW / SMOOTH DYNAMICS ####
#   True correlations evolve gradually via sinusoidal interpolation
#   between a small number of anchor matrices.
#   Observation noise is low, so the sample correlation tracks truth well.
# -----------------------------------------------------------------------------

gen_smooth <- function(p = 10, T = 10, n_obs = 100, noise_sd = 0.05,
                       n_anchors = 3, seed = 1) {
  set.seed(seed)
  
  # Build anchor correlation matrices
  anchors <- lapply(seq_len(n_anchors), function(i) random_corr(p,))
  
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
  obs <- lapply(true_corrs, function(C) {
    # Cn <- add_corr_noise(C, noise_sd)
    obs_from_corr(C, n_obs)
  })
  
  list(
    scenario   = "smooth",
    p          = p,
    T          = T,
    n_obs      = n_obs,
    noise_sd   = noise_sd,
    true_corrs = true_corrs,
    obs  = obs
  )
}


# --------------true_corrs# -----------------------------------------------------------------------------
# SCENARIO 2: HIGH NOISE / LOW SIGNAL
#   True correlations still evolve smoothly, but observations are noisy:
#   (a) small n_obs  → high sampling variance in sample correlations
#   (b) large noise_sd → direct perturbation of the correlation matrix
#   The signal (true temporal change) is modest relative to noise.
# -----------------------------------------------------------------------------

gen_high_noise <- function(p = 8, T = 60, n_obs = 25, noise_sd = 0.25,
                           n_anchors = 3, signal_strength = 0.3, seed = 2) {
  set.seed(seed)
  
  # Anchors placed close together → weak true signal
  anchors <- lapply(seq_len(n_anchors), function(i) {
    base <- random_corr(p, spread = 0.4)
    interp_corr(base, random_corr(p, spread = 0.4), signal_strength)
  })
  
  # True correlation path (same smooth interpolation as Scenario 1)
  true_corrs <- vector("list", T)
  for (t in seq_len(T)) {
    phase  <- (t - 1) / (T - 1) * (n_anchors - 1)
    seg    <- min(floor(phase), n_anchors - 2)
    alpha  <- phase - seg
    alpha_s <- (1 - cos(pi * alpha)) / 2
    true_corrs[[t]] <- interp_corr(anchors[[seg + 1]], anchors[[seg + 2]], alpha_s)
  }
  
  # Observed: few observations + large additive noise
  obs_corrs <- lapply(true_corrs, function(C) {
    Cn <- add_corr_noise(C, noise_sd)
    obs_from_corr(Cn, n_obs)
  })
  
  list(
    scenario        = "high_noise",
    p               = p,
    T               = T,
    n_obs           = n_obs,
    noise_sd        = noise_sd,
    signal_strength = signal_strength,
    true_corrs      = true_corrs,
    obs_corrs       = obs_corrs
  )
}