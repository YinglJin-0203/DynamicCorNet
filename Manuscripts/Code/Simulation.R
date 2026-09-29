# -----------------------------------------------------------------------------
#### Set up ####
# -----------------------------------------------------------------------------


library(MASS)      # mvrnorm
library(Matrix)    # nearPD

source("Code/dyn_mds.R")
source("Code/lambda_sweep.R")
source("Code/dMDS_Helpers.R")
source("Code/get_similarity.R")
source("Code/lcurve_corner_dist.R")
source("Code/lcurve_corner_menger.R")


# -----------------------------------------------------------------------------
#### Function to simulate data ####
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

# -----------------------------------------------------------------------------
#### Functions to generated correlation matrices ####
# -----------------------------------------------------------------------------

#' Ensure a matrix is a valid correlation matrix (PSD)
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
  X <- mvrnorm(n_obs, mu = rep(0, p), Sigma = C)
  X   # return sample 
}


# -----------------------------------------------------------------------------
#### Full simulation ####
# -----------------------------------------------------------------------------

# number of observation points
nTvec <- seq(10, 100, by = 10)
# t <- commandArgs(trailingOnly = TRUE) 
# t <- as.numeric(t)
t <- 1
nT <- nTvec[t]


# number of variables
Pvec <- seq(5, 50, by = 5)

pb <- txtProgressBar(0, length(Pvec), 0, style = 3)
comp_time_vec <- numeric(length(Pvec))

for(k in seq_along(Pvec)){
  P <- Pvec[k]
  
  cat("=== Generating datasets ===\n")
  ds_smooth <- gen_smooth(p = P, T = nT)
  # write_rds(ds_smooth, paste0("Manuscripts/Data/sim_smooth_p", P, "_T", nT, ".rds"))
  
  cat("\n=== Fit Dynamic MDS ===\n")
  # search grid
  lambdas <- seq(0, 10, length.out = 100)
  # search 
  comp_time <- system.time({
    sim_obs_cors <- lapply(ds_smooth$obs, cor, method="spearman")
    sweep_smooth <- lambda_sweep(sim_obs_cors, lambdas)
  })
  
  comp_time_vec[k] <- comp_time[3] 
  
  # save simulation output
  # write_rds(sweep_smooth, paste0("Manuscripts/Output/out_smooth_p", P, "_T", nT, ".rds"))
  
  setTxtProgressBar(pb, k)
  
}

close(pb)


# save computation time
# comp_time_file <- "Manuscripts/Output/comp_time.csv"
# comp_time_df <- data.frame(T = nT, P = Pvec,  comp_time = comp_time_vec)
# if (file.exists(comp_time_file)) {
#   write.table(comp_time_df, file = comp_time_file, append = TRUE, 
#               sep       = ",",
#               row.names = FALSE,
#               col.names = !file.exists(comp_time_file))
# } else {
#   write.csv(comp_time_df, comp_time_file, row.names = FALSE)
# }




