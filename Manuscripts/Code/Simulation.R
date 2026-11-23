# -----------------------------------------------------------------------------
#### Set up ####
# -----------------------------------------------------------------------------


library(MASS)      # mvrnorm
library(Matrix)    # nearPD
library(tidyverse)

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


#### L curves ####
file_names <- list.files("Manuscripts/Output")[-1]

best_lam_df <- list()
plot_list <- list()
output_list <- list()

for(i in seq_along(file_names)){
  
  Pi <- as.numeric(str_extract(file_names[i], "(?<=p)\\d+"))
  Ti <- as.numeric(str_extract(file_names[i], "(?<=T)\\d+"))
  
  output_i <- read_rds(paste0("Manuscripts/Output/", file_names[i]))
  output_i$P <- Pi
  output_i$T <- Ti
  
  output_list[[i]] <- output_i
  
  
  # best lambda
  # stress0 <- output_i$stress[output_i$lambda==0]
  # movement0 <- output_i$movement[output_i$lambda==0]
  # output_norm <- output_i %>%
  #   mutate(stress = stress/stress0,
  #          movement = movement/movement0)
  
  # best_lam <- data.frame(
  #   method = c("Menger", "Max dist"),
  #   lambda = c(lcurve_corner_menger(output_i, plot = F)$lambda_star,
  #   # move_ave = lcurve_corner_ma(output_i)$lambda_star,
  #   # splines = lcorner_splines(output_i, plot = F, spar = 0.5)$lambda_star,
  #              lcurve_corner_dist(output_i)$lambda_star)
  # )
  # 
  # best_lam <- lcurve_corner_angle(output_i)$lambda_star
  best_lam <- c(Menger = lcurve_corner_menger(output_i, plot = F)$lambda_star, 
                "Max dist" = lcurve_corner_dist(output_i)$lambda_star)
  best_lam_df[[i]] <- c(P = Pi, T = Ti, best_lam)
  # points <- output_i %>% filter(lambda %in% best_lam$lambda) %>% 
  #   left_join(best_lam)
  
}

p_lev <- paste("P =", c(5, 10, 15, 25, 50))
t_lev <- paste("T =", c(10, 20, 30, 50, 100))

all_output <- bind_rows(output_list) %>%
  filter(T %in% c(10, 20, 30, 50, 100) & P %in% c(5, 10, 15, 25, 50)) %>%
  mutate(P = paste("P =", P), T = paste("T =", T)) %>%
  mutate(P = factor(P, levels = p_lev), T = factor(T, levels = t_lev))
all_lambda <- bind_rows(best_lam_df) %>%
  filter(T %in% c(10, 20, 30, 50, 100) & P %in% c(5, 10, 15, 25, 50)) %>%
  mutate(P = paste("P =", P), T = paste("T =", T)) %>%
  pivot_longer(3:4, names_to = "method", values_to = "lambda") %>%
  left_join(all_output) %>%
  mutate(P = factor(P, levels = p_lev), T = factor(T, levels = t_lev))

p <- all_output %>%
    ggplot(aes(x=log(movement), y=log(stress)))+
    geom_point(size = 0.5) +
    geom_line() +
    geom_point(data = all_lambda, aes(x = log(movement), y = log(stress)),
               size = 2, col = "red") +
    ggrepel::geom_text_repel(
      data          = all_lambda,
      aes(x = log(movement), y = log(stress), label = method),
      color         = "red",
      size          = 5,
      lineheight    = 0.85,
      arrow         = arrow(length = unit(0.01, "npc"), type = "closed"),
      min.segment.length = 0,   # always draw the arrow, even if close
      segment.colour = "red",
      segment.size = 0.3,
      force         = 0.5,
      force_pull = 1,
      max.overlaps  = Inf,
      # push labels toward the upper-right relative to their point
      direction     = "both",
      nudge_x       = 0,
      nudge_y       = 0,
      hjust         = 0,
      vjust         = 0
    ) +
    # force equal visual scaling on both axes so the right angle
    # is actually rendered as a right angle
    # coord_fixed(ratio = 1) +
    labs(x = "Penalty (log-scaled)", y = "Stress (log-scaled)")+
    theme_minimal()+
  theme(axis.title = element_text(size = 20), strip.text = element_text(size = 20))+
  # facet_grid(rows = vars(P), cols = vars(T), scales = "free")
  # facet_wrap(~P+T, scales="free", labeller = label_both)
  ggh4x::facet_grid2(rows = vars(P), cols = vars(T),
              scales = "free", independent = "all")
p
ggsave("Manuscripts/Figures/sim_lcurves.pdf", p, 
       width = 24, height = 20, units = "in", dpi = 300)
  
  



per_page <- 24
n_plots  <- length(plot_list)                     # 110
n_pages  <- ceiling(n_plots / per_page)           # 5

for(i in 1:n_pages){
  idx <- ((i - 1) * per_page + 1):min(i * per_page, n_plots)
  p <- ggpubr::ggarrange(
    plotlist = plot_list[idx],
    ncol = 4,
    nrow = 6,
    labels = "AUTO",       # Adds uppercase A, B, C... labels
    common.legend = TRUE,  # Uses the legend from the first plot
    legend = "bottom"
  )
  
}

per_page <- 24
n_plots  <- length(plot_list)                     # 110
n_pages  <- ceiling(n_plots / per_page)           # 5



lapply(seq_len(n_pages), function(i) {
  idx <- ((i - 1) * per_page + 1):min(i * per_page, n_plots)
  
  p <- ggpubr::ggarrange(
    plotlist      = plot_list[idx],
    ncol          = 4,
    nrow          = 6,
    common.legend = TRUE,
    legend        = "bottom"
  )
  
  ggsave(paste0("Manuscripts/Figures/sim_lcurves", i, ".pdf"), p, 
         width = 14, height = 20, units = "in", dpi = 300)
  
})


