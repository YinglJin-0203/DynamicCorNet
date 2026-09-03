# -----------------------------------------------------------------------------
# Whole simulation
# -----------------------------------------------------------------------------

rm(list=ls())

#### Set up  ####
source("Simulation/Code/SimHelpers.r")
source("Simulation/Code/SimData.R")
source("Simulation/Code/SimFit.R")

library(tidyverse)

# number of observation points
nTvec <- seq(10, 100, by = 10)
t <- commandArgs(trailingOnly = TRUE)
t <- as.numeric(t)
nT <- nTvec[t]

nT <- 10

# number of variables
Pvec <- seq(5, 50, by = 5)


pb <- txtProgressBar(0, length(Pvec), 0, style = 3)
#comp_time_vec <- numeric(length(Pvec))

P <- 10
for(k in seq_along(Pvec)){
  P <- Pvec[k]
  
  cat("=== Generating datasets ===\n")
  ds_smooth <- gen_smooth(p = P, T = nT)
  write_rds(ds_smooth, paste0("Simulation/Data/sim_smooth_p", P, "_T", nT, ".rds"))
  # write_rds(ds_smooth, paste0("ata/sim_smooth_p", P, "_T", nT, ".rds"))
  # ds_noise  <- gen_high_noise()
  
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
  write_rds(sweep_smooth, paste0("Simulation/Output/out_smooth_p", P, "_T", nT, ".rds"))
  
  setTxtProgressBar(pb, k)

}

close(pb)


# save computation time
# save computation time
comp_time_file <- "Simulation/Output/comp_time.csv"
comp_time_df <- data.frame(T = nT, P = Pvec,  comp_time = comp_time_vec)
if (file.exists(comp_time_file)) {
  write.table(comp_time_df, file = comp_time_file, append = TRUE, 
              sep       = ",",
              row.names = FALSE,
              col.names = !file.exists(comp_time_file))
} else {
  write.csv(comp_time_df, comp_time_file, row.names = FALSE)
}


# comp_time_vec

# sweep_smooth %>%
#   pivot_longer(2:4) %>%
#   ggplot(aes(x=lambda, y=value, col=name, group=name))+
#   geom_line()+
#   geom_point(size = 0.5)
# 
# which.min(sweep_smooth$objective)
# corner_smooth <- lcurve_corner(sweep_smooth)
# 
# 
# cat(sprintf("\nSelected lambda* = %.4f\n", corner_smooth$lambda_star))
# 
# cat("\n=== Scenario 2: High Noise / Low Signal ===\n")
# sweep_noise <- lambda_sweep(ds_noise, lambdas)
# corner_noise <- lcurve_corner(sweep_noise)
# cat(sprintf("\nSelected lambda* = %.4f\n", corner_noise$lambda_star))
# 
# # Plots
# plot_sweep(sweep_smooth, "Slow/Smooth Dynamics", corner_smooth$lambda_star)
# plot_sweep(sweep_noise,  "High Noise / Low Signal", corner_noise$lambda_star)
# 
# # Movement ratio diagnostic
# R_smooth <- sweep_smooth$movement[which.min(abs(sweep_smooth$lambda - corner_smooth$lambda_star))] /
#   sweep_smooth$movement[which.min(sweep_smooth$lambda)]
# R_noise  <- sweep_noise$movement[which.min(abs(sweep_noise$lambda  - corner_noise$lambda_star))]  /
#   sweep_noise$movement[which.min(sweep_noise$lambda)]
# 
# cat(sprintf("\nMovement ratio R(lambda*): Smooth = %.2f | High Noise = %.2f\n",
#             R_smooth, R_noise))
