rm(list = ls())

# This script is to investigae
# 1. Whether magnitude chanage can be reflected
# 2. How to determine penalization parameter
# 3. How does initialization affect the algorithm

#### Set up ####
set.seed(123)  # for reproducibility

# packages
library(tidyverse)
theme_set(theme_minimal())
library(here)
library(smacof)
library(splines2)
library(mgcv)
# library(Rcpp)
# library(igraph)

source(here("Code/DynDissmMat.R"))
source(here("Code/Helpers/DynMDSHelpers.R"))
source(here("Code/DynamicMDS.R"))


#### Generate/load data ####

# # load more complicated simulated data
# df_sim <- read_rds("Simulation/Data/df_sim_vary_cor_t10.rds")
# 
# # scalars
# P <- ncol(df_sim) - 2 # total number of variables
# N <- length(unique(df_sim$id)) # sample size
# nt <- length(unique(df_sim$time)) # number of measurements
# tgrid <- sort(unique(df_sim$time))

# or generate simple data

## investigate the effect of dimension
Pvec <- seq(5, 50, by = 5) # investigate the effect of variable 
# P <- 5
N <- 100 # sample size
nt <- 3 # time grid
cor_vec <- c(0.2, 0.5, 0.9)


# value to parallel over
k <- commandArgs(trailingOnly = TRUE)
p <- as.numeric(k)

# same within-group correlation
# no between group correlation

# gen_data_list <- list()
# 
# for(p in seq_along(Pvec)){
  
P <- Pvec[p]
  
df_sim <- expand.grid(id = 1:N, time = 1:nt)
df_sim[ , paste0("X", 1:P)]<- NA

for(t in 1:nt){
  cor_t <- cor_vec[t]
  cor_mat_t <- diag(1, P, P)
  cor_mat_t[row(cor_mat_t)!=col(cor_mat_t)] <- cor_t
  
  df_sim[df_sim$time==t, paste0("X", 1:P)] <- MASS::mvrnorm(n=N, mu=rep(0, P), Sigma = cor_mat_t)
}
  
#   gen_data_list[[p]] <- df_sim
#   
# }
# 
# lapply(gen_data_list, dim)

# check
# gen_data_list[[4]] %>% group_by(time) %>% select(-id) %>% group_map(~cor(.x, method = "spearman"))

#### Model fit ####

# parameters to investigate: P

# containers
# L <- length(Pvec)
# 
# lambda_center <- numeric(L)
# lambda_optimal <- numeric(L)
# converge0 <- numeric(L)
# converge_num <- numeric(L)
# 
# comp_time <- numeric(L)
# min_loss_vec <- numeric(L)
# loss_list <- list()
# stress_list <- list()
# kendall_list <- list()
# stability_list <- list()
# 
# pb <- txtProgressBar(0, L, 0, style = 3)

# for(p in seq_along(gen_data_list)){
  
  # dissimilarity
  # P <- Pvec[p]
dis_list_dyn <- DynDissimMat(df_sim %>% select(-id))
  
  ##### initialization (scale tuned to correlation) #####
t1 <- Sys.time()
init_par <- runif(sum(sapply(dis_list_dyn, ncol))* 2, -1/sqrt(2), 1/sqrt(2))

# determine the range of lambda:
# fit non-penalized MDS, use the stress/penalty as center
coords0 <- stats::optim(
  par = init_par,
  fn = stress_DynMDS,
  diss_list = dis_list_dyn,
  P = sapply(dis_list_dyn, ncol),
  ndim = 2, 
  lambda = 0, 
  scale = T,
  method = "BFGS",
  control = list(maxit = 500)
)

  # converge0[p] <- coords0$con
coords0 <- reshape_configs(coords0$par, sapply(dis_list_dyn, ncol), 2, nt)
  
  # stress
stress0 <- mapply(function(coord, diss) {
  diss_vec <- diss[lower.tri(diss, diag = F)]
  
  sum((diss_vec - as.vector(dist(coord, diag = F, upper=F)))^2) /
    sum(diss_vec^2)
}, coords0, dis_list_dyn)

  
  # stability penalization
stability0 <- sum(sapply(Map(`-`, coords0[-1], coords0[-length(coords0)]), function(x)sum(x^2)))
  # center of search grid
lambda_center_p <- sum(stress0)/stability0
    # lambda_center[p] <- lambda_center_p
    
  # search grid
lambda_grid <- outer(lambda_center_p, 10^seq(-5, 5, by = 1), "*")

  ##### grid seach #####
    
ngrid <- length(lambda_grid)
    # containers
stress_mat <- matrix(NA, nt, ngrid)
kendall_mat <- matrix(NA, nt, ngrid)
stability_vec <- numeric(ngrid)
loss_vec <- numeric(ngrid)
converge_vec <- numeric(ngrid)
    
for(l in 1:ngrid){
  
  lambda_l <- lambda_grid[l]
  final_coords <- stats::optim(
        par = init_par,
        fn = stress_DynMDS,
        diss_list = dis_list_dyn,
        P = sapply(dis_list_dyn, ncol),
        ndim = 2, 
        lambda = lambda_l, 
        scale = T,
        method = "BFGS",
        control = list(maxit = 500)
      )
    
    loss_vec[l] <- final_coords$value
    converge_vec[l] <- final_coords$convergence
    
    # coordinates
    coords <- reshape_configs(vec = final_coords$par, 
                              P = sapply(dis_list_dyn, ncol),
                              ndim = 2, Tmax = nt)
    # stress
    stress_mat[ , l] <- mapply(function(coord, diss) {
      diss_vec <- diss[lower.tri(diss, diag = F)]
      
      sum((diss_vec - as.vector(dist(coord, diag = F, upper=F)))^2) /
        sum(diss_vec^2)
    }, coords, dis_list_dyn)
    
    # kendall
    kendall_mat[ ,l] <- mapply(function(coord, diss) {
      diss_vec <- diss[lower.tri(diss, diag = F)]
      dist_vec <-  as.vector(dist(coord, diag = F, upper=F))
      cor(diss_vec, dist_vec, method = "kendall")
    }, coords, dis_list_dyn)
    
    # stability penalization
    stability_vec[l] <- sum(sapply(Map(`-`, coords[-1], coords[-length(coords)]), function(x)sum(x^2)))
}
t2 <- Sys.time()
    # comp_time[p] <- difftime(t2, t1, units = "secs") # time diff in seconds
comp_time <- difftime(t2, t1, units = "secs") # time diff in seconds
    # setTxtProgressBar(pb, p)
    
    # stress_list[[p]] <- stress_mat
    # kendall_list[[p]] <- kendall_mat
    # stability_list[[p]] <- stability_vec
    # loss_list[[p]] <- loss_vec
    
    # best lambda
    # lambda_optimal[p] <- lambda_grid[which.min(loss_vec)]
    # min_loss_vec[p] <- min(loss_vec)
    # 
    # converge
    # converge_num[p] <- sum(converge_vec)


# }

# close(pb)

#### save results ####

out_p <- list(
  stress = stress_mat,
  kendall = kendall_mat,
  stability = stability_vec,
  loss = loss_vec,
  converge = converge_vec,
  lambda_grid = lambda_grid
)

filename <- paste0("Simulation/Output/output_P", Pvec[p], ".rds")
saveRDS(out_p, file = filename)

#### examine results ####

# coverge
# converge0
# converge_num # there is a computational problem at P = 15? 
# 
# # initial
# lambda_center # didn't change much
# 
# # lambda selection
# lambda_optimal # changed a lot
# png("Simulation/SimulationFigs/lambda_p_optimal_lambda.png")
# plot(Pvec, log(lambda_optimal), type = "b", 
#      xlab = " ", ylab = " ", main = "Log(lambda) ~ P")
# dev.off()
# 
# # minimum loss
# min_loss_vec
# png("Simulation/SimulationFigs/lambda_p_optimal_loss.png")
# plot(Pvec[1:5], min_loss_vec[1:5],
#      type = "b", 
#      xlab = " ", ylab = " ", main = "Minimum loss ~ P") # increase with P
# dev.off()
# 
# # computation time
# comp_time
# png("Simulation/SimulationFigs/lambda_p_comp_time.png")
# plot(Pvec, comp_time, type = "b",
#      xlab = " ", ylab = " ", main = "Computation (seconds) ~ P")
# dev.off()
# 
# # stress: change in rougly the same shape across grid
# png("Simulation/SimulationFigs/lambda_p_stress.png")
# sapply(stress_list, colSums) %>%
#   data.frame() %>%
#   mutate(lambda_diff = 10^seq(-5, 5, by = 1)) %>%
#   pivot_longer(X1:X4) %>%
#   mutate(name=factor(name, levels = paste0("X", 1:4), 
#                      labels = paste0("P = ", Pvec[1:4]))) %>%
#   ggplot(aes(x=lambda_diff, y=value, group=name, col=name))+
#   geom_line()+
#   geom_point()+
#   scale_x_log10()+
#   labs(x = " ", y=" ", title = "Stress ~ diff of log(lambda)", col = "")
# dev.off()
# 
# # kendall: Also roughly the same shape, but less P is more variable
# sapply(kendall_list, colSums) %>%
#   data.frame() %>%
#   mutate(lambda_diff = 10^seq(-5, 5, by = 1)) %>%
#   pivot_longer(X1:X4) %>%
#   mutate(name=factor(name, levels = paste0("X", 1:4), 
#                      labels = paste0("P = ", Pvec[1:4]))) %>%
#   ggplot(aes(x=lambda_diff, y=value, group=name, col=name))+
#   geom_line()+
#   scale_x_log10()
# 
# # loss
# loss_list
# png("Simulation/SimulationFigs/lambda_p_loss.png")
# do.call(rbind, loss_list) %>% t() %>%
#   as.data.frame() %>%
#   mutate(lambda_diff = 10^seq(-5, 5, by = 1)) %>%
#   pivot_longer(V1:V4) %>%
#   mutate(name=factor(name, levels = paste0("V", 1:4), 
#                      labels = paste0("P = ", Pvec))) %>%
#   ggplot(aes(x=lambda_diff, y=value, group=name, col=name))+
#   geom_line()+
#   geom_point()+
#   scale_x_log10()+
#   labs(x = " ", y=" ", title = "Loss ~ diff of log(lambda)", col = "")
# dev.off()
# 
# # stability: also same shape, increase with P
# # I wonder if I penalize "Average change" would it be better? 
# png("Simulation/SimulationFigs/lambda_p_stability.png")
# do.call(rbind, stability_list) %>% t() %>%
#   as.data.frame() %>%
#   mutate(lambda_diff = 10^seq(-5, 5, by = 1)) %>%
#   pivot_longer(V1:V4) %>%
#   mutate(name=factor(name, levels = paste0("V", 1:4), 
#                      labels = paste0("P = ", Pvec))) %>%
#   ggplot(aes(x=lambda_diff, y=value, group=name, col=name))+
#   geom_line()+
#   geom_point()+
#   scale_x_log10()+
#   labs(x = " ", y=" ", title = "Stability ~ diff of log(lambda)", col = "")
# dev.off()
# 
# # plot
# 
# plot_list <- list()
# 
# for(i in seq_along(Pvec[1:5])){
# 
#   xgrid <- 10^seq(-5, 5, by = 1)
#   stress_i <- colSums(stress_list[[i]])
#   stability_i <- stability_list[[i]]
#   loss_i <- loss_list[[i]]
#   
#   id_min_loss <- which.min(loss_list[[i]]) # minimum lambda
#   id_eq_wt <- which.min(abs(stress_i-stability_i)) # if stress = stability
#   
#   # plot
#   plot_list[[i]]<-data.frame(
#     lambda_diff = rep(xgrid, 3),
#     value = c(stress_i, stability_i, loss_i),
#     name = rep(c("stress","stability", "loss"), each = ngrid)
#   ) %>%
#     ggplot(aes(x=lambda_diff, y=value, col=name, group=name))+
#     geom_line()+
#     geom_point()+
#     geom_vline(xintercept = xgrid[id_min_loss])+
#     geom_vline(xintercept = xgrid[id_eq_wt])+
#     scale_x_log10()+
#     labs(x = "Diff of log(lambda)", y=" ", title = paste("P =", Pvec[i]),
#          col = " ", linetype = " ")
#  
# }
# 
# 
# plot_list[[1]]
# plot_list[[2]]
# plot_list[[3]]
# plot_list[[4]]
# 
# for(i in 1:4){
#   filename <- paste0("Simulation/SimulationFigs/lambda_p_select_P", Pvec[i], ".png")
#   png(filename)
#   print(plot_list[[i]])
#   dev.off()
# }
# 
# png("Simulation/SimulationFigs/lambda_p_select.png", width = 15, height = 3)
# gridExtra::grid.arrange(grobs = plot_list, nrow = 1)
# dev.off()
