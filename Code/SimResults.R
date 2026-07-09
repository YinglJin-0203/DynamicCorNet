
rm(list=ls())

library(tidyverse)
library(here)
library(igraph)
theme_set(theme_minimal())

source("Simulation/Code/SimHelpers.r")
source("Simulation/Code/SimFit.R")

#### Data ####

# sparse example
exp_data <- read_rds("Simulation/Data/sim_smooth_p50_T100.rds")
# correlation between X1 and X5
obs_cor <- sapply(exp_data$obs, function(x){cor(x[, 1], x[, 30], method = "spearman")})
tru_cor <- sapply(exp_data$true_corrs, function(x)x[1, 30])
png("Simulation/SimulationFigs/example_cor_large.png")
data.frame(time = 1:length(tru_cor), obs_cor, tru_cor) %>%
  pivot_longer(2:3) %>%
  ggplot(aes(x=time, y=value, col=name, group=name))+
  geom_line()+
  geom_point()+
  labs(x="time",y="cor", col=" ", title = "P = 50, T = 100")
dev.off()

#### Fit results ####

# compare best lambda selectiong methods?

file_names <- list.files("Simulation/Output")[-1]

best_lam_df <- list()

plot_list <- list()

for(i in seq_along(file_names)){
  
  Pi <- as.numeric(str_extract(file_names[i], "(?<=p)\\d+"))
  Ti <- as.numeric(str_extract(file_names[i], "(?<=T)\\d+"))
  
  output_i <- read_rds(paste0("Simulation/Output/", file_names[i]))
  
  # best lambda
  # stress0 <- output_i$stress[output_i$lambda==0]
  # movement0 <- output_i$movement[output_i$lambda==0]
  # output_norm <- output_i %>%
  #   mutate(stress = stress/stress0,
  #          movement = movement/movement0)
  
  best_lam <- c(
    menger = menger_corner(output_i, plot = F)$lambda_star,
    move_ave = lcurve_corner_ma(output_i)$lambda_star,
    splines = lcorner_splines(output_i, plot = F, spar = 0.5)$lambda_star,
    perp_dist = lcurve_corner_dist(output_i)$lambda_star)
          
  
  # best_lam <- lcurve_corner_angle(output_i)$lambda_star
  best_lam_df[[i]] <- c(P = Pi, T = Ti, best_lam)
  
  # plot
  plot_list[[i]] <- output_i %>%
    pivot_longer(2:4) %>%
    group_by(name) %>%
    ggplot(aes(x=lambda, y=value, col=name, group=name))+
    geom_point(size = 0.5)+
    geom_line()+
    geom_vline(xintercept = best_lam)+
    facet_wrap(~name, nrow = 1, scales = "free")+
    labs(x = "lambda", y=" ", title = paste("P =", Pi, ", T = ", Ti),
         col = " ")

  
}

png("Simulation/SimulationFigs/pick_lambda_compare.png", height=500, width = 1500)
bind_rows(best_lam_df) %>% 
  pivot_longer(3:6) %>%
  ggplot(aes(x=value))+
  geom_histogram(bins=50)+
  facet_wrap(~name)
dev.off()

png("Simulation/SimulationFigs/stress_move_lambda_random.png", height=500, width = 900)
do.call(ggpubr::ggarrange, c(plotlist = plot_list[sample(length(plot_list), 4)], 
                             nrow = 2, ncol=2, common.legend=T))
dev.off()

plot_list[[2]]


png("Simulation/SimulationFigs/best_lambda_menger_tile.png")

best_lam_mat <- as.data.frame(best_lam_mat[, 1:3])
colnames(best_lam_mat) <- c("P", "T", "best_lam")
best_lam_mat %>%
  filter(T != 15) %>%
  ggplot(aes(x=T, y=P, fill=best_lam))+
  geom_tile()+
  labs(title = "Menger", fille = "lambda")
dev.off()

plot_list[[60]]

#### network plots ####

# method 2: L curve normalized angle

par(mfrow = c(1, 1))
exp_output <- read_rds("Simulation/Output/out_smooth_p15_T10.rds")
plot(log(exp_output$movement), log(exp_output$stress))
best_lam <- menger_corner(exp_output)$lambda_star # best lambda: 9.9

exp_data <- read_rds("Simulation/Data/sim_smooth_p15_T10.rds")
obs_cor_p5 <- lapply(exp_data$obs, cor, method="spearman")
dyn_fit <- dyn_mds(obs_corrs = obs_cor_p5,
        lambda = best_lam, d = 2)
layout <- dyn_fit$embeddings

png("Simulation/SimulationFigs/dynamic_network_p15.png", height=500, width = 1500)
par(mfrow = c(2, 5))
for(i in seq_along(layout)){
  layout_i <- layout[[i]]
  
  # edges
  cor_mat_i <- obs_cor_p5[[i]]
  rownames(cor_mat_i) <- colnames(cor_mat_i) <- paste0("X", 1:15)
  edges <- which(abs(cor_mat_i) > 0.5 & upper.tri(cor_mat_i), arr.ind = T)
  from      <- rownames(cor_mat_i)[edges[, 1]]
  to        <- colnames(cor_mat_i)[edges[, 2]]
  wt <- cor_mat_i[edges]
  
  # initialize 
  net_i <- igraph::graph_from_data_frame(
    d = data.frame(from=from, to=to, weight = wt),
    directed = F, vertices = data.frame(name = rownames(cor_mat_i))
  )
  
  # visual elements of edges
  E(net_i)$width <- abs(E(net_i)$weight) * 8          # scale thickness to |value|
  E(net_i)$color <- ifelse(E(net_i)$weight > 0, "steelblue", "tomato")  # sign → color
  

  plot(net_i, layout = layout_i, main = paste0("t = ", i), vertex.size = 30,
       vertex.label.cex = 1.5,
       vertex.color = adjustcolor("lightblue", alpha.f = 0.5),
       vertex.frame.color = "lightblue",
       edge.curved = 0.2, 
       margin = c(0, 0, 0, 0))
  
}

legend("bottomleft",
       legend = c("Positive", "Negative"),
       col    = c("steelblue", "tomato"),
       lwd    = 3, bty = "n"
)
dev.off()
## This led to very, very small change in graphical layout that are not even noticeable
## The graph is basically still
## this happens even after normalizing with lambda = 0
## L curve: maximum curvature in log(Stress)~log(movement)
## may not have a clear elbow - latch onto the far right end
## unequal axis scaling can distort curvature
## discrete axis artifacts

exp_output2 <- read_rds("Simulation/Output/out_smooth_p35.rds")
par(mfrow = c(1, 1))
stress0 <- exp_output$stress[exp_output$lambda==0]
movement0 <- exp_output$movement[exp_output$lambda=0]
exp_output_norm <- exp_output %>%
  mutate(
    stress_norm = 
  )
plot(exp_output$movement, exp_output$stress)
points(log(exp_output2$movement), log(exp_output2$stress), col="red")

#### Computation time ####
comp_time <- read.csv("Simulation/Output/comp_time.csv")

range(comp_time$comp_time) # 8 second to 2 minutes
head(comp_time)
ggplot(comp_time, aes(x=T, y=comp_time, col=P, group = P))+
  geom_point()+
  geom_line()
