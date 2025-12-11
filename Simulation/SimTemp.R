# This script is used to simulate multivariate temporal data and 
# implement dynamic and splines MDS on them
# generated data and graph coordinates are saved

set.seed(123)  # for reproducibility

library(tidyverse)
library(here)
library(smacof)
library(splines2)
library(mgcv)

# function to simulate dataset
source(here("helpers/SimGroupT.R"))

# functions to generate graph
source(here("helpers/AdjacencyMat.R"))
source(here("helpers/DynamicMDS.R"))
source(here("helpers/SplinesMDS.R"))


#### Scalars ####
P <- 15 # total number of variables
pGroup  <- rep(5, 3) # variable in each group
G <- 3 # number of groups
N <- 100 # sample size

#### varying scalar: measure density ####
# K <- commandArgs(trailingOnly = TRUE)
# K <- as.numeric(K)
# nt <- 10+10*(K-1)
nt <- 50
tvec <- seq(0, 1, length.out = nt)

#### Generate data ####

##### constant between- and within-group correlation #####

# # group level
# cor_group <- matrix(0.5, G, G)
# diag(cor_group) <- 1
# # cor_group
# 
# # group-var level
# rho <- c(0.3, 0.5, 0.7)
# cor_group_var <- list()
# for(g in seq_along(rho)){
#   cor_g <- matrix(rho[g], nrow=pGroup[g], ncol=pGroup[g])
#   diag(cor_g) <- 1
#   cor_group_var[[g]] <- cor_g
# }
# cor_group_var <- Matrix::bdiag(cor_group_var)
# # cor_group_var
# 
# # generate data
# df_sim_temp1 <- lapply(tvec, GenDataT, N=100, P=15, G=3, pGroup=rep(5, 3),  cor_g = cor_group, 
#                        cor_pg = cor_group_var, sigma=1, f=0, id = 1:N)
# df_sim_temp1 <- bind_rows(df_sim_temp1)
# 
# # save data
# saveRDS(df_sim_temp1, file = paste0("data/SimData/df_sim_temp1_t", nt, ".rds"))
# 

##### changing within-group and constant between-group correlation #####

# group level
cor_group <- matrix(0.5, G, G)
diag(cor_group) <- 1
# cor_group

# within group
cor_group_var <- list()

rho_max <- 0.7
for(t in seq_along(tvec)){
  
  # group 1 not changing
  cor_g1 <- matrix(0.3, nrow=pGroup[1], ncol=pGroup[1])
  diag(cor_g1) <- 1
  
  # group 2 increasing
  cor_g2 <- matrix((rho_max/(nt-1))*(t-1), pGroup[2], pGroup[2])
  diag(cor_g2) <- 1
  
  # group 3 decreasing
  cor_g3 <- 0.7-cor_g2
  diag(cor_g3) <- 1
  
  cor_t <- Matrix::bdiag(cor_g1, cor_g2, cor_g3)
  cor_group_var[[t]] <- cor_t
}

# cor_group_var

# generate
df_sim_temp2 <- lapply(1:nt, 
                       function(t){
                         GenDataT(N=100, t = tvec[t], P = 15, G = 3, pGroup = rep(5, 3),
                                  cor_g = cor_group, cor_pg = cor_group_var[[t]],
                                  sigma=1, f=0, id = 1:N)
                       })
df_sim_temp2 <- bind_rows(df_sim_temp2)

# save data
# saveRDS(df_sim_temp2, file = paste0("data/SimData/df_sim_temp2_t", nt, ".rds"))

#### Calculated coordinates ####

##### Dynamic MDS #####


# step 1: adjacency matrix
adj_mat1 <- GetAdjMat(df_sim_temp2 %>% dplyr::select(-id), mds_type = "Dynamic")
# graph_list1 <- graph_dyn_net(adj_mat1, mds_type = "Dynamic")

# step 2: calculate layout
# time1 <- system.time({
prof1 <- profvis::profvis({
  coord_list1 <- DynamicMDS(adj_mat1, lambda = 5)
})

#### Calculate coordinates for Splines MDS ####
# step 1: adjacency matrix
adj_mat2 <- GetAdjMat(df_sim_temp2 %>% dplyr::select(-id), mds_type = "Splines")
# graph_list2 <- graph_dyn_net(adj_mat2, mds_type = "Splines")

# step 2: calculate layout
# time2 <- system.time({
prof2 <- profvis::profvis({
  coord_list2 <- SplinesMDS(adj_mat2, lambda = 10, K = 30, P = 15, tvec = 1:nt)
})

#### Output ####
# computation time
out_vec <- data.frame(Density = nt, Dyn = time1[3], Spl = time2[3])
# print(out_vec)

if(K==1){
  write.csv(out_vec, file = "Output/SimTemp/DenseCompTime.csv")
} else{
  write.table(out_vec, file = "Output/SimTemp/DenseCompTime.csv", sep = ",", append = T, col.names = F)
}

# coordinates
filename <- paste0("Output/SimTemp/coord_t", nt, ".RData")
save(coord_list1, coord_list2, file = filename)
