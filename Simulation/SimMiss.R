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
source(here("Code/SimGroupT.R"))

# functions to generate graph
source(here("Code/AdjacencyMat.R"))
source(here("Code/Helpers/DynMDSHelpers.R"))
source(here("Code/DynamicMDS.R"))
source(here("Code/Helpers/SplMDSHelpers.R"))
source(here("Code/SplinesMDS.R"))


#### Scalars ####
P <- 15 # total number of variables
pGroup  <- rep(5, 3) # variable in each group
G <- 3 # number of groups
N <- 100 # sample size

#### varying scalar: measure density ####
# nt <- 10
nt <- 100
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
df_sample <- df_sim_temp2
write.csv(df_sample, file="../SampleData/df_sample.csv", row.names = F)
##### introduce irregular missing #####

# group 1 random missing at 1/2 of the time

t_miss <- sample(tvec, size = nt*0.5, replace = F)
df_sim_temp2[df_sim_temp2$time %in% t_miss, paste0("X", 1:5)] <- NA

#### Calculated coordinates ####

##### Dynamic MDS #####

# step 1: adjacency matrix
adj_mat1 <- GetAdjMat(df_sim_temp2 %>% dplyr::select(-id), mds_type = "Dynamic")
lapply(adj_mat1, dim)

# step 2: calculate layout
time1 <- system.time({
  # prof1 <- profvis::profvis({
  coord_list1 <- DynamicMDS(adj_mat1, lambda = 10)
})
# prof1 # time: 16370 when nt = 10; 

#### Calculate coordinates for Splines MDS ####
# step 1: adjacency matrix
adj_mat2 <- GetAdjMat(df_sim_temp2 %>% dplyr::select(-id), mds_type = "Splines")
lapply(adj_mat2, dim)

# step 2: calculate layout
time2 <- system.time({
  # prof2 <- profvis::profvis({
  coord_list2 <- SplinesMDS(adj_mat2, lambda = 7, P = 15, tvec = tvec)
})
# prof2 # time: 26350


#### Output ####

# temporal network figures
library(igraph)
library(RColorBrewer)
display.brewer.pal(7, "Set2")
color_vec <- brewer.pal(3, "Set2")
library(magick)
library(tidyverse)
theme_set(theme_minimal())

# if needed, load output
load("Simulation/coord_miss_t200.RData")
tvec <- seq(0, 1, length.out = 200)
P <- 15

coord_list1[c(1, 51, 101, 155, 192)]

t <- 192
for(t in c(1, 51, 101, 155)){
  
  coord_t1 <- coord_list1[[t]] # dynamic
  coord_t2 <- coord_list2[[t]] # splines
  
  # plot dynamic
  graph_t <- make_empty_graph(n = nrow(coord_t1), directed = FALSE)
  if(nrow(coord_t1)==10){
    V(graph_t)$color <- color_vec[rep(2:3, each = 5)] # color by true group
    V(graph_t)$label = 6:15
  } else {
    V(graph_t)$color <- color_vec[rep(1:3, each = 5)] # color by true group
    V(graph_t)$label = 1:15
    }
  pic1 <- paste0("images/SimFigures/DynMDS_miss2_t", t, ".jpeg")
  jpeg(filename = pic1, height = 500, width = 500)
  plot(graph_t,
       layout = as.matrix(coord_t1),
       # vertex.frame.color=ifelse(is.na(V(graph_t)$color), "grey", NA),
       vertex.frame.color=NA,
       vertex.label.cex=1,
       vertex.size = 20, 
       vertex.color = V(graph_t)$color, 
       margin = 0, main = paste0("Dynamic, Time = ", round(tvec[t], 2)))
  dev.off()
  # fig_list1[[t]] <- image_read(pic1)
  
  # plot spline
  graph_t <- make_empty_graph(n = P, directed = FALSE)
  V(graph_t)$color <- color_vec[rep(1:3, each = 5)] # color by true group
  V(graph_t)$label = 1:15
  if(nrow(coord_t1)==10){
    V(graph_t)$color[1:5] <- NA # color by true group
  }
  pic2 <- paste0("images/SimFigures/SplMDS_miss2_t", t, ".jpeg")
  jpeg(filename = pic2, height = 500, width = 500)
  plot(graph_t,
       layout = as.matrix(coord_t2),
       vertex.frame.color=ifelse(is.na(V(graph_t)$color), "grey", NA),
       vertex.label.cex=1,
       vertex.size = 20, 
       vertex.color = V(graph_t)$color, 
       margin = 0, main = paste0("Splines, Time = ", round(tvec[t], 2)))
  dev.off()
  # fig_list2[[t]] <- image_read(pic2)
}
