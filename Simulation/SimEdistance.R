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
# source(here("Code/AdjacencyMat.R"))
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
nt <- 10
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
df_sim <- lapply(1:nt, 
                       function(t){
                         GenDataT(N=100, t = tvec[t], P = 15, G = 3, pGroup = rep(5, 3),
                                  cor_g = cor_group, cor_pg = cor_group_var[[t]],
                                  sigma=1, f=0, id = 1:N)
                       })
df_sim <- bind_rows(df_sim)

# save data
# saveRDS(df_sim_temp2, file = paste0("data/SimData/df_sim_temp2_t", nt, ".rds"))

#### calculate dissimilarity ####
df <- read.csv("data/IFEDDemoData.csv")
df %>% select(-ID, -Age.at.exam) %>% 
  group_by(Week) %>% 
  group_map(~{Rfast::Dist(t(.x), method = "euclidean")})


mat1 <- matrix(1:6, nrow = 3, ncol=2)
Rfast::Dist(t(mat1))
outer(mat1, t(mat1), function(x, y){sqrt(sum((x-y)^2))})


# calculate euclidean distance as dissimilarity
dist_list <- df_sim %>% 
  select(-id) %>%
  group_by(time) %>% 
  group_map(~{Rfast::Dist(t(.x))})

dim(dist_list[[10]])

dist_list <- lapply(dist_list,
                     function(x){colnames(x) <- paste0("X", 1:15)})
dim(dist_list[[1]])[1]

#### Calculated coordinates ####

##### Dynamic MDS #####

# step 2: calculate layout
time1 <- system.time({
# prof1 <- profvis::profvis({
  coord_list1 <- DynamicMDS(dist_list, lambda = 10)
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

## initialize 
initg <- make_empty_graph(n = P, directed = FALSE)

## plot
tvec[c(1, 6, 10)]

for(t in c(1, 6, 10)){
  
  graph_t <- initg
  coord_t1 <- coord_list1[[t]] # dynamic
  coord_t2 <- coord_list2[[t]] # splines
  
  V(graph_t)$color <- color_vec[rep(1:3, each = 5)] # color by true group
  if(tvec[t] %in% t_miss){
    V(graph_t)$color[1:5] <- NA
    coord_t1 <- rbind(matrix(NA, nrow=5, ncol=2), coord_t1)
    } # for missing variables
  
  # plot dynamic
  pic1 <- paste0("images/SimFigures/DynMDS_t", t, ".jpeg")
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
  pic2 <- paste0("images/SimFigures/SplMDS_", t, ".jpeg")
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

# anime1 <- image_animate(image_join(fig_list1), delay = 200)
# anime2 <- image_animate(image_join(fig_list2), delay = 200)
# 
# # Save to file
# image_write(anime1, path = "Figure/SimTemp/DynMDS.gif")
# image_write(anime2, path = "Figure/SimTemp/SplMDS.gif")





# network figures
