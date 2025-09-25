
library(igraph)
df <- read.csv("data/data_demo.csv")
source("helpers/AdjacencyMat.R")
source("helpers/IntegralDynMDS.R")
head(df)

# # maybe introduce unmeasured variables?
# df <- df %>% mutate_at(vars(X11:X15), function(x){ifelse(.$time%%2==0, NA, x)})
# 
# # IFED dataset
# df <- read.csv("data/AppData.csv")
# df <- df %>% rename(id = Participant, time = Week)
# df <- df %>% select(id, time, HC_A, RL_A, AGD_FA_A, AGD_ANT_A, BudBead_A, WAZ, HAZ, WHZ, BMIZ)

# iris data
data()

#### simulated data ####
adj_list <- GetAdjMat(df %>% select(-id), mds_type = "Dynamic")
# dis_list <- lapply(adj_list, function(x){1-x})
# Integral_stress_DynMDS(runif((15*15-15)/2), dis_list, 15)
# 
# # "overall" adjacency matrix
# IntDis <- IntergralDynamicMDS(adj_list)
# IntAdj <- 1-IntDis
# IntAdj[1:6, 1:6]

# what if I just use a simple average
# equivalent to stress when observation is complete
AveAdj <- apply(simplify2array(adj_list), c(1, 2), mean, na.rm = T)
dim(AveAdj)
AveDis <- 1-AveAdj
AveAdj[1:6, 1:6]



# coordinates
coords <- mds(AveAdj)$conf

# cluster on coordinates
clust <- hclust(dist(coords))
jpeg("images/coord_clust.jpeg")
plot(clust, main = "On coordinates", xlab = "", sub = "")
dev.off()

clust <- cutree(clust, k=3)


int_net <- graph_from_adjacency_matrix(AveAdj, 
                                     mode = "undirected", weighted = T, diag=F)

jpeg("images/IntNetIFED_hclust_coord.jpeg")
plot(int_net, layout=coords, vertex.color = clust, edge.color = NA, main = "On coordinates")
dev.off()


