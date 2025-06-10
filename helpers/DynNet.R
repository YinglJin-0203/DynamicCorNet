# This script write functions to produces the dynamic network visualization



library(igraph)
# library(network)
# library(ndtv)
# library(networkDynamic)
# library(clusterGeneration)
library(here)
# theme_set(theme_minimal())

#### Network Plot ####


graph_dyn_net <- function(adj_mat, cor_th = 0.3){
  
  netgraph_list <- list()
  Tmax <- length(adj_mat)

  # compute each slice
  for(t in 1:Tmax){
    # thresholding is only for visualization purposes
    # need to modify later as the adjacency function
    adj_mat_t <- adj_mat[[t]]
    adj_mat_t[which(adj_mat_t < cor_th)] <- 0
    # graph
    net_t <- graph_from_adjacency_matrix(adj_mat_t, 
                                         mode = "undirected", weighted = T, diag=F)
    # edge properties
    # let's not color the edges by sigh for now
    E(net_t)$width <- E(net_t)$weight*5
    
    # node properties
    # node color
    # node size
    
    netgraph_list[[t]] <- net_t 
  }
  
  return(netgraph_list)
  
}



