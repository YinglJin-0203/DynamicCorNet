# This script write functions to produces the dynamic network visualization
# using igraph package. This is for edgelist only



library(igraph)
# library(network)
# library(ndtv)
# library(networkDynamic)
# library(clusterGeneration)
library(here)
# theme_set(theme_minimal())

#### Network Plot ####


graph_dyn_net <- function(adj_mat, cor_th = 0.3, mds_type = "Splines"){
  
  netgraph_list <- list()
  Tmax <- length(adj_mat)

  # compute each slice
  for(t in 1:Tmax){
    
    adj_mat_t <- adj_mat[[t]]
    
    # thresholding is only for visualization purposes
    # need to modify later as the adjacency function
    adj_mat_t[which(adj_mat_t < cor_th)] <- 0
    # distinguish the variables with missing values
    miss_node <- colnames(adj_mat_t)[is.na(diag(adj_mat_t))] 
    
    # graph
    net_t <- graph_from_adjacency_matrix(adj_mat_t, 
                                         mode = "undirected", weighted = T, diag=F)
    # edge properties
    # let's not color the edges by sigh for now
    E(net_t)$width <- E(net_t)$weight*5
    
    # node color: unfilled circle for missing variables
    V(net_t)$color <- ifelse(!V(net_t)$name %in% miss_node, rgb(0.2, 0.4, 0.8, alpha=0.4), NA)
    
    # node properties
    # node color
    # node size
    
    netgraph_list[[t]] <- net_t 
  }
  
  return(netgraph_list)
  
}



