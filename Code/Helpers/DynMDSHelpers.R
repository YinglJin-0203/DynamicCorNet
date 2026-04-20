#### Format coordinates ####
# Flatten coordinates into a single vector
# configs is a list of coordinates, each element corresponds to a time point
flatten_configs <- function(configs) {
  do.call(c, lapply(configs, as.vector))
}

# Reshape vector back into list of coordinates
# because optim only works with 1D parameter
# P is a vector, each element corresponding to the number of nodes at each time point
reshape_configs <- function(vec, P, ndim, Tmax) {
  configs <- vector("list", Tmax)
  
  cumP <- c(0, cumsum(P))
  for (t in 1:Tmax) {
    start <- (cumP[t]*ndim) + 1
    end <- cumP[t+1] * ndim
    configs[[t]] <- matrix(vec[start:end], nrow = P[t], ncol = ndim)
  }
  return(configs)
}


#### Stress function ####

# Stress function with temporal penalty
# vec: coordinates at this iteration
# P is a vector, each element corresponding to the number of nodes at each time point
stress_DynMDS <- function(vec, diss_list, P, ndim, Tmax, lambda) {
  configs <- reshape_configs(vec, P, ndim, Tmax)
  stress <- 0
  # Kruskal stress
  for (t in 1:Tmax) {
    if(dim(diss_list[[t]])[1]==0){
      stress <- stress+0
    }
    else{
      dists <- as.matrix(dist(configs[[t]]))
      delta <- diss_list[[t]]
      stress <- stress + sum((dists - delta)^2)
    }
  }
  
  # penalization
  node_list <- lapply(diss_list, colnames)
  config_list <- configs
  for(i in 2:length(node_list)){
    if(is.null(node_list[[i]])){
      node_list[[i]] <- node_list[[i-1]]
      config_list[[i]] <- config_list[[i-1]]
    }
  }
  for (t in 2:length(node_list)) {
    if(dim(diss_list[[t]])[1]>0){
      # identify nodes that exists in both slices
      node_t <- node_list[[t]]
      node_t_1 <- node_list[[t-1]]
      rid_t <- node_t %in% node_t_1
      rid_t_1 <- node_t_1 %in% node_t
      diff <- config_list[[t]][rid_t] - config_list[[t-1]][rid_t_1]
      stress <- stress + lambda * sum(diff^2)
    }
  }
  return(stress)
}