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
  ## all variables
  nodes <- unique(unlist(sapply(diss_list, colnames)))
  config_list <- configs
  for(i in 1:length(diss_list)){
    rownames(config_list[[i]]) <- rownames(diss_list[[i]])
  }
  ## penalization
  penalization <- 0
  for(var in nodes){
    coord_i <- lapply(config_list, function(x) {
                  if (var %in% rownames(x)) x[var, ] else NULL
              })
    coord_i <- Filter(Negate(is.null), coord_i)
    coord_i <- do.call(rbind, coord_i)
    penal_i <- sum(apply(coord_i, 2, diff)^2)
    penalization = penalization+penal_i
    
  }
  # loss
  stress <- stress + lambda * penalization
  return(stress)
}
