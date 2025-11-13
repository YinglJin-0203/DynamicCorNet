# This scripts writes the function for dynamic multidimensional scaling
# but computes a layout that is "integrated" over time
# which doesn't change across time


# library(smacof)
# library(mgcv)
# library(tidyverse)
# 
# try <- rnorm(4)
# matrix(try, 2, 2)

#### Estimate by projection ###

##### Stress function ##### 

# Stress function 
# in this case, stress can only be calculated on variables that are measured at all the time points
# also, don't really need to add penalization
# vec: coordinates to estimate. (P^2-P)/2 elements
#   because adjacency matrix must be symmetric
# P is a vector, each element corresponding to the number of nodes at each time point
#   corresonding to 1) number of variables that have observations across all time points; or 
#                   2) all variables  
Integral_stress_DynMDS <- function(vec, diss_list, P) {
  
  # create an adjacency matrix
  IntDis <- matrix(0, nrow = P, ncol = P)
  IntDis[upper.tri(IntDis)] <- vec
  IntDis <- IntDis+t(IntDis)
  
  # Kruskal stress
  Tmax <- length(diss_list)
  stress <- 0
  for (t in 1:Tmax) {
    if(dim(diss_list[[t]])[1]==0){
      stress <- stress+0
    }
    else{
      delta <- diss_list[[t]]
      stress <- stress + sum((IntDis - delta)^2, na.rm = T)
    }
  }
  return(stress)
}

#### Dynamic MDS ####
# data: a dataframe with id, time 
# N: number of subjects
# Tmax: number of measurement points
# tid: measurement grid
# P: number of variables
# cor_method: correlation to calculate. pearson or spearman
# afunc: adjacency function. identify, sigmoid, power
# P is a vector, each element corresponding to the number of nodes at each time point

IntergralDynamicMDS <- function(adj_mat){
  
  # dissimilarity matrix
  # package WGCNA doesn't work. I need to write my own TOM function later
  dis_mat <- lapply(adj_mat, function(x){1-x})
  
  # key scalars
  P <- max(sapply(adj_mat, ncol))
  
  # optimization BFGS
  result <- optim(
    par = runif((P^2-P)/2), # dissimilarity must be between 0 and 1
    fn = Integral_stress_DynMDS,
    diss_list = dis_mat,
    P = P,
    method = "L-BFGS-B",
    control = list(maxit = 1000),
    lower = 0, upper=1
  ) # I will need a convergence warning here
  
  # return both adjacency matrix, dissimilarity matrix and coordinates
  IntDis <- matrix(0, nrow = P, ncol = P)
  IntDis[upper.tri(IntDis)] <- result$par
  IntDis <- IntDis+t(IntDis)
  return(IntDis)
} 



#### test ####
# coords <- DynamicMDS(df, 100, 100, 1:100, 11)
