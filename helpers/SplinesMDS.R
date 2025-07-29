# This script writes functions for splines MDS


library(smacof)
library(splines2)


#### Key scalars ####
# N <- length(unique(df_sim$id)) # sample size
# # time 
# Tmax <- max(df_sim$time)
# tid <- unique(df_sim$time)
# P <- ncol(df_sim %>% dplyr::select(starts_with("X"))) # variable

#### Stress function ####

# xi_vec: flattened coefficient vector with length of P*K*2
# lambda: tuning parameters
# K: spline bases df
# init_coord: P by 2 matirx of initial coordinates 
#             (inital offset so that the points don't all start with zero)

stress_SplMDS <- function(xi_vec, dissim_list, P, K=20, tvec=tvec, lambda=1, init_coord){
  
  # distance on lower dimension
  xi1 = matrix(xi_vec[1: (P*K)], nrow = P)
  xi2 = matrix(xi_vec[(P*K+1): (P*K*2)], nrow = P)
  
  # basis functions matrix, T by K
  # Xmat <- bs(tid, df = K)
  # basic functions and second derivative
  Xmat <- bSpline(tvec, df = K, degree = 3, derivs = 0)
  Xmat2dev <- bSpline(tvec, df = K, degree = 3, derivs = 2)
  
  # coordinates, P by T
  # center around initial coordinates
  c1 <- init_coord[,1] + xi1 %*% t(Xmat)
  c2 <- init_coord[,2] + xi2 %*% t(Xmat)
  
  # pariwise euclidean distance
  tid <- seq_along(tvec)
  coord_list <- lapply(tid, function(t){cbind(c1[, t], c2[, t])})
  dist_list <- lapply(coord_list, dist, diag=T, upper=T)
  dist_list <- lapply(dist_list, as.matrix)
  
  # Kruskal stress
  stress <- mapply(function(diss, dist){
    # high-d dissimilarity
    diss = diss[lower.tri(diss)]
    # low-d distance
    dist = dist[lower.tri(dist)]
    stress_t = sqrt(sum((diss-dist)^2, na.rm = T)/sum(diss^2, na.rm = T))
    return(stress_t)}, 
    dissim_list, dist_list)
  stress <- sum(stress[!is.infinite(stress)], na.rm=T)
  
  # if(any(is.infinite(stress))){
  #   warning(paste("Variables showed perfect similarity at t =", which(is.infinite(stress))))
  # }
  
  # penalization
  penal_c1 <- sapply(1:P, function(x){t(xi1[x, ]) %*% t(Xmat2dev) %*% Xmat2dev %*% xi1[x, ]})
  penal_c2 <- sapply(1:P, function(x){t(xi2[x, ]) %*% t(Xmat2dev) %*% Xmat2dev %*% xi2[x, ]})
  penal <- lambda*(sum(penal_c1)+sum(penal_c2))
  
  loss = stress+penal
  
  return(loss)
}

# stress_SplMDS(dis_try, rnorm(P*K*2), P, K, tid, 10, init_coord)

#### Splines MDS ####

# For splines MDS, layout can be calculated at time slices even when a variables is missing
# adj_mat: adjacency matrix
# lambda: tuning parameter
# K: degree of freedom of spline basis
# P: total number of variables/maximum number of nodes across all time slices
# 

SplinesMDS <- function(adj_mat, lambda, K, P, tvec){
  
  tid <- seq_along(tvec) # time index
  
  # dissimilarity matrix
  # package WGCNA doesn't work. I need to write my own TOM function later
  dis_mat <- lapply(adj_mat, function(x){1-x})
  
  # initial layout at the first time slice
  # variables missing at the first slice: random position
  init_dis <- dis_mat[[1]]
  init_dis[is.na(init_dis)] <- runif(sum(is.na(init_dis)))
  init_coord <- smacofSym(init_dis, ndim=2, init = "random")
  init_coord <- init_coord$conf
  
  # Optimize with regard to 
  final_xi_vec <- optim(
    par = rnorm(P*K*2), 
    fn = stress_SplMDS,
    dissim_list = dis_mat,
    P = P, K = K, tvec = tvec, lambda = lambda, init_coord = init_coord,
    method = "BFGS", control = list(maxit=500))
  
  # calculate coordinates
  xi1 <- matrix(final_xi_vec$par[1: (P*K)], nrow = P)
  xi2 <- matrix(final_xi_vec$par[(P*K+1): (P*K*2)], nrow = P)
  
  # Xmat <- bs(tid, df = 20)
  Xmat <- bSpline(tvec, df = K, degree = 3, derivs = 0)
  c1 <- init_coord[,1] + xi1 %*% t(Xmat)
  c2 <- init_coord[,2] + xi2 %*% t(Xmat)
  
  coords <- lapply(tid, function(x){return(data.frame(c1 = c1[ , x], 
                                                      c2 = c2[ , x]))})
  
  return(coords)
}

# SplinesMDS(adj_try, 10, 30, 32, 1:12)




