# This script includes stress function used for MDS


#### Stress function ####

# xi_vec: flattened coefficient vector with length of P*K*2
# lambda: tuning parameters
# K: spline bases df
# init_coord: P by 2 matirx of initial coordinates 
#             (inital offset so that the points don't all start with zero)

stress_function <- function(dissim_list, xi_vec, P=10, K=20, tid=tid, lambda=1, init_coord){
  
  # distance on lower dimension
  xi1 = matrix(xi_vec[1: (P*K)], nrow = P)
  xi2 = matrix(xi_vec[(P*K+1): (P*K*2)], nrow = P)
  
  # basis functions matrix, T by K
  # Xmat <- bs(tid, df = K)
  # basic functions and second derivative
  Xmat <- bSpline(tid, df = K, degree = 3, derivs = 0)
  Xmat2dev <- bSpline(tid, df = K, degree = 3, derivs = 2)
  
  # coordinates, P by T
  # center around initial coordinates
  c1 <- init_coord[,1] + xi1 %*% t(Xmat)
  c2 <- init_coord[,2] + xi2 %*% t(Xmat)
  
  # pariwise euclidean distance
  coord_list <- lapply(tid, function(t){cbind(c1[, t], c2[, t])})
  dist_list <- lapply(coord_list, dist, diag=T, upper=T)
  dist_list <- lapply(dist_list, as.matrix)
  
  # Kruskal stress
  stress <- mapply(function(diss, dist){
    # high-d dissimilarity
    diss = diss[lower.tri(diss)]
    # low-d distance
    dist = dist[lower.tri(dist)]
    stress_t = sqrt(sum((diss-dist)^2)/sum(diss^2))
    return(stress_t)}, 
    dissim_list, dist_list)
  stress <- sum(stress)
  
  # penalization
  penal_c1 <- sapply(1:P, function(x){t(xi1[x, ]) %*% t(Xmat2dev) %*% Xmat2dev %*% xi1[x, ]})
  penal_c2 <- sapply(1:P, function(x){t(xi2[x, ]) %*% t(Xmat2dev) %*% Xmat2dev %*% xi2[x, ]})
  penal <- lambda*(sum(penal_c1)+sum(penal_c2))
  
  loss = stress+penal
  
  return(loss)
}
