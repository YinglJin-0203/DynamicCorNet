# This script writes functions for splines MDS


#### Key scalars ####
# N <- length(unique(df_sim$id)) # sample size
# # time 
# Tmax <- max(df_sim$time)
# tid <- unique(df_sim$time)
# P <- ncol(df_sim %>% dplyr::select(starts_with("X"))) # variable

#### Stress function ####

# t: time point
# xi1, xi2: P by K matrix
# dissim_list: time-stratified dissimilarity matrix
# P: number of variables
# K: DF of spline basis
# init_coord: initialization of layout. P by 2 matrix
# lambda: penalization parameter

SplMDS_stress_t <- function(t, xi1, xi2, dissim_list, P, init_coord, lambda, Xmat, Xmat2dev){
  
  # coordinates, P by T
  # center around initial coordinates
  K <- ncol(Xmat)
  c1 <- init_coord[,1] + xi1 %*% Xmat[t, ]
  c2 <- init_coord[,2] + xi2 %*% Xmat[t, ]
  # c1 <- c1[,1]
  # c2 <- c2[,1]
  
  # pariwise euclidean distance
  # dist_t <-as.matrix(dist(cbind(c1,c2), diag = T, upper = T)) # dist() is a o(n^2) operation
  # dist_t <- sqrt(outer(c1, c1, "-")^2 + outer(c2, c2, "-")^2)
  # dist_t <- Rfast::Dist(cbind(c1, c2))
  dist_t <- fields::rdist(cbind(c1, c2), cbind(c1, c2), compact = T)
  dist_t <- dist_t[lower.tri(dist_t)]
  # dist_t <- Matrix::tril(dist_t, 0)
  
  # Kruskal stress
  diss_t <- dissim_list[[t]]
  # diss_t <- Matrix::tril(diss_t, 0)
  diss_t <- diss_t[lower.tri(diss_t)]
  stress_t = sqrt(sum((diss_t-dist_t)^2)/sum(diss_t^2))
  
  # penalization
  penal_c1 <- sum((xi1 %*% Xmat2dev[t, ])^2)
  penal_c2 <- sum((xi2 %*% Xmat2dev[t, ])^2)
  penal_t <- lambda*(penal_c1+penal_c2)
  
  loss = stress_t+penal_t
  
  return(loss)
}

#### Overall stress ####
# xi_vec: 1D vector with length 2*P*K
# t_vec: all unique time points

stress_SplMDS <- function(xi_vec, tid_vec, dissim_list, P, init_coord, lambda, Xmat, Xmat2dev){
  
  K <- ncol(Xmat)
  xi1 = matrix(xi_vec[1: (P*K)], nrow = P)
  xi2 = matrix(xi_vec[(P*K+1): (P*K*2)], nrow = P)
  
  
  stress <- sapply(tid_vec, SplMDS_stress_t, 
                   xi1=xi1, xi2=xi2, dissim_list = dissim_list, 
                   P = P,
                   init_coord=init_coord,
                   lambda = lambda,
                   Xmat = Xmat, 
                   Xmat2dev = Xmat2dev)
  stress <- sum(stress, na.rm = T)
  
  return(stress)
}

#### Splines MDS ####

# For splines MDS, layout can be calculated at time slices even when a variables is missing
# adj_mat: adjacency matrix
# lambda: tuning parameter
# K: degree of freedom of spline basis
# P: total number of variables/maximum number of nodes across all time slices
# 

SplinesMDS <- function(adj_mat, lambda, P, tvec){
  
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
  
  # spline basis 
  internal_t <- tvec[-c(1, length(tvec))]
  internal_knot <- ifelse(length(internal_t) < 30, internal_t, 
                          seq(0, 1, length.out = 32)[-c(1, 32)])
  Xmat <- bSpline(tvec, knots = internal_knot, degree = 2, derivs = 0)
  Xmat2dev <- bSpline(tvec, knots = internal_knot, degree = 2, derivs = 2)
  K <- ncol(Xmat)
  # Xmat <- mSpline(tvec, df = K, degree = 2, derivs = 0)
  # Xmat2dev <- mSpline(tvec, df = K, degree = 2, derivs = 2)
  
  # Optimize with regard to 
  final_xi_vec <- optim(
    par = rnorm(P*K*2), 
    fn = stress_SplMDS,
    dissim_list = dis_mat,
    P = P, tid_vec = tid, lambda = lambda, init_coord = init_coord,
    method = "BFGS", control = list(maxit=500),
    Xmat = Xmat, Xmat2dev = Xmat2dev)
  
  # calculate coordinates
  xi1 <- matrix(final_xi_vec$par[1: (P*K)], nrow = P)
  xi2 <- matrix(final_xi_vec$par[(P*K+1): (P*K*2)], nrow = P)
  
  # Xmat <- bs(tid, df = 20)
  c1 <- init_coord[,1] + xi1 %*% t(Xmat)
  c2 <- init_coord[,2] + xi2 %*% t(Xmat)
  
  coords <- lapply(tid, function(x){return(data.frame(c1 = c1[ , x], 
                                                      c2 = c2[ , x]))})
  
  return(coords)
}
