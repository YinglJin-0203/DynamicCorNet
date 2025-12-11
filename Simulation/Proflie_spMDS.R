adj_mat <- adj_mat2
lambda <- 10
K <- 30 
P <- 15

mat1 <- matrix(seq(1:9), 3, 3)
Matrix::tril(mat1, 1)

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
  
  # spline basis 
  Xmat <- bSpline(tvec, df = K, degree = 3, derivs = 0)
  Xmat2dev <- bSpline(tvec, df = K, degree = 3, derivs = 2)
  
  
  profvis::profvis({
   stress_SplMDS(xi_vec = rnorm(P*K*2), 
                 dissim_list = dis_mat,
                 P = P, K = K, tid_vec = tid, lambda = lambda, init_coord = init_coord,
                 Xmat = Xmat, Xmat2dev = Xmat2dev
                 )
  })
  
  # Optimize with regard to 
  profvis::profvis({
  final_xi_vec <- optim(
    par = rnorm(P*K*2), 
    fn = stress_SplMDS,
    dissim_list = dis_mat,
    P = P, K = K, tid_vec = tid, lambda = lambda, init_coord = init_coord,
    method = "BFGS", control = list(maxit=500),
    Xmat = Xmat, Xmat2dev = Xmat2dev)
  })
  
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

# fields::rdist()