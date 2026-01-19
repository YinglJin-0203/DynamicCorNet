# This script writes functions for splines MDS




# For splines MDS, layout can be calculated at time slices even when a variables is missing
# adj_mat: adjacency matrix
# lambda: tuning parameter
# K: degree of freedom of spline basis
# P: total number of variables/maximum number of nodes across all time slices
# 

#' Calculate temporal network layout with Splines Multidimensional Scaling
#'
#' @param adj_mat a list of adjacency matrices over the follow up time
#' @param lambda penalization parameter
#' @param P number of variables
#' @param tvec unique follow up time points
#'
#' @returns
#' @export
#'
#' @examples
SplinesMDS <- function(dis_mat, lambda, P, tvec){
  
  tid <- seq_along(tvec) # time index
  
  # initial layout at the first time slice
  # variables missing at the first slice: random position
  init_dis <- dis_mat[[1]]
  init_dis[is.na(init_dis)] <- runif(sum(is.na(init_dis)))
  init_coord <- mds(init_dis, ndim=2)
  init_coord <- init_coord$conf
  
  # spline basis 
  ## knots
  # if(length(tvec) < 32){
  #   internal_knot <- tvec[-c(1, length(tvec))]
  # } else{
  #   internal_knot <- seq(min(tvec), max(tvec), length.out = 32)
  #   internal_knot <- internal_knot[-c(1, 32)]
  # }
  # design matrix
  # Xmat <- bSpline(tvec, knots = internal_knot, degree = 2, derivs = 0, Boundary.knots = range(tvec))
  # Xmat2dev <- bSpline(tvec, knots = internal_knot, degree = 2, derivs = 2, Boundary.knots = range(tvec))
  # K <- ncol(Xmat)
  tknots <- seq(min(tvec), max(tvec), by = min(diff(tvec)))
  Xmat <- bSpline(tknots, df = 5, degree = 2, derivs = 0)
  Xmat2dev <- bSpline(tknots, df = 5, degree = 2, derivs = 2)
  K <- ncol(Xmat)
  
  # Optimize with regard to 
  final_xi_vec <- optim(
    par = rnorm(P*K*2), 
    fn = stress_SplMDS,
    dissim_list = dis_mat,
    P = P, tid_vec = tid, lambda = lambda, init_coord = init_coord,
    method = "BFGS", control = list(maxit=500),
    Xmat = Xmat, Xmat2dev = Xmat2dev)
  
  # output coefficients and design matrix
  xi1 <- matrix(final_xi_vec$par[1: (P*K)], nrow = P)
  xi2 <- matrix(final_xi_vec$par[(P*K+1): (P*K*2)], nrow = P)
  
  coefs <- list(init_coord = init_coord, 
                xi1 = xi1, 
                xi2 = xi2,
                Xmat = Xmat)
  # Xmat <- bs(tid, df = 20)
  # c1 <- init_coord[,1] + xi1 %*% t(Xmat)
  # c2 <- init_coord[,2] + xi2 %*% t(Xmat)
  # coords <- lapply(tid, function(x){return(data.frame(c1 = c1[ , x], 
  #                                                     c2 = c2[ , x]))})
  
  return(coefs)
}
