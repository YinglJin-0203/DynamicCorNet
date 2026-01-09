#' Calculate temporal network layout with dynamic multidimensional scaling
#'
#' @param dis_mat: dissimilarity matrix
#' @param lambda: penalization parameter
#'
#' @returns
#' @export
#'
#' @examples
DynamicMDS <- function(dis_mat, lambda=5){
  
  # dissimilarity matrix
  # key scalars
  P <- sapply(dis_mat, ncol)
  Tmax <- length(dis_mat)
  
  # initialization
  init_coord <- list()
  for(t in 1:Tmax){
    if(dim(dis_mat[[t]])[1]>2){
      init_coord[[t]] <- mds(dis_mat[[t]], ndim = 2)$conf
    }
    else{
      init_coord[[t]] <- matrix(rnorm(dim(dis_mat[[t]])[1]*2), ncol = 2)
    }
  }
  
  # optimization BFGS
  result <- optim(
      par = flatten_configs(init_coord),
      fn = stress_DynMDS,
      diss_list = dis_mat,
      P = P,
      ndim = 2,
      Tmax = Tmax,
      lambda = lambda,
      method = "BFGS",
      control = list(maxit = 1000)
    ) # I will need a convergence warning here
  
  coords = reshape_configs(result$par, P=P, ndim=2, Tmax = Tmax)
  
  # return both adjacency matrix, dissimilarity matrix and coordinates
  return(coords)
} 

# test
DynamicMDS(test, 6)


