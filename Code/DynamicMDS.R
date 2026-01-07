#' Calculate temporal network layout with dynamic multidimensional scaling
#'
#' @param adj_mat 
#' @param lambda 
#'
#' @returns
#' @export
#'
#' @examples
DynamicMDS <- function(adj_mat, lambda=10){
  
  # dissimilarity matrix
  # package WGCNA doesn't work. I need to write my own TOM function later
  # dis_mat <- lapply(adj_mat, function(x){1-x})
  dis_mat <- adj_mat
  
  # key scalars
  P <- sapply(adj_mat, ncol)
  Tmax <- length(adj_mat)
  
  # initialization
  init_coord <- list()
  for(t in 1:Tmax){
    if(dim(dis_mat[[t]])[1]>2){
      init_coord[[t]] <- smacofSym(dis_mat[[t]], ndim = 2)$conf
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



