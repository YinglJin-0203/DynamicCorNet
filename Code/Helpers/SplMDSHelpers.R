##### Stress function at one time point #####

#' Calculate stress at one time point
#'
#' @param t current time index
#' @param xi1 
#' @param xi2 
#' @param diss_t P by P matrix, dissimilarity matrix at time t
#' @param P number of variables
#' @param init_coord P by 2 matrix indicating the initialized coordinates
#' @param lambda 
#' @param Xmat 
#' @param Xmat2dev 
#' @param lower_idx optional vector of lower-triangle indices for `diss_t`
#'
#' @returns Stress loss at a certain time points

SplMDS_stress_t <- function(t, xi1, xi2, diss_t, P, 
                            init_coord, lambda, Xmat, Xmat2dev,
                            lower_idx = NULL){
  
  x_t <- Xmat[t, ]
  c1 <- init_coord[,1] + as.vector(xi1 %*% x_t)
  c2 <- init_coord[,2] + as.vector(xi2 %*% x_t)
  
  # pairwise euclidean distance (lower triangle only)
  coord_t <- cbind(c1, c2)
  dist_t <- as.vector(stats::dist(coord_t, method = "euclidean", diag = FALSE, upper = FALSE))
  
  # Kruskal stress
  if (is.null(lower_idx)) {
    lower_idx <- which(lower.tri(diss_t, diag = FALSE))
  }
  diss_vec <- diss_t[lower_idx]
  stress_t <- sqrt(sum((diss_vec - dist_t)^2) / sum(diss_vec^2))
  
  # penalization
  x_t_2dev <- Xmat2dev[t, ]
  penal_c1 <- sum((xi1 %*% x_t_2dev)^2)
  penal_c2 <- sum((xi2 %*% x_t_2dev)^2)
  penal_t <- lambda * (penal_c1 + penal_c2)
  
  loss <- stress_t + penal_t
  
  return(loss)
}


##### Over all stress #####

#' Calculate stress over all time points
#'
#' @param xi_vec 
#' @param tid_vec a vector of time index, not the original time
#' @param dissim_list 
#' @param P 
#' @param init_coord 
#' @param lambda 
#' @param Xmat 
#' @param Xmat2dev 
#'
#' @returns
#' @export
#'
#' @examples
stress_SplMDS <- function(xi_vec, tid_vec, dissim_list, P, 
                          init_coord, lambda, Xmat, Xmat2dev){
  
  K <- ncol(Xmat)
  xi1 <- matrix(xi_vec[1:(P * K)], nrow = P)
  xi2 <- matrix(xi_vec[(P * K + 1):(P * K * 2)], nrow = P)
  
  lower_idx <- which(lower.tri(matrix(FALSE, nrow = P, ncol = P), diag = FALSE))
  n_t <- length(tid_vec)
  stress <- numeric(n_t)
  
  for (i in seq_len(n_t)) {
    stress[i] <- SplMDS_stress_t(
      t = tid_vec[i],
      xi1 = xi1,
      xi2 = xi2,
      diss_t = dissim_list[[i]],
      P = P,
      init_coord = init_coord,
      lambda = lambda,
      Xmat = Xmat,
      Xmat2dev = Xmat2dev,
      lower_idx = lower_idx
    )
  }
  
  return(sum(stress, na.rm = TRUE))
}
