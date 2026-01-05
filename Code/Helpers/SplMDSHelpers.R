##### Stress function at one time point #####

#' Calculate stress at one time point
#'
#' @param t current time index
#' @param xi1 
#' @param xi2 
#' @param diss_t P by P matrix, disimilarity list at time t
#' @param P number of variables
#' @param init_coord P by 2 matrix indicating the initialized coordinates
#' @param lambda 
#' @param Xmat 
#' @param Xmat2dev 
#'
#' @returns
#' @export
#'
#' @examples
SplMDS_stress_t <- function(t, xi1, xi2, diss_t, P, 
                            init_coord, lambda, Xmat, Xmat2dev){
  
  # calculates coordinates at t (P by 2)
  K <- ncol(Xmat) # splines dimension
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
  # diss_t <- dissim_list[[t]]
  # diss_t <- Matrix::tril(diss_t, 0)
  diss_t <- diss_t[lower.tri(diss_t)]
  stress_t <- sqrt(sum((diss_t-dist_t)^2)/sum(diss_t^2))
  
  # penalization
  penal_c1 <- sum((xi1 %*% Xmat2dev[t, ])^2)
  penal_c2 <- sum((xi2 %*% Xmat2dev[t, ])^2)
  penal_t <- lambda*(penal_c1+penal_c2)
  
  loss = stress_t+penal_t
  
  return(loss)
}

# test
# xi_vec <- rnorm(P*K*2)
# xi1 = matrix(xi_vec[1: (P*K)], nrow = P)
# xi2 = matrix(xi_vec[(P*K+1): (P*K*2)], nrow = P)
# SplMDS_stress_t(9, xi1, xi2,
#                 diss_t = dis_mat[[9]],
#                 init_coord = init_coord, lambda = 7,
#                 Xmat = Xmat, Xmat2dev = Xmat2dev)

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
  xi1 = matrix(xi_vec[1: (P*K)], nrow = P)
  xi2 = matrix(xi_vec[(P*K+1): (P*K*2)], nrow = P)
  
  
  stress <- mapply(SplMDS_stress_t, 
                   tid_vec, dissim_list,
                   MoreArgs = list(
                   xi1=xi1, xi2=xi2, 
                   P = P,
                   init_coord=init_coord,
                   lambda = lambda,
                   Xmat = Xmat, 
                   Xmat2dev = Xmat2dev))
  stress <- sum(stress, na.rm = T)
  
  return(stress)
}
