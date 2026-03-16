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
#' @param use_rcpp logical; if TRUE and Rcpp is available, use C++ stress kernel
#'
#' @returns
#' @export
#'
#' @examples
SplinesMDS <- function(dis_mat, lambda, P, tvec, use_rcpp = TRUE) {
  n_time <- length(tvec)
  tid <- seq_len(n_time) # time index
  
  if (length(dis_mat) != n_time) {
    stop("`dis_mat` and `tvec` must have the same length")
  }
  
  # initial layout at the first time slice
  # variables missing at the first slice: random position
  init_dis <- dis_mat[[1]]
  if (anyNA(init_dis)) {
    init_dis[is.na(init_dis)] <- runif(sum(is.na(init_dis)))
  }
  init_coord <- mds(init_dis, ndim = 2)$conf
  
  # spline basis at observed time points
  t_input <- as.numeric(tvec)
  spline_df <- min(5L, n_time)
  Xmat <- bSpline(t_input, df = spline_df, degree = 2, derivs = 0)
  Xmat2dev <- bSpline(t_input, df = spline_df, degree = 2, derivs = 2)
  K <- ncol(Xmat)
  
  # Precompute lower-triangle dissimilarities once for optimization loop
  lower_idx <- which(lower.tri(matrix(FALSE, nrow = P, ncol = P), diag = FALSE))
  diss_vec_list <- lapply(dis_mat, function(x) x[lower_idx])

  # Optimize with regard to spline coefficients
  init_par <- stats::rnorm(P * K * 2)
  final_xi_vec <- stats::optim(
    par = init_par,
    fn = stress_SplMDS,
    dissim_list = dis_mat,
    P = P,
    tid_vec = tid,
    lambda = lambda,
    init_coord = init_coord,
    method = "BFGS",
    control = list(maxit = 500),
    Xmat = Xmat,
    Xmat2dev = Xmat2dev,
    diss_vec_list = diss_vec_list,
    use_rcpp = use_rcpp
  )
  
  # output coefficients and design matrix
  xi1 <- matrix(final_xi_vec$par[1:(P * K)], nrow = P)
  xi2 <- matrix(final_xi_vec$par[(P * K + 1):(P * K * 2)], nrow = P)
  rownames(xi1) <- rownames(xi2) <- rownames(dis_mat[[1]])
  
  coefs <- list(
    init_coord = init_coord,
    xi1 = xi1,
    xi2 = xi2,
    Xmat = Xmat
  )
  
  return(coefs)
}
