##### Stress function at one time point #####

.splmds_env <- new.env(parent = baseenv())

.splmds_init_rcpp <- function() {
  if (!is.null(.splmds_env$cpp_ready)) {
    return(invisible(.splmds_env$cpp_ready))
  }
  
  .splmds_env$cpp_ready <- FALSE
  
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  
  code <- '
  #include <Rcpp.h>
  using namespace Rcpp;
  
  // [[Rcpp::export]]
  double splmds_stress_from_coords_cpp(NumericVector c1,
                                       NumericVector c2,
                                       NumericVector diss_vec) {
    const int n = c1.size();
    const int expected = n * (n - 1) / 2;
    if (diss_vec.size() != expected) {
      stop("diss_vec length does not match lower-triangle size");
    }
    
    double num = 0.0;
    double den = 0.0;
    int k = 0;
    for (int i = 1; i < n; ++i) {
      const double x_i = c1[i];
      const double y_i = c2[i];
      for (int j = 0; j < i; ++j) {
        const double dx = x_i - c1[j];
        const double dy = y_i - c2[j];
        const double dist = std::sqrt(dx * dx + dy * dy);
        const double d = diss_vec[k++];
        const double r = d - dist;
        num += r * r;
        den += d * d;
      }
    }
    return std::sqrt(num / den);
  }
  '
  
  ok <- TRUE
  tryCatch({
    Rcpp::cppFunction(code = code, env = .splmds_env)
  }, error = function(e) {
    ok <<- FALSE
  })
  .splmds_env$cpp_ready <- ok
  invisible(ok)
}

.splmds_stress_from_coords <- function(c1, c2, diss_vec, use_rcpp = TRUE) {
  if (use_rcpp) {
    .splmds_init_rcpp()
    if (isTRUE(.splmds_env$cpp_ready) && exists("splmds_stress_from_coords_cpp", envir = .splmds_env, inherits = FALSE)) {
      f_cpp <- get("splmds_stress_from_coords_cpp", envir = .splmds_env, inherits = FALSE)
      return(f_cpp(c1, c2, diss_vec))
    }
  }
  
  dist_t <- as.vector(stats::dist(cbind(c1, c2), method = "euclidean", diag = FALSE, upper = FALSE))
  sqrt(sum((diss_vec - dist_t)^2) / sum(diss_vec^2))
}

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
#' @param lower_idx optional vector of lower-triangle indices for `diss_t`
#' @param use_rcpp logical; if TRUE and Rcpp is available, use C++ stress kernel
#'
#' @returns
#' @export
#'
#' @examples
SplMDS_stress_t <- function(t, xi1, xi2, diss_t, P, 
                            init_coord, lambda, Xmat, Xmat2dev,
                            lower_idx = NULL, use_rcpp = TRUE) {
  
  # calculates coordinates at t (P by 2)
  x_t <- Xmat[t, ]
  c1 <- init_coord[, 1] + as.vector(xi1 %*% x_t)
  c2 <- init_coord[, 2] + as.vector(xi2 %*% x_t)
  
  # Kruskal stress
  if (is.null(lower_idx)) {
    lower_idx <- which(lower.tri(diss_t, diag = FALSE))
  }
  diss_vec <- diss_t[lower_idx]
  stress_t <- .splmds_stress_from_coords(c1 = c1, c2 = c2, diss_vec = diss_vec, use_rcpp = use_rcpp)
  
  # penalization
  x_t_2dev <- Xmat2dev[t, ]
  penal_c1 <- sum((xi1 %*% x_t_2dev)^2)
  penal_c2 <- sum((xi2 %*% x_t_2dev)^2)
  penal_t <- lambda * (penal_c1 + penal_c2)
  
  return(stress_t + penal_t)
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
#' @param diss_vec_list optional precomputed lower-triangle dissimilarities
#' @param use_rcpp logical; if TRUE and Rcpp is available, use C++ stress kernel
#'
#' @returns
#' @export
#'
#' @examples
stress_SplMDS <- function(xi_vec, tid_vec, dissim_list, P, 
                          init_coord, lambda, Xmat, Xmat2dev,
                          diss_vec_list = NULL, use_rcpp = TRUE) {
  K <- ncol(Xmat)
  xi1 <- matrix(xi_vec[1:(P * K)], nrow = P)
  xi2 <- matrix(xi_vec[(P * K + 1):(P * K * 2)], nrow = P)
  
  n_t <- length(tid_vec)
  lower_idx <- which(lower.tri(matrix(FALSE, nrow = P, ncol = P), diag = FALSE))
  
  # Precompute all coordinates at once: P x n_t matrices
  c1_all <- xi1 %*% t(Xmat[tid_vec, , drop = FALSE])
  c2_all <- xi2 %*% t(Xmat[tid_vec, , drop = FALSE])
  c1_all <- sweep(c1_all, 1L, init_coord[, 1], `+`)
  c2_all <- sweep(c2_all, 1L, init_coord[, 2], `+`)
  
  # Precompute penalty terms across all time points
  z1 <- xi1 %*% t(Xmat2dev[tid_vec, , drop = FALSE])
  z2 <- xi2 %*% t(Xmat2dev[tid_vec, , drop = FALSE])
  penal_vec <- lambda * (colSums(z1^2) + colSums(z2^2))
  
  # Stress term per time point
  if (is.null(diss_vec_list)) {
    diss_vec_list <- lapply(dissim_list, function(x) x[lower_idx])
  }
  
  if (use_rcpp) {
    .splmds_init_rcpp()
  }
  
  stress <- numeric(n_t)
  for (i in seq_len(n_t)) {
    diss_vec <- diss_vec_list[[i]]
    stress[i] <- .splmds_stress_from_coords(
      c1 = c1_all[, i],
      c2 = c2_all[, i],
      diss_vec = diss_vec,
      use_rcpp = use_rcpp
    ) + penal_vec[i]
  }
  
  return(sum(stress, na.rm = TRUE))
}
