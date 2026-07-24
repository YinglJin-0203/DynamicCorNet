#' Calculate a dissimilarity matrix from a multivariate dataset
#'
#' @param data A data frame or matrix of numeric variables (rows = observations,
#'   columns = variables)
#' @param method Character string: one of "pearson", "spearman", or "euclidean"
#' @param use Method for handling missing values when computing correlations,
#'   passed to cor() (default "pairwise.complete.obs")
#'
#' @return A dissimilarity matrix (class "dist" for euclidean, matrix for correlation methods)
get_similarity <- function(data, method = c("pearson", "spearman", "euclidean"),
                               use = "pairwise.complete.obs") {
  
  method <- match.arg(method)
  
  # Ensure data is a matrix of numeric values
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    stop("All variables in 'data' must be numeric.")
  }
  
  # correlation
  if (method %in% c("pearson", "spearman")) {
    
    # Compute correlation matrix between variables (columns)
    cor_mat <- cor(data, method = method, use = use)
    
    return(cor_mat)
    
  } else if (method == "euclidean") {
    
    # Center and scale each variable (column)
    data_scaled <- scale(data, center = TRUE, scale = TRUE)
    
    # Pairwise Euclidean distance between observations (rows)
    dist_mat <- as.matrix(dist(t(data_scaled), method = "euclidean"))
    
    # Min-max normalize the distances to [0, 1]
    # (excluding the diagonal, which is always 0, from affecting the min)
    off_diag <- dist_mat[upper.tri(dist_mat)]
    min_d <- min(off_diag, na.rm = T)
    max_d <- max(off_diag, na.rm = T)
    
    if (max_d == min_d) {
      warning("All pairwise distances are identical; returning zero matrix.")
      dist_norm <- matrix(0, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
    } else {
      dist_norm <- (max_d - dist_mat) / (max_d - min_d)
      diag(dist_norm) <- 1
    }
    
    dimnames(dist_norm) <- dimnames(dist_mat)
    
    return(dist_norm)
  }
}
