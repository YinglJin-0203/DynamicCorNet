# function to calculate dissimilarity matrix
#' @param data a data frame with a column "time"
#' @param method dissimilarity measure, euclidean distance, spearman or pearson correlation
#' @param cor_method if adjacency measure is correlation, what kind of correlation
#' @param mds_type what multidimensional scaling methods to use

SplDissimMat <- function(data, method = "spearman") {
  method_options <- c("spearman", "pearson", "euclidean")
  mid <- pmatch(method, method_options)
  
  if (is.na(mid)) {
    stop("`method` must be one of: spearman, pearson, euclidean")
  }
  if (!"time" %in% names(data)) {
    stop("`data` must contain a `time` column")
  }
  
  feature_names <- setdiff(names(data), "time")
  if (length(feature_names) == 0L) {
    stop("`data` must contain at least one feature column besides `time`")
  }
  
  # preprocess once
  time_fac <- as.factor(data$time)
  split_idx <- split(seq_len(nrow(data)), time_fac)
  feature_mat <- as.matrix(data[feature_names])
  storage.mode(feature_mat) <- "double"
  feature_mat <- scale(feature_mat, center = TRUE, scale = TRUE)
  
  # dissimilarity computation per time point
  if (mid <= 2L) {
    dis_mat <- lapply(split_idx, function(idx) {
      x <- feature_mat[idx, , drop = FALSE]
      cor_mat <- stats::cor(x, method = method, use = "pairwise.complete.obs")
      1 - abs(cor_mat)
    })
  } else {
    dis_mat <- lapply(split_idx, function(idx) {
      x <- feature_mat[idx, , drop = FALSE]
      as.matrix(stats::dist(t(x)))
    })
  }
  
  return(dis_mat)
}

# test <- SplDissimMat(df, method = "euclidean")
# lapply(test, function(x)(sum(is.na(x))))
# lapply(test, dim)
# View(test[[12]])
