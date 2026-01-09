
# function to calculate dissimilarity matrix
#' @param data a data frame with a column "time"
#' @param method dissimilarity measure, euclidean distance, spearman or pearson correlation
#' @param cor_method if adjacency measure is correlation, what kind of correlation
#' @param mds_type what multidimensional scaling methods to use

SplDissimMat <- function(data, method = "spearman"){
  
  method_options <- c("spearman", "pearson", "euclidean")
  mid <- pmatch(method, method_options)
  
  # preprocess
  data <- data %>%
    group_by(time) %>%
    mutate_all(scale, center = T, scale = T)
  
  # correlation
  if(mid<=2){
    cor_mat <- data %>%
      group_map(~{cor(.x, method = method, use = "pairwise.complete.obs")})
    dis_mat <- lapply(cor_mat, function(x){1-abs(x)})
  }
  else if(mid==3){
    dis_mat <- data %>%
      group_map(~{dist(.x %>% t())})
    dis_mat <- lapply(dis_mat, as.matrix)
  }
  

  return(dis_mat)
}

# test <- SplDissimMat(df, method = "euclidean")
# lapply(test, function(x)(sum(is.na(x))))
# lapply(test, dim)
# View(test[[12]])
