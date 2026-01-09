
# function to calculate dissimilarity matrix
#' @param data a data frame with a column "time"
#' @param method dissimilarity measure, euclidean distance, spearman or pearson correlation
#' @param cor_method if adjacency measure is correlation, what kind of correlation
#' @param mds_type what multidimensional scaling methods to use

DynDissimMat <- function(data, method = "spearman"){
  
  method_options <- c("spearman", "pearson", "euclidean")
  mid <- pmatch(method, method_options)
  
  # preprocess
  data <- data %>%
    group_by(time) %>%
    mutate_all(scale, center = T, scale = T)
  
  # correlation
  if(mid<=2){
    cor_mat <- data %>%
      group_map(~{cor(.x %>% select(where(~ sum(!is.na(.)) > 1)),
                      method = method, use = "pairwise.complete.obs")})
    dis_mat <- lapply(cor_mat, function(x){1-abs(x)})
  }
  else if(mid==3){
    dis_mat <- data %>%
      group_map(~{dist(.x %>% select(where(~ sum(!is.na(.)) > 0)) %>% t())})
    dis_mat <- lapply(dis_mat, as.matrix)
  }
  
  # remove NA due to incomplete pair
  dis_mat <- lapply(dis_mat, function(m){
      m[upper.tri(m)] <- 0
      keep_id <- rowSums(is.na(m))==0
      mat <- m[keep_id, keep_id]
      mat <- mat+t(mat)
      return(mat)
    })
  
  return(dis_mat)
}


df <- read.csv("Data/IFEDDemoData.csv")
df <- df %>% rename(time=Week) %>% select(-ID, -Age.at.exam)
test <- DynDissimMat(df, method = "euclidean")
# lapply(test, function(x)(sum(is.na(x))))
# lapply(test, dim)
# View(test[[12]])
