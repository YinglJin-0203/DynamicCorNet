
# function to calculate adjacency matrix
#' @param data a data frame with a column "time"
#' @param adj adjacency measure, corelation or distance
#' @param cor_method if adjacency measure is correlation, what kind of correlation
#' @param mds_type what multidimensional scaling methods to use

GetAdjMat <- function(data, adj = "Correlation", cor_method = "pearson", mds_type="Spline"){
  
  # dissimilarity
  if(adj=="Correlation"){
    cor_mat <- data %>% 
      group_by(time) %>%
      group_map(~{cor(.x, method = cor_method, use = "pairwise.complete.obs")})
    adj_mat <- lapply(cor_mat, abs)
  }
  if(adj=="Distance"){
    # dist_mat <- data %>%
    #   group_by(time) %>%
    #   mutate_at(vars(!time), scale, center = T, scale = T) %>% 
    #   group_map(~{fields::rdist(.x, method = cor_method, use = "pairwise.complete.obs")})
    
  }
  
  # Dynamic MDS: remove NA
  if(mds_type=="Dynamic"){
    adj_mat <- lapply(adj_mat, function(m){
      m[upper.tri(m)] <- 0
      keep_id <- rowSums(is.na(m))==0
      mat <- m[keep_id, keep_id]
      mat <- mat+t(mat)
      return(mat)
    })
  }
  
  return(adj_mat)
}

