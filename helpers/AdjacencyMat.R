
# function to calculate adjacency matrix


GetAdjMat <- function(data, sim = "correlation", cor_method = "pearson",
                      afunc="absolute"){
  # split data by time
  # similarity matrix
  cor_mat <- data %>% 
    dplyr::select(-id) %>%
    group_by(time) %>%
    group_map(~{cor(.x, method = cor_method, use = "pairwise.complete.obs")})
  
  # remove NA in the correlation matrix
  # these NAs happen because 
  # 1. one or both variables is all missing, or
  cor_mat <- lapply(cor_mat, function(x){
    col_name <- colnames(x)[!is.na(diag(x))] 
    return(x[col_name, col_name])})
  # 2. there are no pairwise complete observations
  cor_mat <- lapply(cor_mat, function(x){x[complete.cases(x), complete.cases(x)]})
  
  sim_mat <- lapply(cor_mat, abs)
  
  # adjacency matrix: should depend on adjacency function
  adj_mat <- sim_mat
  
  return(adj_mat)
}
