
# function to calculate adjacency matrix


GetAdjMat <- function(data, sim = "correlation", cor_method = "pearson",
                      afunc="absolute"){
  # split data by time
  # similarity matrix
  cor_mat <- data %>% 
    dplyr::select(-id) %>%
    group_by(time) %>%
    group_map(~{cor(.x, method = cor_method, use = "pairwise.complete.obs")})
  
  sim_mat <- lapply(cor_mat, abs)
  
  # adjacency matrix: should depend on adjacency function
  adj_mat <- sim_mat
  
  return(adj_mat)
}
