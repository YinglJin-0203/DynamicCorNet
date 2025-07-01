
# function to calculate adjacency matrix


GetAdjMat <- function(data, sim = "correlation", cor_method = "pearson",
                      afunc="absolute", mds_type="Spline"){
  # split data by time
  # similarity matrix
  cor_mat <- data %>% 
    #dplyr::select(-id) %>%
    group_by(time) %>%
    group_map(~{cor(.x, method = cor_method, use = "pairwise.complete.obs")})
  
 
  if(mds_type=="Dynamic"){
    # remove NA in the correlation matrix
    # these NAs happen because 
    # 1. one or both variables is all missing, or
    cor_mat <- lapply(cor_mat, function(x){
      col_name <- colnames(x)[!is.na(diag(x))] 
      return(x[col_name, col_name])})
    # 2. there are no pairwise complete observations
    cor_mat <- lapply(cor_mat, function(x){x[complete.cases(x), complete.cases(x)]})
    # similarity matrix
    sim_mat <- lapply(cor_mat, abs)
    # remove if all correlation is either 1 or -1
    # this happens when sample size is two small(n=1 or 2)
    # invalid_sim_mat <- sapply(sim_mat, function(x){all(round(x, 10)==1)})
    for(i in seq_along(sim_mat)){
      if(all(round(sim_mat[[i]], 10)==1)){
        sim_mat[[i]]<- matrix(0,0,0)
      }
    }
  }
  else{# similarity matrix
    sim_mat <- lapply(cor_mat, abs)}  
  
  # adjacency matrix: should depend on adjacency function
  adj_mat <- sim_mat
  
  return(adj_mat)
}

