
#' Calculate an "average" adjacency matrix over a set of discrete time points

#' @param df a data frame with a column "time"
#' @param adj adjacency measure, corelation or distance
#' @param cor_method if adjacency measure is correlation, what kind of correlation


GetIntAdjMat <- function(df, adj = "Correlation", cor_method = "pearson", weight = T){
  adj_list <- GetAdjMat(df, adj = adj, cor_method = cor_method, mds_type = "Splines")
  adj_array <- simplify2array(adj_list)
  # weight 
  if(weight==T){
    tuniq <- sort(unique(df$time))
    wt <- c(diff(tuniq), 0) # what about the last time point?
  }
  else{
    wt <- rep(1, length(tuniq))
  }
  # average
  adj_ave <- apply(adj_array, c(1, 2), function(x){sum(wt*x, na.rm = T)/sum(wt)})
  return(adj_ave)
}
