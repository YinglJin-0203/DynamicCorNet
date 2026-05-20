
#' Calculate an "average" adjacency matrix over a set of discrete time points

#' @param df a data frame with a column "time"
#' @param weight does integration needs to be weighted by time
#' @param method if adjacency measure is correlation, what kind of correlation


IntDissmMat <- function(df, method = "spearman", weight = T){
  diss_list <- SplDissimMat(df, method = method)
  diss_array <- simplify2array(diss_list)
  t_uniq <- sort(unique(df$time))
  # weight 
  if(weight==T){
    wt <- diff(t_uniq)
    wt <- c(wt, 0) 
    wt <- wt/sum(wt)
  }
  else{
    wt <- rep(1, length(t_uniq))
  }
  # average
  AveDiss <- apply(diss_array, c(1, 2), function(x){sum(wt*x, na.rm = T)})
  return(AveDiss)
}

# WTdiss <- IntDissmMat(df, weight = T)
# diss <-  IntDissmMat(df, weight = F)
