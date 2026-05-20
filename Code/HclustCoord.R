
#' Fit hierarchical clustering models for the output of Dynamic or Spliens MDS
#'
#' @param layout 
#' @param method 
#'
#' @returns
#' @export
#'
#' @examples
HclustCoord <- function(layout, method = "splines"){
  
  if(method == "splines"){
    # layout parameters
    xi1 <- layout$xi1
    xi2 <- layout$xi2
    init_coords <- layout$init_coord
    Xmat <- layout$Xmat
    
    coord_list <- lapply(1:nrow(Xmat), 
                         function(x){
                           cbind(
                             init_coords[, 1] + xi1 %*% Xmat[x,],
                             init_coords[, 2] + xi2 %*% Xmat[x,]
                           )
                         })
  } else {
    coord_list <- layout
  }
  
  clust_list <- lapply(coord_list, function(c){
    hclust_t <- hclust(dist(c))
    return(hclust_t)
  })
  
  return(clust_list)
  }


# df <- read.csv("Data/IFEDDemoData.csv")
# df <- df %>% select(-ID, -Age.at.exam) %>%
#   rename(time=Week)
# # diss_list <- DynDissimMat(df)
# diss_list <- SplDissimMat(df)
# # try_lo <- DynamicMDS(diss_list)
# try_lo <- SplinesMDS(diss_list, 7, 12, sort(unique(df$time)))
# HclustCoord(try_lo, method = "splines")
