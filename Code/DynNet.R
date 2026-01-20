#
#' Function to create layout from a dataset
#'
#' @param df: dataset, must have columns "id" and "time"
#' @param mds_type: dynamic or spliens MDS
#' @param cor_type: what is the similarity/dissimilarity measure to use? 
#'
#' @returns
#' @export
#'
#' @examples
#' 
#' 
TempNetGraph <- function(df, mds_type, cor_type){
  
  tuniq <- sort(unique(df[, "time"]))
  # type of MDS
  if(is.null(MDS)){
    mds_type <- ifelse(length(tuniq)<20, "dynamic", "splines")
  }
  
  
}



