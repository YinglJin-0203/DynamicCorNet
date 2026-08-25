# =============================================================================
# Helper functions for dMDS algorithm
# =============================================================================

#' Convert correlation matrix to dissimilarity (1 - |r|, then scaled)
#' C: observed correlation matrix
corr_to_dist <- function(C) as.dist(1 - abs(C))

#' Weighted MDS stress (raw stress, no normalization)
raw_stress <- function(X, Delta) {
  D <- as.matrix(dist(X))
  sum((D - Delta)^2)
}

#' Clean and sanity checks ####
clean_sparse_columns <- function(df, min_obs = 10) {
  # For each column, count non-missing observations
  obs_counts <- sapply(df, function(col) sum(!is.na(col)))
  
  # Identify columns with fewer than min_obs non-missing values
  sparse_cols <- names(obs_counts)[obs_counts < min_obs]
  
  if (length(sparse_cols) > 0) {
    message("Setting these columns entirely to NA (", 
            paste(sparse_cols, collapse = ", "), 
            ") due to having fewer than ", min_obs, " observations.")
    df[sparse_cols] <- NA
  }
  
  return(df)
}