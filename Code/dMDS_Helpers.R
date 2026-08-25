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