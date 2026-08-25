#' Title: L-curve corner detection by the maximum perpendicular distance to the line connecting both ends
#'
#' @param sweep_df: dataframe with stress(loss), movement(penalization) and lambda grid
#'
#' @returns 
#' 
#' 
lcurve_corner_dist <- function(sweep_df) {
  s <- log(sweep_df$stress)
  m <- log(sweep_df$movement)
  
  # Normalize both axes to [0,1] before computing distance
  s_n <- (s - min(s)) / (max(s) - min(s))
  m_n <- (m - min(m)) / (max(m) - min(m))
  
  # Line from first to last point
  p1 <- c(s_n[1], m_n[1])
  p2 <- c(s_n[length(s_n)], m_n[length(m_n)])
  d  <- p2 - p1
  
  # Perpendicular distance from each point to the line
  dist_to_line <- abs(d[2] * s_n - d[1] * m_n + (p2[1]*p1[2] - p2[2]*p1[1])) /
    sqrt(sum(d^2))
  
  idx <- which.max(dist_to_line)
  list(lambda_star = sweep_df$lambda[idx], idx = idx,
       dist_to_line = dist_to_line)
}
