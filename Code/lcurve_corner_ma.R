#' Title: L-curve corner search using analytically curvature along moving-average smoothed curve
#'
#' @param sweep_df: dataframe with stress(loss), movement(penalization) and lambda grid
#'
#' @returns a list of three elements
#' \item{lambda_star}{a scalar. the value of regularization parameter corresponding to L-curve corner}
#' \item{idx}{a integer. the position of lambda_star on the search grid}
#' \item{kappa}{a vector. the analytical curvature at each of the smoothed L-curve}
#' 
#' @export
#'
#' @examples
lcurve_corner_ma <- function(sweep_df) {
  s <- log(sweep_df$stress)
  m <- log(sweep_df$movement)
  lam <- log(sweep_df$lambda)
  
  # Smooth before differentiating (3-point moving average)
  smooth3 <- function(x) {
    n <- length(x)
    c(x[1], (x[1:(n-2)] + x[2:(n-1)] + x[3:n]) / 3, x[n])
  }
  s <- smooth3(s); m <- smooth3(m)
  
  # Numerical derivatives w.r.t. log(lambda)
  ds <- diff(s) / diff(lam); ds <- c(ds[1], ds)
  dm <- diff(m) / diff(lam); dm <- c(dm[1], dm)
  d2s <- diff(ds) / diff(lam); d2s <- c(d2s[1], d2s)
  d2m <- diff(dm) / diff(lam); d2m <- c(d2m[1], d2m)
  
  # Curvature
  kappa <- abs(d2s * dm - ds * d2m) / (ds^2 + dm^2)^1.5
  kappa[!is.finite(kappa)] <- 0
  
  idx <- which.max(kappa)
  list(
    lambda_star = sweep_df$lambda[idx],
    idx         = idx,
    kappa       = kappa
  )
}