#' One full dynamic MDS fit for a given lambda
#'
#' @param obs_sim    list of observed similarity matrices. Each matrix corresponds to observed similarity at a time point
#' @param lambda     temporal smoothness penalty weight
#' @param d          embedding dimension, which is the dimension of the visualization space
#' @param maxit      maximum number of iterations
#' @param tol        convergence tolerance
#' @param X_init     optional. list of T initial embedding
#' 
#' @returns A list of six elements:
#' \item{embeddings}{a list of T embeddings corresponding to position coordinates at each time point}
#' \item{stress}{the final value of stress}
#' \item{movement}{the final value of movement}
#' \item{objective}{the final value of the objective function}
#' \item{lambda}{the value of regularization parameter}
#' \item{iter}{number of iterations needed to converge}

dyn_mds <- function(obs_sim, lambda = 1, d = 2,
                    maxit = 500, tol = 1e-5, X_init = NULL) {
  T  <- length(obs_sim)
  p  <- nrow(obs_sim[[1]])
  
  # similarity to dissimilarity
  Deltas <- lapply(obs_sim, function(C) as.matrix(corr_to_dist(C)))
  
  # Initialize embedding via classical MDS on time-averaged dissimilarity
  if (is.null(X_init)) {
    Delta_mean <- Reduce("+", Deltas) / T
    cmd  <- cmdscale(as.dist(Delta_mean), k = d)
    X    <- lapply(seq_len(T), function(t) cmd + matrix(rnorm(p * d, sd = 0.01), p, d))
    # add change of layout across time to avoid zero gradient (saddle point)
  } else {
    X <- X_init
  }
  
  obj_prev <- Inf
  
  for (iter in seq_len(maxit)) {
    
    X_new <- vector("list", T)
    
    for (t in seq_len(T)) {
      Dt  <- Deltas[[t]]
      Xt  <- X[[t]]
      Dhat <- as.matrix(dist(Xt))
      Dhat[Dhat < 1e-10] <- 1e-10 # preventing division by zero
      
      # SMACOF Guttman transform (majorization step)
      B <- -Dt / Dhat
      diag(B) <- -rowSums(B)
      Xmaj <- (1/p) * (B %*% Xt)
      
      # Temporal coupling: weighted average with neighbours
      neigh <- matrix(0, p, d)
      w <- 0
      if (t > 1) { neigh <- neigh + X[[t-1]]; w <- w + 1 }
      if (t < T) { neigh <- neigh + X[[t+1]]; w <- w + 1 }
      
      # Closed-form update balancing SMACOF gradient and temporal penalty
      X_new[[t]] <- (Xmaj + lambda * neigh) / (1 + lambda * w)
    }
    
    X <- X_new
    
    # Compute objective
    stress   <- sum(mapply(raw_stress, X, Deltas))
    movement <- sum(sapply(2:T, function(t) sum((X[[t]] - X[[t-1]])^2)))
    obj      <- stress + lambda * movement
    
    if (is.finite(obj_prev) && abs(obj_prev - obj) / (abs(obj_prev) + 1e-10) < tol) break
    obj_prev <- obj
  }
  
  stress   <- sum(mapply(raw_stress, X, Deltas))
  movement <- sum(sapply(2:T, function(t) sum((X[[t]] - X[[t-1]])^2)))
  
  list(
    embeddings = X,
    stress     = stress,
    movement   = movement,
    objective  = stress + lambda * movement,
    lambda     = lambda,
    iters      = iter
  )
}
