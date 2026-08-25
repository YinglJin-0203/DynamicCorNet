#' One full dynamic MDS fit for a given lambda
#'
#' @param obs_corrs  list of T observed correlation matrices
#' @param lambda     temporal smoothness penalty weight
#' @param d          embedding dimension
#' @param maxit      max iterations
#' @param tol        convergence tolerance
#' @param X_init     optional list of T initial embeddings

dyn_mds <- function(obs_corrs, lambda = 1, d = 2,
                    maxit = 500, tol = 1e-5, X_init = NULL) {
  T  <- length(obs_corrs)
  p  <- nrow(obs_corrs[[1]])
  
  # Target dissimilarity matrices
  Deltas <- lapply(obs_corrs, function(C) as.matrix(corr_to_dist(C)))
  
  # Initialise embeddings via classical MDS on time-averaged dissimilarity
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
