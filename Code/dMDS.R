# -----------------------------------------------------------------------------
# DYNAMIC MDS CORE
#   Minimises: sum_t stress_t  +  lambda * sum_t ||X_t - X_{t-1}||^2_F
#   via alternating majorization (SMACOF-style) with temporal coupling.
# -----------------------------------------------------------------------------

#' Convert correlation matrix to dissimilarity (1 - |r|, then scaled)
#' C: observed correlation matrix
corr_to_dist <- function(C) as.dist(1 - abs(C))

#' Weighted MDS stress (raw stress, no normalization)
raw_stress <- function(X, Delta) {
  D <- as.matrix(dist(X))
  sum((D - Delta)^2)
}

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

# -----------------------------------------------------------------------------
# LAMBDA SWEEP  — compute stress / movement / objective across lambda grid
# -----------------------------------------------------------------------------
#' obs_corrs: list of correlation matrices from observed dataset

lambda_sweep <- function(obs_corrs, lambdas = 10^seq(-3, 3, by = 0.5), d = 2) {
  n_lambda <- length(lambdas)
  results  <- vector("list", n_lambda)
  
  # Warm start: fit largest lambda first, pass solution to next
  lam_sorted <- sort(lambdas, decreasing = TRUE)
  X_init <- NULL
  
  for (i in seq_along(lam_sorted)) {
    lam <- lam_sorted[i]
    fit <- dyn_mds(obs_corrs, lambda = lam, d = d, X_init = X_init)
    results[[i]] <- list(
      lambda    = lam,
      stress    = fit$stress,
      movement  = fit$movement,
      objective = fit$objective
    )
    X_init <- fit$embeddings   # warm start
    # cat(sprintf("lambda = %7.4f | stress = %8.2f | movement = %8.2f | obj = %8.2f\n | iter = %8.2f\n",
    #             lam, fit$stress, fit$movement, fit$objective, fit$iter))
  }
  
  df <- do.call(rbind, lapply(results, as.data.frame))
  df[order(df$lambda), ]
}


