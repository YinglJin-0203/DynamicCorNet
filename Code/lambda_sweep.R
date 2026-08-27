#' Search for the appropriate penalization parameter for dMDS
#'
#' @param obs_sim    list of observed similarity matrices. Each matrix correspondings to observed similarity at a time point
#' @param lambda     vector of values of penalization parameter to search
#' @param d          embedding dimension
#' 
#' @returns a data frame with four columns: lambda (the search grid of regularization parameter), stress, movement and objective
#' 
lambda_sweep <- function(obs_sim, lambdas = 10^seq(-3, 3, by = 0.5), d = 2) {
  n_lambda <- length(lambdas)
  results  <- vector("list", n_lambda)
  
  # Warm start: fit largest lambda first, pass solution to next
  lam_sorted <- sort(lambdas, decreasing = TRUE)
  X_init <- NULL
  
  for (i in seq_along(lam_sorted)) {
    lam <- lam_sorted[i]
    fit <- dyn_mds(obs_sim, lambda = lam, d = d, X_init = X_init)
    results[[i]] <- list(
      lambda    = lam,
      stress    = fit$stress,
      movement  = fit$movement,
      objective = fit$objective
    )
    X_init <- fit$embeddings   
  }
  
  df <- do.call(rbind, lapply(results, as.data.frame))
  df[order(df$lambda), ]
}


