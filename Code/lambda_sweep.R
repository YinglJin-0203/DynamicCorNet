
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


