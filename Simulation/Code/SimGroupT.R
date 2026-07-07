
#' Simulate multivariate temporal data with a low/medium density with variable groups for one subject
#' Both within-group and between-group correlation can be specified and can change over time
#' Regular grid, no missing
#' 
#' @param N sample size
#' @param t current time point
#' @param P Total number of variables
#' @param G Total number of groups
#' @param pGroup Number of variables in each group. A vector of length G
#' @param cor_g Group-level random effect correlation at this time point. Should be a G by G matrix
#' @param cor_pg Group-varaible level random effect correlation at this time point, should be a list of G by G elements
#' @param f The over all mean. A scalar or a matrix of T by P 
#' @param id optional id 
#' 
#' @returns
#' @export
#'
#' @examples


GenDataT <- function(N, t, P, G, pGroup,
                     cor_g, cor_pg, sigma=1, f = 0,
                     id){
  
  # sanity checks
  # if(length(cor_within)!=G){stop("Error: Dimension of within-group correlaiton does not match number of groups")}
  # if(length(pGroup)>G){stop("Error: Number of groups or number of.")}
  ## within group correlation
  ## between group correlation
  
  # container
  df_i <- data.frame(id = 1:N, time = t)
  
  # add fixed effect
  df_i[ , paste0("X", 1:P)] <- f
  
  # generate between-group random effect
  ag <- MASS::mvrnorm(n=N, mu=rep(0, G), Sigma = cor_g)
  ag <- ag[, rep(1:G, times=pGroup)]
  df_i[, paste0("X", 1:P)] <- df_i[, paste0("X", 1:P)] + ag
  
  # generate within-group random effect
  bpg <- MASS::mvrnorm(n=N, mu=rep(0, P), Sigma = cor_pg)
  df_i[, paste0("X", 1:P)] <- df_i[, paste0("X", 1:P)] + bpg
  
  # add random error
  epsilon <- matrix(rnorm(P*nrow(df_i), mean=0, sd=sigma), nrow=nrow(df_i), ncol=P)
  df_i[, paste0("X", 1:P)] <- df_i[, paste0("X", 1:P)] + epsilon
  
  # add id
  df_i[ , "id"] <- id
  
  return(df_i)
  
}

# test
# GenDataT(100, 0.1, 15, 3, rep(5, 3), cor_g = cor_between[1:3, 1:3,1], cor_pg = cor_within)
