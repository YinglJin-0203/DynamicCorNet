# This script is written specifically for data simulation

#### Function ####
# P: number of variables
# G: number of groups
# pGroup: number of variables within each group 
          # please note that single isolated variables are not regarded as groups
# tid: time index
# cor_within: structured as a list, where each element corresponds to a group, 
#             and each element is a 3d array, the third dimension as time
# cor_between: 3D array
# sigma: random error variance
# ft: fixed effect of each variables

GenData_i <- function(P, G, pGroup, tid, cor_between, cor_within, sigma=1,
                      f_pg = function(t){return(rep(0, P))}){
  
  # sanity checks
  if(length(cor_within)!=G | length(pGroup)!=G | dim(cor_between)[1]!=G){
    stop("Error: Dimension of correlation matrix does not match number of groups")}
  
  # container
  df_i <- data.frame(time=tid)
  
  # fixed effect
  ft_vec <- t(sapply(tid, f_pg))
  
  # generate between-group random effect
  a_g <- t(sapply(tid, function(t){MASS::mvrnorm(n=1, mu=rep(0, G), Sigma = cor_between[,,t])}))
  a_g <- a_g[, rep(1:G, times=pGroup)]
  
  # generate within-group random effect
  b_pg <- sapply(tid, 
                 function(t){
                   cor_t <- lapply(cor_within, function(x){x[,,t]})
                   b_pg_t <- lapply(cor_t, 
                                    function(x){
                                      MASS::mvrnorm(n=1, mu=rep(0, nrow(x)), Sigma = x)})
                   b_pg_t <- unlist(b_pg_t)
                   return(b_pg_t)}
  )
  b_pg <- t(b_pg)
  
  # add them up 
  ft_vec[, 1:sum(pGroup)] <- ft_vec[, 1:sum(pGroup)]+a_g+b_pg
  
  # add random error
  ft_vec <- ft_vec+matrix(rnorm(P*length(tid), mean=0, sd=sigma), nrow=length(tid), ncol=P)
  
  # generate isolated variable
  df_i <- cbind(df_i, ft_vec)
  
  return(df_i)
  
}


#### Data generation ####

# scalars 
N <- 100
Tmax <- 20
maxcor <- 0.99 # the maximum correlation strenght 
# Correlation

# Case 1: 
# No between-group correlation
# Within-grouop: 
#   1. positive, increasing 
#   2. positive, decreasing 
#   3. positive, picewise constant
# 3 groups, each has 3 variables

# Case 2: 
  # positive, decreasing correlation between groups 1 and 2
  # Within-grouop: 
  #   1. positive, increasing 
  #   2. positive, decreasing 
  #   3. positive, picewise constant
  # 3 groups, each has 3 variables

# Case 3: 
# negative, decreasing correlation between groups 1 and 2
# Within-grouop: 
#   1. positive, increasing 
#   2. positive, decreasing 
#   3. positive, picewise constant
# 3 groups, each has 3 variables

## group 1
cor_g1 <- lapply(1:Tmax,
                 function(t){
                   cormat <- diag(1, nrow = 3, ncol=3)
                   cormat[row(cormat)!=col(cormat)] <- 0.05*(t-1)
                   return(cormat)
                 })
cor_g1 <- simplify2array(cor_g1)
apply(cor_g1, 3, matrixcalc::is.positive.semi.definite)

## group 2
cor_g2 <- lapply(1:Tmax,
                 function(t){
                   cormat <- diag(1, nrow = 3, ncol = 3)
                   cormat[row(cormat)!=col(cormat)] <- -0.025*(21-t)
                   return(cormat)
                 })
cor_g2 <- simplify2array(cor_g2)
apply(cor_g2, 3, matrixcalc::is.positive.semi.definite)

## group 3
cor_g3 <- lapply(1:Tmax,
                 function(t){
                   cormat <- diag(1, nrow = 3, ncol=3)
                   cormat[row(cormat) != col(cormat)] <- ifelse(t<=10, 0.1, 0.5)
                   return(cormat)
                 })
cor_g3 <- simplify2array(cor_g3)
apply(cor_g1, 3, matrixcalc::is.positive.semi.definite)

## between group correlation
## decreasing positive between group 1 and 2
cor_G <- lapply(1:Tmax,
                function(t){
                  cormat <- diag(1, nrow = 3, ncol = 3)
                  cormat[row(cormat) != col(cormat)] <- 0
                  cormat[1,2] <- cormat[2,1] <- ifelse(0.03*(21-t)>0, 0.03*(21-t), 0.01)
                  return(cormat)
                })

cor_G <- simplify2array(cor_G)
apply(cor_G, 3, matrixcalc::is.positive.semi.definite)

## generation
df_sim <- list()
for(i in 1:N){
  df_sim[[i]] <- GenData_i(P=10, G=3, pGroup=c(3, 3, 3), tid=1:20, cor_between = cor_G, 
                           cor_within = list(cor_g1, cor_g2, cor_g3), sigma = 0.1)
}

df_sim <- bind_rows(df_sim, .id="id")
colnames(df_sim) <- c("id", "time", paste0("X", 1:10))

#### Checks ####

## empirical correlation
## By slice
df_sim %>% 
  group_by(time) %>%
  group_map(.f=~{cor(.x %>% select(starts_with("X")))[1, 1:6]})
## Overall
heatmap(cor(df_sim %>% select(starts_with("X"))))

#### save data ####
write.csv(df_sim, file = here("data/SimCase3.csv"), row.names = F)
rm(list=ls())

