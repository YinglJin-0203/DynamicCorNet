# Here is the code used for data simulation
set.seed(123)  # for reproducibility

library(tidyverse)

source(here("Code/SimGroupT.R"))

#### Scalars ####
P <- 15 # total number of variables
pGroup  <- rep(5, 3) # variable in each group
G <- 3 # number of groups
N <- 100 # sample size

#### varying scalar: measure density ####
# nt <- 10
nt <- 100
tvec <- seq(0, 1, length.out = nt)


#### simulation ####
##### changing within-group and constant between-group correlation #####

# group level
cor_group <- matrix(0.5, G, G)
diag(cor_group) <- 1
# cor_group

# within group
cor_group_var <- list()

rho_max <- 0.7
for(t in seq_along(tvec)){
  
  # group 1 not changing
  cor_g1 <- matrix(0.3, nrow=pGroup[1], ncol=pGroup[1])
  diag(cor_g1) <- 1
  
  # group 2 increasing
  cor_g2 <- matrix((rho_max/(nt-1))*(t-1), pGroup[2], pGroup[2])
  diag(cor_g2) <- 1
  
  # group 3 decreasing
  cor_g3 <- 0.7-cor_g2
  diag(cor_g3) <- 1
  
  cor_t <- Matrix::bdiag(cor_g1, cor_g2, cor_g3)
  cor_group_var[[t]] <- cor_t
}

# cor_group_var

# generate
df_sim_temp2 <- lapply(1:nt, 
                       function(t){
                         GenDataT(N=100, t = tvec[t], P = 15, G = 3, pGroup = rep(5, 3),
                                  cor_g = cor_group, cor_pg = cor_group_var[[t]],
                                  sigma=1, f=0, id = 1:N)
                       })
df_sim_temp2 <- bind_rows(df_sim_temp2)

# save data
# saveRDS(df_sim_temp2, file = paste0("data/SimData/df_sim_temp2_t", nt, ".rds"))
df_sample <- df_sim_temp2
write.csv(df_sample, file="../SampleData/df_sample.csv", row.names = F)
