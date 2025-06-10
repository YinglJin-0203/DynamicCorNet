# Now we have a outcome of interest,
# we want to fix the outcome of interest at the center of the surface
# and let exposures evolve around it

rm(list=ls())
library(smacof)
library(here)
library(tidyverse)
# library(splines)
library(splines2)
library(RColorBrewer)
library(network)
library(ndtv)
library(networkDynamic)
library(igraph)

set.seed(422)

list.files(here("Codes/plot"))
source(here("Codes/Plot/Stress.R"))

# example data
# data(intelligence) # use the example data
# load(here("Data/SimDataLin.RData"))
# load(here("Data/SimBetweenGroup.RData"))
load(here("Data/SimData/Sim1_Y.RData")) # The dataset with a continuous outcome
colnames(df_sim)

#### Key scalars ####
head(df_sim)
N <- length(unique(df_sim$id)) # sample size
Tmax <- max(df_sim$time)
tid <- unique(df_sim$time)

# number of variables to visualize, including confounder, risk factors and outcomes
P <- ncol(df_sim %>% dplyr::select(starts_with("X"), Y)) # exposure variable


#### EWAS example ####
# a toy example use static network plot

##### step 1: all variables included #####
cormat <- df_sim %>% filter(time==1) %>%
  dplyr::select(paste0("X", 1:6), "Y") %>%
  cor()
# dmat <- 1-abs(cormat) # in this scale method, variables that are negatively correlated are also "similar"
                      # and thus will be put closer together
P <- dim(cormat)[1]-1 # number of covariates

## network plot
cormat_th <- cormat
cormat_th[abs(cormat_th) < 0.3] <- 0# thresholding is only for visualization purposes
netg <- graph_from_adjacency_matrix(cormat_th, 
                                     mode = "undirected",
                                     weighted = T,
                                     diag=F)
# plot(netg)

### edge properties
E(netg)$weight <- abs(E(netg)$weight)
E(netg)$width <- E(netg)$weight*5
E(netg)$color <- ifelse(E(netg)$weight >0, "#99CCFF", "#FF6699")

### vertex propertie
# brewer.pal(5, "Accent")
V(netg)$color <- c(rep("#BEAED4", P), "#FDC086")  # Vertex colors
V(netg)$label <- c(paste0("X", 1:P), "Y")                      
V(netg)$shape <- c(rep("circle", P), "square")  # Vertex shapes
V(netg)$size <- c(rep(15, P), 20) 

jpeg(filename=here("Figure/Static/TopExp3_t1.jpeg"))
plot(netg, margin = 0)
dev.off()

##### Step 2: plot regression model results #####

# If we'd like to do variable selection within group
# we could borrow information from MLE regression model
fit1 <- nlme::lme(Y ~ X1+X2+X3+X4+X5+X6+time, random = ~ 1|id, data = df_sim)
pvals1 <- summary(fit1)$tTable[2:7, "p-value"]
# coef1 <- summary(fit1)$tTable[c("X1", "X2", "X3"), "Value"]
  
## network plot
cormat_th <- df_sim %>% dplyr::select(paste0("X", 1:6), "Y") %>% cor()
cormat_th[abs(cormat_th) < 0.3] <- 0# thresholding by correlation magnitude
cormat_th["Y", 1:6][pvals1>0.05] <- 0
cormat_th[1:6, "Y"][pvals1>0.05] <- 0 # thresholding by multiple regression p value

netg <- graph_from_adjacency_matrix(cormat_th, 
                                    mode = "undirected",
                                    weighted = T,
                                    diag=F)
# plot(netg)

### edge properties
E(netg)$weight <- abs(E(netg)$weight)
E(netg)$width <- E(netg)$weight*5
E(netg)$color <- ifelse(E(netg)$weight >0, "#99CCFF", "#FF6699")

### vertex propertie
# brewer.pal(5, "Accent")
P <- 6
V(netg)$color <- c(rep("#BEAED4", P), "#FDC086")  # Vertex colors
V(netg)$label <- c(paste0("X", 1:P), "Y")                      
V(netg)$shape <- c(rep("circle", P), "square")  # Vertex shapes
V(netg)$size <- c(rep(15, P), 20) 

jpeg(filename=here("Figure/Static/TopExp3_pval.jpeg"))
plot(netg, margin = 0)
dev.off()
