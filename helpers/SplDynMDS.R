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

source(here("Codes/Plot/Stress.R"))

# example data
# data(intelligence) # use the example data
# load(here("Data/SimDataLin.RData"))
# load(here("Data/SimBetweenGroup.RData"))
# load(here("Data/SimBetweenGroupY.RData")) # The dataset with a continuous outcome

load(here("Data/SimData/Sim2.RData"))
df_sim$time <- as.numeric(df_sim$time)
colnames(df_sim) <- c("id", "time", paste0("X", 1:10))
head(df_sim)

#### Key scalars ####
N <- length(unique(df_sim$id)) # sample size
# time 
Tmax <- max(df_sim$time)
tid <- unique(df_sim$time)
P <- ncol(df_sim %>% dplyr::select(starts_with("X"))) # variable


#### High dimensional dissimilarity ####

# dissimilarity matrix
# coverted from correlation matrix
dissim_list <- lapply(1:Tmax, function(t) {
  cor_t <- df_sim %>% filter(time == t) %>% dplyr::select(starts_with("X")) %>% cor()
  dissim_t <- sim2diss(cor_t, method = "corr")
  return(dissim_t)
})

lapply(dissim_list, dim)

# test function
## initial coordinates: MDS at the specific time
# init_coord <- matrix(0, P, 2)
# when there is an outcome, I would like to pin the outcome variable at the center
# so I will calculate the coordinates of exposures first
init_coord <- smacofSym(dissim_list[[1]], ndim=2, 
                        init = "random")
init_coord <- init_coord$conf
plot(init_coord[,1], init_coord[,2], 
     col = rep(brewer.pal(4, "Accent"), times = c(3, 3, 3, 1)),
     pch = 20)
# init_coord <- rbind(init_coord, apply(init_coord, 2, mean))
# rownames(init_coord) <- c(rownames(init_coord)[1:10], "Y")

stress_function(dissim_list, xi_vec = rnorm(P*20*2), P = 10, K = 20, tid = tid, init_coord = init_coord)


#### Optimize ####
# lambda = 0, 1, 10
system.time({
          final_xi_vec <- optim(
            par = rnorm(P*30*2), 
            fn = stress_function,
            dissim_list = dissim_list,
            P = P, K = 30, tid = tid, lambda = 10, init_coord = init_coord,
            method = "BFGS", control = list(maxit=500))
})
# took 10 minutes when lambda = 0

K = 30
# calculate coordinates
xi1 <- matrix(final_xi_vec$par[1: (P*K)], nrow = P)
xi2 = matrix(final_xi_vec$par[(P*K+1): (P*K*2)], nrow = P)

# Xmat <- bs(tid, df = 20)
Xmat <- bSpline(tid, df = K, degree = 3, derivs = 0)

c1 <- init_coord[,1] + xi1 %*% t(Xmat)
c2 <- init_coord[,2] + xi2 %*% t(Xmat)

coord_sim2 <- lapply(tid, function(x){return(data.frame(c1 = c1[ , x], 
                                                            c2 = c2[ , x]))})
save(coord_sim2, file = here("Data/Coords/coord_sim2.RData"))

#### Preliminary checks ####
slice <- 100
plot(coord_sim2[[slice]][,1], coord_sim2[[slice]][,2], 
     col = rep(brewer.pal(4, "Accent"), times = c(3, 3, 3, 1)),
     pch = rep(c(20, 17), times = c(10, 1)))

#### Network Plot ####

# first, create the edge list at each time point

igraph_edgelist <- list()

for(t in 1:100){
  cormat_t <- df_sim %>% filter(time == t) %>% 
    dplyr::select(starts_with("X")) %>%
    cor()
  # cormat_t <- abs(cormat_t)
  cormat_t[abs(cormat_t)<0.3] <- 0 
  # thresholding is only for visualization purposes
  # it is not used for coordinate calculation
  net_t <- graph_from_adjacency_matrix(cormat_t, 
                                       mode = "undirected",
                                       weighted = T,
                                       diag=F)
  # edge properties
  E(net_t)$color <- ifelse(E(net_t)$weight >0, "#99CCFF", "#FF6699")
  E(net_t)$weight <- abs(E(net_t)$weight)
  E(net_t)$width <- E(net_t)$weight*5
  # add color of vertices based on generation scheme
  # V(net_t)$color <- c(rep("tomato", 3), rep("skyblue", 3), rep("darkgoldenrod", 4))
  
  # E(net_t)$weight <- E(net_t)$weight^2
  # plot(net_t)
  
  igraph_edgelist[[t]] <- net_t 
}

plot(igraph_edgelist[[100]])
cormat_t

# second, make the graph dynamic
DynNet2 <- network.initialize(P, directed = F)
set.vertex.attribute(DynNet2, "vertex.color",
                     rep(brewer.pal(4, "Accent"), times = c(3, 3, 3, 1)))
# set.vertex.attribute(DynNet2, "vertex.sides",
#                      rep(c(3, 4), times = c(10, 1)))


plot(DynNet2, vertex.col=DynNet2 %v% "vertex.color")

# edge properties to evaluate


for(t in 1:100){
  this_elist <- as_data_frame(igraph_edgelist[[t]], "edges")
  # this_color <- V(igraph_edgelist[[t]])$color
  
  # format
  this_elist$from <- factor(this_elist$from,
                            levels = c(paste0("X", 1:10), "Y"),
                            labels = 1:11)
  this_elist$from <- as.numeric(this_elist$from)
  
  this_elist$to <- factor(this_elist$to,
                          levels = c(paste0("X", 1:10), "Y"),
                          labels = 1:11)
  this_elist$to <- as.numeric(this_elist$to)
  
  # add edge to plot
  add.edges.active(DynNet2, tail=this_elist$from, head = this_elist$to,
                   onset=t, terminus = t+1,
                   names.eval = lapply(1:nrow(this_elist), function(x)c("width", "color")),
                   vals.eval= lapply(1:nrow(this_elist), function(x)data.frame(this_elist$width[x], this_elist$color[x]))
  )
  # add vertices layout to plot
  activate.vertex.attribute(DynNet2, "x", coord_sim2[[t]][, 1], onset=t, terminus = t+1)
  activate.vertex.attribute(DynNet2, "y", coord_sim2[[t]][, 2], onset=t, terminus = t+1)
}


# put into animation
compute.animation(DynNet2, 
                  animation.mode = "useAttribute",
                  slice.par=list(start=1, end=100, interval=1, 
                                 aggregate.dur=1, rule='all'))

# Show time evolution through static images at different time points:
render.d3movie(DynNet2, usearrows = F, 
               # vertex visual properties
               displaylabels = T, label.pos=6, 
               vertex.cex = 2,
               vertex.col = DynNet2 %v% "vertex.color",
               # vertex.sides = DynNet2 %v% "vertex.sides",
               # edge visual properties
               edge.lwd = function(slice){2*slice %e% "width"},
               edge.col = function(slice){slice %e% "color"},
               # other
               launchBrowser=T, filename=here("Figure/SplMDS_Sim2.html"),
               render.par = list(show.time=F))

#### Compare coordinate smoothness ####

DynMDScoord_p10 %>% 
  lapply(as.data.frame) %>%
  lapply(function(x){x <- x %>% mutate(var = paste0("X", 1:10))}) %>% 
  bind_rows(.id="time") %>%
  mutate(time=as.numeric(time), 
         var=factor(var, levels = paste0("X", 1:10))) %>% 
  filter(var %in% c("X1", "X4", "X7")) %>%
  ggplot()+
  geom_line(aes(x=time, y = V1, col = "c1"))+
  geom_line(aes(x=time, y = V2, col = "c2"))+
  facet_wrap(~var)+
  labs(x="time", y="", title = "Dynamic MDS, lambda = 10", col = "")
ggsave(filename = here("Figure/DynMDSCoord_p10.jpeg"), width = 9, height=3)






