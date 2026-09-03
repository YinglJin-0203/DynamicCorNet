library(here)
library(tidyverse)
library(gridExtra)
library(smacof)
library(splines2)
library(RColorBrewer)
library(mgcv)
library(igraph)
library(magick)

source("Code/dyn_mds.R")
source("Code/lambda_sweep.R")
source("Code/dMDS_Helpers.R")
source("Code/get_similarity.R")
source("Code/lcurve_corner_dist.R")
source("Code/lcurve_corner_menger.R")


theme_set(theme_minimal())

set.seed(825)

#### descriptives #####

df <- read.csv("SampleData/IFEDDemoData.csv")
df <- df%>% rename(id=ID, time= Week)
df <- df %>% group_by(time) %>% group_modify(~clean_sparse_columns(.x, min_obs = 10))

obs_cors <- df %>%
  group_by(time) %>%
  group_map(~{get_similarity(.x, use = "pairwise.complete.obs", method = "spearman")})
t_uniq <- sort(unique(df$time))

# LOCF
filled_obs_cor <- obs_cors
filled_obs_cor[[1]][is.na(filled_obs_cor[[1]])] <- 1e-5
for (i in 2:length(filled_obs_cor)) {
    na_mask <- is.na(filled_obs_cor[[i]])
    filled_obs_cor[[i]][na_mask] <- filled_obs_cor[[i - 1]][na_mask]
}
filled_obs_cor

# grid search
lambdas <- seq(0, 10, length.out = 100) # lambda search grid
sweep_smooth <- lambda_sweep(filled_obs_cor, lambdas)
maxdist_lam <- lcurve_corner_dist(sweep_smooth)

# layout
dmds_fit <- dyn_mds(obs_sim = filled_obs_cor, lambda = maxdist_lam$lambda_star, d = 2)  

# t_uniq <- c(1, 4, 16, 24)

par(mfrow=c(4, 3), mar = c(0, 0, 2, 0), oma = c(1, 1, 1, 1))
for(i in seq_along(t_uniq)){
  
  layout_i <- dmds_fit$embeddings[[i]]
  cor_i <- obs_cors[[i]]
  df_i <- df[df$time == t_uniq[i], ]
  
  # edges — exclude NAs before filtering by threshold
  edges <- which(abs(cor_i) > 0.5 & upper.tri(cor_i) & !is.na(cor_i), arr.ind = T)
  from <- rownames(cor_i)[edges[, 1]]
  to   <- colnames(cor_i)[edges[, 2]]
  wt   <- cor_i[edges]
  
  # build edge data frame — empty data frame if no edges
  edge_df <- if (length(wt) > 0) {
    data.frame(from = from, to = to, weight = wt)
  } else {
    data.frame(from = character(0), to = character(0), weight = numeric(0))
  }
      
  # initialize
  net_i <- igraph::graph_from_data_frame(
    d        = edge_df,
    directed = F,
    vertices = data.frame(name = rownames(cor_i))
  )
      
  # visual elements of edges
  E(net_i)$width <- abs(E(net_i)$weight) * 8
  E(net_i)$color <- ifelse(E(net_i)$weight > 0, "steelblue", "tomato")
  
  # visual elements of vertices
  miss_var_id <- sapply(df_i[, rownames(cor_i)], function(x) all(is.na(x)))
  V(net_i)$color       <- ifelse(miss_var_id, NA, "lightgrey")
  V(net_i)$frame.color <- "lightgrey"
  
  plot(net_i, layout = layout_i,
       vertex.size        = 20,
       vertex.label.cex   = 1,
       vertex.color       = V(net_i)$color,
       vertex.frame.color = V(net_i)$frame.color,
       edge.curved        = 0.2,
       main = paste0("t = ", t_uniq[[i]]), 
       margin = 0)
  
}

 
dev.off()
