
##### Set up ####

library(here)
library(tidyverse)
library(gridExtra)
library(smacof)
library(splines2)
library(RColorBrewer)
library(mgcv)
library(igraph)
library(magick)
library(ggforce)

source("Code/dyn_mds.R")
source("Code/lambda_sweep.R")
source("Code/dMDS_Helpers.R")
source("Code/get_similarity.R")
source("Code/lcurve_corner_dist.R")
source("Code/lcurve_corner_menger.R")


theme_set(theme_minimal())

set.seed(825)

#### data cleanr #####
df <- read.csv("SampleData/IFEDDemoData.csv")
df$ID <- as.factor(df$ID)
# df <- df %>% select(-X)
# df <- df %>% select(-id)
write.csv(df, "SampleData/IFEDDemoData.csv", row.names = F)


#### descriptives #####

df <- read.csv("SampleData/IFEDDemoData.csv")
df <- df%>% rename(id=ID, time= Week)
df <- df %>% group_by(time) %>% group_modify(~clean_sparse_columns(.x, min_obs = 10))

obs_cors <- df %>%
  select(-id) %>%
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

# t_uniq <- c(4, 16, 24, 32)
# t_uniq4 <- c(3, 7, 9, 11)
# t_uniq
# i <- 10

par(mfrow=c(4, 3), mar = c(2, 0, 4, 0))
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
       margin = c(0, 0, 0, 0),
       main = paste0("Week ", t_uniq[[i]]))
  # title(paste0("Week ", t_uniq[[i]]), line = 1) 
  
}

 
dev.off()


##### Descriptives #####

# trajectories
df %>% 
  pivot_longer(3:14) %>%
  filter(!is.na(value)) %>%
  ggplot() +
  geom_line(aes(x=time, y=value, group = id), linewidth = 0.5, na.rm = T, alpha = 0.5) + 
  geom_smooth(aes(x=time, y=value), na.rm = T)+
  facet_wrap(~name, scales = "free", strip.position = "left")+
  labs(x="Week", y= " ")
ggsave("")  

##### correlation #####

# assume df has columns: time, id, and variable columns
# compute pairwise correlations at each time point
vars <- setdiff(colnames(df), c("Week", "ID"))

# get all unique pairs
pairs <- combn(vars, 2, simplify = FALSE)
N <- length(unique(df$ID))

# compute correlation for each pair at each time point
cor_df <- map_dfr(pairs, function(pair) {
  df %>%
    group_by(Week) %>%
    summarise(
      correlation = cor(.data[[pair[1]]], .data[[pair[2]]],
                        use = "pairwise.complete.obs", 
                        method = "spearman"),
      group = paste(pair[1], "vs", pair[2]),
      var1 = pair[1],
      var2 = pair[2],
      Npair = sum(complete.cases(.data[[pair[1]]], .data[[pair[2]]])),
      .groups = "drop"
    ) %>%  mutate(Npct = Npair/N)
})

# plot each pair separately onto two pages
n_pairs <- length(unique(cor_df$group))
plots_per_page <- ceiling(n_pairs / 2)  # split evenly across 2 pages

# page 1
p1 <- cor_df %>%
  mutate(group = str_wrap(group, width = 20)) %>%
  filter(!is.na(correlation)) %>%
  filter(Npair >= 10) %>%
  arrange(var1, var2) %>% 
  ggplot() +
  geom_point(aes(x=Week, y=correlation, alpha = Npct)) +
  geom_line(aes(x=Week, y=correlation)) +
  facet_wrap_paginate(~ group, ncol = 5, nrow = 7, page = 1) +  # adjust nrow/ncol as needed
  labs(x = "Week", y = "",  alpha = "Proportion of complete pairs")+
  theme(legend.position = "bottom")

# page 2
p2 <- p1 + facet_wrap_paginate(~ group, ncol =  5, nrow = 7, page = 2)
p2

p1

cor_df %>%
  filter(!is.na(correlation)) %>%
  filter(Npair >= 10) %>%
  arrange(var1, var2) %>% 
  group_by(group) %>%
  ggplot() + 
  geom_point(aes(x=Week, y=correlation, alpha = Npct))+
  geom_line(aes(x=Week, y=correlation))+
  facet_wrap(~group, ncol = 4, nrow = 17)


