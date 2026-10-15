
#### Set up ####


library(tidyverse)
library(igraph)
library(ggforce)


source("Code/dyn_mds.R")
source("Code/lambda_sweep.R")
source("Code/dMDS_Helpers.R")
source("Code/get_similarity.R")
source("Code/lcurve_corner_dist.R")
source("Code/lcurve_corner_menger.R")


theme_set(theme_minimal())

set.seed(825)
# load data
df <- read.csv("Manuscripts/Data/IFEDDemoData.csv")
df$ID <- as.factor(df$ID)

#### Dynamic MDS ####

df <- df %>% rename(id=ID, time= Week)
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


# plot
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


##### Supplemental graphs #####

# trajectories
colnames(df)
df_type <- data.frame(name = colnames(df)[-c(1:2)],
           procedure = c(rep("Physical", 5), rep("Specimen", 3), rep("Ultrasound", 4)))
df %>% 
  pivot_longer(3:14) %>%
  filter(!is.na(value)) %>%
  left_join(df_type) %>% 
  # filter(procedure == "Physical") %>%
  filter(procedure == "Ultrasound") %>%
  ggplot() +
  geom_line(aes(x=time, y=value, group = id), linewidth = 0.5, na.rm = T, alpha = 0.5) + 
  geom_smooth(aes(x=time, y=value), na.rm = T, method = "loess")+
  facet_wrap(~name, nrow = 1, scale = "free")+
  labs(x="Week", y= " ")


# pairwise correlation
# compute pairwise correlations at each time point
vars <- setdiff(colnames(df), c("time", "id"))

# get all unique pairs
pairs <- combn(vars, 2, simplify = FALSE)
N <- length(unique(df$id))

# compute correlation for each pair at each time point
cor_df <- map_dfr(pairs, function(pair) {
  df %>%
    group_by(time) %>%
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

# ---- 2. Plot: matrix of trajectories ----
# Order factor levels so the row/col layout is consistent
cor_df <- cor_df %>%
  mutate(var1 = factor(var1, levels = vars),
         var2 = factor(var2, levels = vars))

cor_df %>% 
  filter(!is.na(correlation)) %>%
  filter(Npair >= 10) %>%
  ggplot(aes(x = time, y = correlation)) +
  geom_point(aes(x=time, y=correlation, alpha = Npct)) +
  geom_line(aes(x=time, y=correlation)) +
  geom_hline(yintercept = c(-0.4, 0.4), col="red", linewidth = 0.2, linetype = "dashed") +
  facet_grid(rows = vars(var1), cols = vars(var2)) +
  theme_minimal(base_size = 8) +
  theme(
    strip.text = element_text(size = 6),
    axis.text = element_text(size = 5),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "bottom"
  ) +
  labs(x = "Week", y = " ", alpha = "Proportion of complete pairs")

# heatmaps
library(reshape2)

times <- sort(unique(df$time))


cor_long_all <- map_dfr(times, function(t) {
  cor_mat <- df %>%
    filter(time == t) %>%
    ungroup() %>%
    select(-id, -time) %>%
    cor(method = "spearman", use = "pairwise.complete.obs")
  
  
  reshape2::melt(cor_mat) %>%
    setNames(c("Var1", "Var2", "correlation")) %>%
    mutate(time = t)
})

ggplot(cor_long_all, aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 2) +
  scale_fill_gradient2(
    low      = "#E41A1C",
    mid      = "white",
    high     = "#377EB8",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Correlation"
  ) +
  facet_wrap(~ time, ncol = 4) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    axis.title  = element_blank(),
    strip.text  = element_text(size = 8)
  ) +
  coord_fixed() +
  labs(title = "")

