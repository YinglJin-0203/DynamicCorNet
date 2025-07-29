#### small subset for speed concern ####

df <- read.csv("data/AppData2.csv")
colnames(df)
head(df)

rand_id <- sample(unique(df$id), 50)

sub_df <- df %>% filter(id %in% rand_id) %>% 
  select(id, time, ends_with("_A"), ends_with("Z"), ends_with("_N"), starts_with("Tot")) %>%
  select(-RBudBead_A, -LBudBead_A, -WT_A)
write.csv(sub_df, file = "data/AppDataSmall.csv", row.names = F)

#### Load data ####
df <- read.csv("data/AppDataSmall.csv")
df %>% filter(time==1) %>% View()

#### What happened at time point 1?  ####
adj_list <- GetAdjMat(df %>% select(-id))
graph_list <- graph_dyn_net(adj_list)
V(graph_list[[1]])$color
V(graph_list[[1]])

# MDS 
t_uniq <- sort(unique(df$time))
coord_list <- SplinesMDS(adj_list, lambda = 10, K = 10, P = dim(adj_list[[1]])[1], 
           tvec = t_uniq)
coord_list[[1]]

# clustering
group_list <- lapply(adj_list,
         function(x){
           dis_mat <- 1-x
           var_id <- which(!is.na(diag(dis_mat)))
           dis_mat <- dis_mat[var_id, var_id]
           dis_mat <- as.dist(dis_mat)
           hclust_fit <- hclust(dis_mat)
           return(cutree(hclust_fit, k = 3))
         })
# group color
graph_t <- graph_list[[1]]
V(graph_t)

color_t <- group_list[[1]]
color_t

coord_t <- coord_list[[1]]
names(color_t)
V(graph_t)$color <- color_t[V(graph_t)$name]
