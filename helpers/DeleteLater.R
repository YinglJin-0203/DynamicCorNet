#### small subset for speed concern ####

df <- read.csv("data/AppData2.csv")
colnames(df)
head(df)

rand_id <- sample(unique(df$id), 50)

sub_df <- df %>% filter(id %in% rand_id) %>% 
  select(id, time, ends_with("_A"), ends_with("Z"), ends_with("_N"), starts_with("Tot")) %>%
  select(-RBudBead_A, -LBudBead_A, -WT_A)
write.csv(sub_df, file = "data/AppDataSmall.csv", row.names = F)

#### hierarchical clustering ####

# dissimilarity mat
adj_list <- GetAdjMat(sub_df %>% select(-id))

group_list <- lapply(adj_list,
       function(x){
         dis_mat <- 1-x
         var_id <- which(!is.na(diag(dis_mat)))
         dis_mat <- dis_mat[var_id, var_id]
         dis_mat <- as.dist(dis_mat)
         hclust_fit <- hclust(dis_mat)
         return(cutree(hclust_fit, k = 3))
       })

graph_list <- graph_dyn_net(adj_list, cor_th = 0.3)

lapply(1:length(graph_list),
       function(x){
         V(graph_list[[x]])$color <- group_list[[x]][V(graph_list[[x]])]
         return(V(graph))
       })
for(i in seq_along(graph_list)){
V(graph_list[[i]])$color <- group_list[[i]][V(graph_list[[i]])]
print(V(graph_list[[i]])$color)

}

adj_mat <- adj_list[[1]]
dis_mat <- 1-adj_list[[1]]
var_id <- which(!is.na(diag(dis_mat)))
dis_mat <- dis_mat[var_id, var_id]
dis_mat <- as.dist(dis_mat)
##  sanity check
range(dis_mat, na.rm = T)
hclust_fit <- hclust(dis_mat)
plot(hclust_fit)
groups <- cutree(hclust_fit, k = 3)
class(groups)


## color
graph_list <- graph_dyn_net(adj_list, 0.3)
V(graph_list[[1]])$color <- groups[V(graph_list[[1]])]
