
# df <- read.csv("Data/IFEDDemoData.csv")
df <- df %>% select(-ID, -Age.at.exam) %>% 
  rename(time=Week)

t_vec <- sort(unique(df$time))
diss_mats <- SplDissimMat(df)
coords <- SplinesMDS(diss_mats, 7, 12, t_vec)

int_diss <- IntDissmMat(df)
range(int_diss)
int_adj <- 1-int_diss
int_adj[int_adj <= 1] <- 0
int_net <- graph_from_adjacency_matrix(int_adj, mode = "undirected", weighted = T, diag=F)
plot(int_net)
E(int_net)$width <- 7*E(int_net)$weight
net_t <- graph_from_adjacency_matrix(diss_mats[[1]])
clusts <- HclustCoord(coords, method = "splines", 3)

plot(net_t)
clust_t <- clusts[[1]]
node_col_t <- brewer.pal(3, "Accent")[clust_t]
V(net_t)$color <- node_col_t

plot(net_t)

class(group_list[tvec==1][[1]])

dendro_data(group_list[tvec==1])

