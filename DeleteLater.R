
df <- read.csv("Data/IFEDDemoData.csv")
df <- df %>% select(-ID, -Length, -Weight.for.age, -Height.for.age, 
                        -Weight.for.height, -Length, -Bud.bead.diameter,
                        -Age.at.exam) %>% 
  rename(time=Week)

diss_mats <- DynDissimMat(df, method = "euclidean")
diss_t <- diss_mats[[1]]
diss_t <- (max(diss_t)-diss_t)/(max(diss_t)-min(diss_t))
diss_t
View(1-diss_t)
min(diss_t[diss_t > 2.220446e-15])

adj_t <- diss_t
adj_t[adj_t > (1-0.89)] <- 0
diag(adj_t) <- 0


net_t <- graph_from_adjacency_matrix(adj_t, weighted = T, mode = "undirected")
E(net_t)$label <- E(net_t)$weight

wt_range <- range(E(net_t)$weight)
E(net_t)$width = (wt_range[2]-E(net_t)$weight)/(wt_range[2]=wt_range[1])
E(net_t)$width <- 10*E(net_t)$width

plot(net_t)


V(net_t)$name <- nodes_t
V(net_t)$color <- ifelse(V(net_t)$name %in% vars_t, rgb(0.2, 0.4, 0.8, alpha=0.4), NA)
adj_mat[adj_mat]
nodes <- colnames(diss_t)

net_t <- make_empty_graph(n = 8, directed = F)
V(net_t)$name <- nodes

for(i in 2:nrow(diss_t)){
  for(j in 1:(i-1)){
    if(diss_t[i,j] < 10){
      net_t <- add_edges(net_t, c(i, j), weight = diss_t[i,j])
    }
  }
}

E(net_t)
plot(net_t, edge.width = 10/E(net_t)$weight)
nodes
V(net_t)$name

t(cbind(edges$from, edges$to))

net_t <- graph_from_data_frame(edges[, 1:2], directed = F)
E(net_t)$weight <- edges$weight
plot(net_t, edge_width = E(net_t)$weight,
     edge_label = E(net_t)$weight)
net_t %>% 
  add_edges()
matrix(edges[,1:2], ncol=1, byrow=T)

which(lower.tri(diss_t), arr.ind = T)
diss_t[lower.tri(diss_t)]
E(net_t)
View(diss_t)
expand.grid(nodes[-1], nodes[-length(nodes)])

1/diss_t
net_t <- make_empty_graph(8)
nodes <- colnames(diss_t)
V(net_t)$name <- nodes
diss_t[lower.tri(diss_t)]


net_t <- net_t %>% add_edges(match(c("Ovary.volume", "Uterus.volume"), V(net_t)$name), weight = 1)
plot(net_t,
     edge_width = E(net_t)$weight,
     edge_label = E(net_t)$weight)

for(i in 2:nrow(diss_t)){
  for(j in 1:(i-1)){
    if(diss_t[i,j]  < 20){
    net_t %>% add_edges(match(c(nodes[i], nodes[j]), V(net_t)$name))
      }
  }
}

g <- make_empty_graph(directed = F) %>%
  add_vertices(length(nodes), name = nodes)

plot(g)
edge_df <- data.frame(diss_t) %>%
  rownames_to_column("from") %>%
  pivot_longer(!from, names_to = "to", values_to = "weight") %>% 
  filter(from!=to & weight <20)
  distinct(.) %>% 
  filter(from != to)
edge_df$weight[edge_df$weight>25] <- NA

net_t <- graph_from_data_frame(edge_df, directed = F)
plot(net_t)


as.vector(edge_df[ ,1:2])
unlist(edge_df[, 1:2])
add_edges(net_t, edge_df)

tvec <- sort(unique(df$Age.at.exam))

graph_t <- make_empty_graph(n = 5, directed = FALSE)
V(graph_t)$name <- c("A", "B", "C", "D", "E")
graph_t %>% add_edges(match(c("A", "B"), V(graph_t)$name))
plot(graph_t)

E(graph_t) <- edge_attr(graph_t)

min(tvec)
max(tvec)
min(diff(tvec))
df <- df %>% rename(time=Week) %>%
  group_by(time) %>%
  mutate_all(scale, center = T, scale = T) 
  
  group_map(~{cor(.x %>% select(where(~ sum(!is.na(.)) > 2)),
                  method = method, use = "pairwise.complete.obs")})
df <- df %>% select(-ID, -Length, -Height.for.age, -Weight.for.age, -Weight.for.height) 
dist_list <- df %>% group_by(Week) %>% 
  mutate_all(scale, center = T, scale = T) %>%
  group_map(~{dist(t(.x))})

View(dist_list[[3]])

test <- matrix(rnorm(20), 10, 2)
sqrt(sum((test[,1]-test[,2])^2, na.rm=T)*10/7)
test[c(2, 5, 7), 2] <- NA
Rfast::Dist(t(test))
dist(t(test))

df <- df %>% select(-Week)
df <- df %>% mutate_all(scale, center = T, scale = T)
dist_mat <- Rfast::Dist(t(df))

adj_mat <- Geselect()adj_mat <- GetAdjMat(data= df, adj = "Correlation", cor_method = "spearman", mds_type = "Splines")
adj_mat[[15]]

t_uniq <- sort(unique(df[, "time"]))
# scaling>
t_uniq <- t_uniq/max(t_uniq)
system.time({
coord_list <- SplinesMDS(adj_mat, lambda = 7, P = dim(adj_mat[[1]])[1], tvec = t_uniq)
})

1355/3600
