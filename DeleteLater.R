
df <- read.csv("Data/IFEDDemoData.csv")
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
