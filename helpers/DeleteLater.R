#### small subset for speed concern ####

df <- read.csv("data/AppData.csv")
colnames(df)
head(df)

rand_id <- sample(unique(df$id), 20)

df <- df %>% # filter(id %in% rand_id) %>%
  select(id, time, ends_with("_A"), ends_with("Z"), ends_with("_N"), starts_with("Tot")) %>%
  select(-RBudBead_A, -LBudBead_A, -WT_A) %>%
  rename(Participant = id, Week = time)
write.csv(df, file = "data/AppData.csv", row.names = F)

df %>% select(!c("time", "id"))

#### Missing ####
df_sub <- df %>%
  select(id, time, ends_with("_A"), ends_with("Z"), ends_with("_N"), starts_with("Tot")) %>%
  select(-RBudBead_A, -LBudBead_A, -WT_A) %>%
  rename(Participant = id, Week = time)

N <- length(unique(df$id))
week <- sort(unique(df$time)) 
  
pct_miss <- df_sub %>% 
  group_by(Week) %>%
  group_map(~{
    nmiss <- sapply(.x, function(x){sum(is.na(x))})
    as_vector(nmiss)/N
  }) 

pct_miss <- bind_rows(pct_miss) %>% select(-Participant)

library(paletteer)
miss_plot <- data.frame(week, pct_miss) %>% 
  pivot_longer(-week) %>%
  mutate(name=factor(name, levels = colnames(df_sub)[3:19])) %>%
  ggplot()+
  geom_tile(aes(x=week, y=name, fill=value))+
    scale_fill_paletteer_c("pals::coolwarm")+
  scale_x_continuous(breaks = week)+
  labs(x="Week", y="", fill = "Proportion of missing")
ggsave("images/IFED_miss.jpeg", plot = miss_plot, width = 6, height = 3)

#### Complete pair #####
P <- ncol(df)-2
pnames <- colnames(df)[3:(P+2)]
miss_pair_pct <- matrix(NA, P, P)

df_t10 <- df %>% filter(Week==t_uniq[10]) %>% select(-Participant, -Week)
N <- nrow(df_t10)

for(i in 1:(P-1)){
  for(j in (i+1):P){
    
    complete_pair <- !is.na(df_t10[, i]) & !is.na(df_t10[, j])
    miss_pair_pct[i, j] <- miss_pair_pct[j,i] <- sum(complete_pair)

  }
}

colnames(miss_pair_pct) <- rownames(miss_pair_pct) <- pnames
data.frame(miss_pair_pct) %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1) %>%
  mutate(value_plot = ifelse(value <= 2, 1, 0)) %>%
  ggplot()+
  geom_tile(aes(x=var1, y=name, fill=value_plot), show.legend = F)+
  scale_fill_viridis_c(na.value = "white")+
  labs(x="", y="", title = "<=2 complete pairs at week 20", fill = " ")+
  theme(axis.text.x = element_text(angle=45))
ggsave("images/IFED_comp_pair.jpeg", width = 6, height = 6)

#### Load data ####
df <- read.csv("data/AppDataSmall.csv")
df <- read.csv("data/AppSleepData.csv")
set.seed(730)
head(df)

#### Tab 2 ####

GetAdjMat(data= df() %>% select(!c("id", "time")) %>% rename(time = time_id), 
          cor_method = "pearson")

df %>% select(!c("id", "time")) %>% rename(time = time_id) %>% group_by(time) %>%
  group_map(~{cor(.x, method = "pearson", use = "pairwise.complete.obs")})



# only use time points where all variables are measured
head(df)
df <- df %>% filter(time %in% c(4, 8, 16, 24, 32))

#### heatmap and hclust  ####
# explore the effect of correlation type
# also the effect of thresholding
adj_list <- GetAdjMat(df %>% select(-id), cor_method = "spearman")
lapply(adj_list, function(mat)sum(is.na(mat)))
graph_list <- graph_dyn_net(adj_list)

# splines MDS coordinates
t_uniq <- sort(unique(df$time))
coord_list <- SplinesMDS(adj_list, lambda = 10, K = 10, P = dim(adj_list[[1]])[1], 
           tvec = t_uniq)


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

# graph
library(magick)

pearson_list <- list()
spearman_list <- list()
for(input_tid in seq_along(t_uniq)){
  
  graph_t <- graph_list[[input_tid]]
  group_t <- group_list[[input_tid]]
  coord_t <- coord_list[[input_tid]]
  
  V(graph_t)$color <- group_t[V(graph_t)$name]
  
  
  # plot
  # picname <- paste0("images/temp,jpeg")
  jpeg(filename = "images/temp.jpeg", height = 600, width = 600)
  plot(graph_t,
       layout = as.matrix(coord_t),
       vertex.frame.color=ifelse(is.na(V(graph_t)$color), "grey", NA),
       vertex.label.cex=1,
       vertex.size = 20, 
       vertex.color = V(graph_t)$color, 
       margin = 0, main = paste0("Week ", t_uniq[input_tid], ", Spearman"))
  dev.off()
  # pearson_list[[input_tid]] <- image_read("images/temp.jpeg")
  spearman_list[[input_tid]] <- image_read("images/temp.jpeg")
}

pearson_animation <- image_animate(image_join(pearson_list), delay = 500)

# Save to file
image_write(pearson_animation, path = "images/pearson_animation.gif")
image_write(spearman_animation, path = "images/spearman_animation_full.gif")
