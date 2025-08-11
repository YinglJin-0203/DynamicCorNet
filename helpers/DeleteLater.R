#### small subset for speed concern ####

df <- read.csv("data/AppData2.csv")
colnames(df)
head(df)

rand_id <- sample(unique(df$id), 20)

sub_df <- df %>% filter(id %in% rand_id) %>% 
  select(id, time, ends_with("_A"), ends_with("Z"), ends_with("_N"), starts_with("Tot")) %>%
  select(-RBudBead_A, -LBudBead_A, -WT_A) %>%
  rename(Participant = id, Week = time)
write.csv(sub_df, file = "data/AppDataSmall.csv", row.names = F)

df %>% select(!c("time", "id"))

#### Load data ####
df <- read.csv("data/AppDataSmall.csv")
set.seed(730)

#### Tab 2 ####
tbrk <- unique(df$time)
df_sum <- df %>% select(id, time, FSH_N) %>%
  rename(var="FSH_N") %>% 
  group_by(time) %>%
  mutate(med=median(var, na.rm = T)) %>% 
  mutate(miss = ifelse(is.na(var), "Missing", "Present"),
         Nmiss = sum(is.na(var)), 
         Pctmiss = sum(is.na(var))/length(var)) 

df_sum %>% select(time, miss, Nmiss, Pctmiss) %>% 
  distinct(.) %>%
  mutate(lab = paste0(Nmiss, " (", round(100*Pctmiss, 2), "%)")) %>%
  ggplot()+
  geom_col(aes(x=time, y = 50))+
  geom_col(aes(x=time, y=Nmiss), col = "red", fill="red")

df_sum %>% ggplot()+
  geom_line(aes(x=time, y=med), alpha=0.5, linetype = "dashed")+
  geom_boxplot(aes(x=time, y=var, group=time), outlier.size = 0.5, fill = "black", alpha=0.5)+
  geom_jitter(aes(x=time, y=var, group=time), size = 0.5)+
  scale_x_continuous(breaks = df_xlab$time, name = "time",
                     sec.axis = sec_axis(~., name = "Missing", breaks = df_xlab$time, label =df_xlab$Nmiss))+
  theme(axis.text.x.top = element_text(angle=45))
  
cat("aaaaa", "bbbbb", "ccccc", sep = "\n")

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
