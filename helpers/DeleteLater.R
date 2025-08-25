
rm(list = ls())

#### small subset for speed concern ####

df <- read.csv("data/AppDataSmall.csv")
colnames(df)

df %>% select(Participant, Week, Testosterone_N) %>% 
  filter(Participant==520102)

df %>% 
  select(Participant, Week, TotVol_Uterus, TotVol_Thy) %>% 
  group_by(Week) %>%
  group_modify(~{data.frame(cor = cor(.x[, c("TotVol_Uterus", "TotVol_Thy")], 
                                      method = "pearson",
                                      use = "pairwise.complete.obs")[1, 2])})  %>%
  ungroup() %>% 
  filter(complete.cases(.)) %>%
  ggplot()+
  geom_point(aes(x=Week, y=cor), na.rm = T)+
  geom_line(aes(x=Week, y=cor), na.rm = T)#+
  labs(title = "Empirical correlation", x = input$time_var, y = " ")+
  scale_x_continuous(breaks = t_uniq)

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
# container
miss_pair_pct <- matrix(NA, P, P)
colnames(miss_pair_pct) <- rownames(miss_pair_pct) <- pnames
View(miss_pair_pct)

t_uniq <- sort(unique(df$Week))
t_uniq[10]
df_t10 <- df %>% filter(Week==t_uniq[10]) %>% select(-Participant, -Week)
N <- nrow(df_t10)

for(i in seq_along(pnames)){
  for(j in seq_along(pnames)){
    
    complete_pair <- (!is.na(df_t10[, pnames[i]])) & (!is.na(df_t10[, pnames[j]]))
    miss_pair_pct[pnames[i], pnames[j]] <- sum(complete_pair)

  }
}

t_uniq[10]

data.frame(miss_pair_pct) %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2") %>% #filter(var1 == "BudBead_A") %>% View()
  mutate_at(c("var1", "var2"), factor, levels = pnames) %>% 
  mutate(pct=value/N, 
         flag = ifelse(value<=2, "8", NA)) %>%
  ggplot()+
  geom_tile(aes(x=var1, y=var2, fill=pct))+
  # geom_point(aes(x=var1, y=var2, shape=flag))+
  scale_fill_paletteer_c("pals::coolwarm", direction = -1)+
  labs(x="", y="", title = "Proportion of complete pairs at week 28", fill = " ")+
  theme(axis.text.x = element_text(angle=45, vjust = 0.6))
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
