
rm(list = ls())
library(tidyverse)

#### small subset for speed concern ####

df <- read.csv("data/AppDataSmall.csv")
data_demo <- df %>% select(Participant, Week, HC_A, RL_A, AGD_)

xlist <- future({list(1:5)})
then(xlist, function(x)x[[1]])
library(ggalluvial)
library(tidyverse)

data(majors)
as.data.frame(UCBAdmissions)

class(df[,"BMIZ"])

df %>% rename(time = Week) %>% group_by(time) %>%
  group_map(~mean(.x$BMIZ, na.rm = T))


adj_mat <- GetAdjMat(data=df %>% select(-Participant) %>% rename(time = Week), 
          mds_type = "Splines")
group_list <- 
  lapply(adj_mat,
         function(x){
           dis_mat <- 1-x
           # remove fully unmeasured variable
           var_id <- which(!is.na(diag(dis_mat)))
           dis_mat <- dis_mat[var_id, var_id]
           # remove pairs of variables whose dissimilarity cannot be computed
           # because there is <2 complete pairs
           dis_mat <- dis_mat[complete.cases(dis_mat), complete.cases(dis_mat)]
           dis_mat <- as.dist(dis_mat)
           hclust_fit <- hclust(dis_mat)
           return(cutree(hclust_fit, k = 3))
         })
group_list
group_c <- brewer.pal(3, "Accent")[group_list[[12]]]
names(group_c) <- names(group_list[[12]])

t_uniq <- sort(unique(df$Week))
df_flow <- bind_rows(group_list, .id = "time") %>% 
  mutate(time=t_uniq) %>%
  pivot_longer(-time) %>%
  mutate(value = as.factor(value), time = as.numeric(time)) 


df_flow %>% 
  ggplot(aes(x=time, y=name))+
  geom_tile(aes(fill=value))+
  scale_x_continuous(breaks = t_uniq)+
  labs(x="Time", y = " ", fill = "Group")
ggsave("images/group1.jpeg", height = 5, width = 8)

df_flow$time <- factor(df_flow$time, levels = t_uniq)

df_flow %>%
  ggplot()+
  geom_flow(aes(x=time, stratum = value, fill=name, alluvium=name))+
  geom_stratum(aes(x=time, stratum = value, fill=value, alluvium=name))+
  labs(x="Time", y = " ", fill = "Group")
ggsave("images/group2.jpeg", height = 5, width = 8)


library(ggplot2)
df <- data.frame(
  x = rep(1:3, each = 3),
  y = rep(1:3, 3),
  group = rep(c("A","B","C"), 3)
)

lapply(df, class)

ggplot(df) +
  geom_point(aes(x=x, y=y, fill = group, color = group)) #+   # both fill & color map to group
  #guides(fill = guide_legend("Group"), color = guide_legend("Group"))                                        

unique(paste0("X", c(1:5, 11:15, 6:10)))
               