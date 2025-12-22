
library(tidyverse)

df <- read.csv("data/IFEDDemoData.csv")

df_try <- df %>% 
  select(!ID) %>% 
  rename(time = Week)

adj_try <- GetAdjMat(df_try, adj = "Correlation", cor_method = "pearson", mds_type = "Dynamic")

lapply(adj_try, function(x)sum(is.na(x)))

cor_mat <- df_try%>% 
  group_by(time) %>%
  group_map(~{cor(.x, method = "pearson", use = "pairwise.complete.obs")})
adj_mat <- lapply(cor_mat, abs)
View(adj_mat[[10]])
dim(adj_mat[[10]])
complete.cases(adj_mat[[10]])
mat <- adj_mat[[10]]

mat[upper.tri(mat)] <- 0
id <- rowSums(is.na(mat))==0
mat1 <- mat[id, id]
mat2 <- (mat1+t(mat1))

df_t <- df_try %>% group_by(time) %>%
  mutate_at(vars(!time), scale, center = T, scale = T) %>% 
  filter(time==2)

Rfast::Dist(df_t %>% select(-time) %>% t(), method = "euclidean")
  group_map(~{Rfast::Dist(.x, method = "euclidean")})

adj_try10 <- adj_try[[10]]

adj_try10[complete.case]

adj_try_int <- simplify2array(adj_try)
ave_adj <- apply(adj_try_int, c(1, 2), mean, na.rm = T)
adj_try_int[1:5, 1:5, 1:10]
View(ave_adj)

rename(time="Week", id = "ID", var="Thyroid.volume") %>%
  mutate(id = as.factor(id)) %>%
  arrange(time) %>%
  group_by(time) %>%
  summarize(Nmiss = sum(is.na(var)), 
            Pctmiss = sum(is.na(var))/length(var),
            label = paste0(Nmiss, " (", round(100*Pctmiss, 2), "%) ")) %>%
  filter(Nmiss > 20)

df_pair <- df[, c("Week", "ID", "Thyroid.volume", "Ovary.volume")] %>%
  rename(time="Week", id = "ID")

var1 <- "Thyroid.volume"
var2 <- "Ovary.volume"
N <- 136
df_cor <- df_pair %>% 
  group_by(time) %>%
  summarize(cor = cor(.data[[var1]], .data[[var2]], 
                      use = "pairwise.complete.obs"),
            Npair = sum(complete.cases(.data[[var1]], .data[[var2]]))/N)
# plot
df_cor %>% 
  filter(complete.cases(.)) %>%
  ggplot()+
  geom_point(aes(x=time, y=cor, color = Npair, size = Npair))+
  geom_line(aes(x=time, y=cor))+
  labs(title = "Empirical correlation",  y = " ",
       color = "% of complete pairs", size = "% of complete pairs")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(order = 1), 
         size = guide_legend(order = 1))

df_pair %>%
  pivot_longer(c("Thyroid.volume", "Ovary.volume")) %>%
  filter(complete.cases(.)) %>%
  ggplot(aes(x = time, y=value, colour = name, group=interaction(id, name)))+
  geom_line(alpha = 0.7, linewidth=0.5)+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")
  # scale_x_continuous(breaks = t_brk)+
  labs(title = "Variable trend", x = input$time_var, y = " ")+
  theme(legend.position = "bottom")
