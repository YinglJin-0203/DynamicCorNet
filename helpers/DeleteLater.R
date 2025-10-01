
rm(list = ls())
library(tidyverse)
library(ggdendro)
#### small subset for speed concern ####

expand_grid(time = 1:10, id = 1:10) %>% 
  mutate(var1 = rnorm(100), var2 = rnorm(100)) %>% 
  pivot_longer(3:4) %>%
  ggplot(aes(x=time, y=value, fill=name, colour = name, group = interaction(time, name)))+
  geom_boxplot(position = "dodge2")

df <- read.csv("data/AppData.csv")
max(df$Testosterone_N, na.rm = T)
df[which.max(df$Testosterone_N), ] # participant 520102 at week 36

df %>%
  ggplot(aes(x=Week, y=Testosterone_N, group = Week))+
  geom_boxplot(outlier.size = 0.5)
ggsave("images/testosterone.jpeg")
