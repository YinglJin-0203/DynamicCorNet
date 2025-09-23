
rm(list = ls())
library(tidyverse)
library(ggdendro)
#### small subset for speed concern ####

expand_grid(time = 1:10, id = 1:10) %>% 
  mutate(var1 = rnorm(100), var2 = rnorm(100)) %>% 
  pivot_longer(3:4) %>%
  ggplot(aes(x=time, y=value, fill=name, colour = name, group = interaction(time, name)))+
  geom_boxplot(position = "dodge2")
