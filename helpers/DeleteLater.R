df <- read.csv("data/AppData.csv")
library(tidyverse)
AppDataSmall <- df %>% filter(id %in% sample(unique(df$id), size = 30, replace = F))

write.csv(AppDataSmall, file = "data/AppDataSmall.csv", row.names = F)


df[, c("id", "time", "TotVol_Ovary")] %>%
  filter(complete.cases(.)) %>%
  ggplot(aes(x=time, y=as.factor(id)))+
  geom_tile()+
  scale_x_continuous(breaks = tid)+
  theme(axis.text.y = element_blank())
  
sort(unique(df$time))
