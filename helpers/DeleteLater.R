df <- read.csv("data/AppDataSmall2.csv")
library(tidyverse)
AppDataSmall2 <- df %>% filter(id %in% sample(unique(df$id), size = 20, replace = F))
write.csv(AppDataSmall2, file = here("data/AppDataSmall2.csv"), row.names = F)

colnames(df)

var1 <- 'HC_A'
df_sum <- df[, c("time", 'HC_A')]
lapply(df_var, class)
as.factor(df_sum$time)


tb_sum <- df_sum %>%  group_by(time) %>% 
  summarise_at(var1, 
               function(x){
                 sum_vec <- t(c(min(x, na.rm = T), mean(x, na.rm=T), median(x, na.rm=T), 
                                max(x, na.rm = T),
                                sd(x, na.rm=T), sum(is.na(x))))
                return(sum_vec)
               })
colnames(tb_sum) <- c("Time", "Min", "Mean", "Median", "Max", "SD", "Nmiss")
tb_sum
  group_map(~{
    .x %>% select('TotVol_Uterus')
  })

 try_sum <- summary(df$TotVol_Uterus)
as.vector(try_sum)
 as.matrix(try_sum)
 
 
 rm(list = ls())
 