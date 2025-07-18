df <- read.csv("data/AppData2.csv")
library(tidyverse)
AppDataSmall2 <- df %>% filter(id %in% sample(unique(df$id), size = 20, replace = F))
write.csv(AppDataSmall2, file = here("data/AppDataSmall2.csv"), row.names = F)
adj_list <- GetAdjMat(data=subset(df, select = -c(id, time_id)), 
                      cor_method = "pearson",
          mds_type = "Splines")
names(adj_list)

graph_list <- graph_dyn_net(adj_list, cor_th = 0.3)
names(graph_list)

SplinesMDS(adj_list, 5, K = 10, P = 20, 
           tvec = as.numeric(unique(df$time)))

mat_try <- matrix(0, nrow=3, ncol=3)
mat_try[row(mat_try)==col(mat_try)] <- 1
mat_try[row(mat_try)!=col(mat_try)] <- -0.6
matrixcalc::is.positive.semi.definite(mat_try)

tvec <- sort(unique(df$time))
which(tvec==7)

rm(list=ls())

library("colorspace")
hcl_palettes(plot = TRUE)
