df <- read.csv("data/AppDataSmall2.csv")
library(tidyverse)

View(GetAdjMat)
adj_list<- GetAdjMat(data=df, sim="correlation", cor_method = "pearson", afunc="absolute", mds_type = "Dynamic")
diag(adj_list[[1]])

DynamicMDS(adj_list, 5)
