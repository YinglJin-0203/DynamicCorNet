df <- read.csv("data/AppData.csv")
library(tidyverse)
AppDataSmall <- df %>% filter(id %in% sample(unique(df$id), size = 30, replace = F))

mat_try <- matrix(0, nrow=3, ncol=3)
mat_try[row(mat_try)==col(mat_try)] <- 1
mat_try[row(mat_try)!=col(mat_try)] <- -0.6
matrixcalc::is.positive.semi.definite(mat_try)

library("colorspace")
hcl_palettes(plot = TRUE)
