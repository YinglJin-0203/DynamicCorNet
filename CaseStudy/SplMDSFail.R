library(tidyverse)


df <- read.csv("data/IFEDDemoData.csv")
head(df)

df <- df %>% rename(id = ID, time = Week)

# dissimilarity matrix
dis_mat <- SplDissimMat(df %>% select(-id), method = "spearman")
dis_mat
lapply(dis_mat, dim)

# Splines MDS
t_uniq <- sort(unique(df$time))
mds_fit <- SplinesMDS(dis_mat, lambda = 0, P = ncol(df)-2, tvec = t_uniq)

# coordinates
tgrid <-seq(min(t_uniq), max(t_uniq), by = min(diff(t_uniq)))

xi1 <- mds_fit$xi1
xi2 <- mds_fit$xi2
init_coords <- mds_fit$init_coord
Xmat <- mds_fit$Xmat

# for one variable
par(mfrow = c(2, 2))
for(p in 1:4){
  
  coord_p1 <- init_coords[p, 1]+ Xmat %*% xi1[p, ]
  coord_p2 <- init_coords[p, 2]+ Xmat %*% xi2[p, ]
  plot(coord_p1, coord_p2, type = "b")
}

coord_p1 <- init_coords[p, 1]+ Xmat %*% xi1[p, ]
coord_p2 <- init_coords[p, 2]+ Xmat %*% xi2[p, ]
plot(coord_p1, coord_p2, type = "l")

plot(tgrid, coord_p1, type = "l")
plot(tgrid, coord_p2, type = "l")

plot(tgrid, coord_p1, type = "l")

coord_t <- cbind(init_coords[, 1] + xi1 %*% Xmat[tid,],
                 init_coords[, 2] + xi2 %*% Xmat[tid,])
}}

