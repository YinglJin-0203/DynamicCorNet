# in this script
# we explore the effect on grouping structure 
# of missing measurements

rm(list=ls())



#### Run MDS on simulated dataset #####

# load data 
df <- read.csv(here("data/SimMissG1.csv"))
# rea data
df <- read.csv(here("data/AppData.csv"))
df <- df %>% rename(id = Participant, time=Week)
df <- df %>% 
  select(-BudBead_A, -starts_with("Tot"))

# adjacency matrix
adj_mat <- GetAdjMat(data= df %>% select(!c("id")), 
          cor_method = "pearson")

# network graph initialization
graph_list <- graph_dyn_net(adj_mat)

# splines MDS coordinates
t_uniq <- sort(unique(df$time))
coord_list <- SplinesMDS(adj_mat, lambda = 10, K = 20, P = dim(adj_mat[[1]])[1], 
                         tvec = t_uniq)


# clustering
group_list <- lapply(adj_mat,
                     function(x){
                       dis_mat <- 1-x
                       # remove fully unmeasured variable
                       var_id <- which(!is.na(diag(dis_mat)))
                       dis_mat <- dis_mat[var_id, var_id]
                       # when there is <2 complete pairwise obsevation
                       dis_mat <- dis_mat[complete.cases(dis_mat), complete.cases(dis_mat)]
                       dis_mat <- as.dist(dis_mat)
                       hclust_fit <- hclust(dis_mat)
                       return(cutree(hclust_fit, k = 2))
                     })

# graph
library(magick)

netplot_list <- list()

for(input_tid in seq_along(t_uniq)){
  
  graph_t <- graph_list[[input_tid]]
  group_t <- group_list[[input_tid]]
  coord_t <- coord_list[[input_tid]]
  
  V(graph_t)$color <- group_t[V(graph_t)$name]
  
  
  # plot
  picname <- paste0("images/ifed_t", input_tid, "c2.jpeg")
  jpeg(filename = picname, height = 600, width = 600)
  plot(graph_t,
       layout = as.matrix(coord_t),
       vertex.frame.color=ifelse(is.na(V(graph_t)$color), "grey", NA),
       vertex.label.cex=1,
       vertex.size = 20, 
       vertex.color = V(graph_t)$color, 
       margin = 0, main = paste0("Time = ", t_uniq[input_tid]))
  dev.off()
  netplot_list[[input_tid]] <- image_read(picname)
  # spearman_list[[input_tid]] <- image_read("images/temp.jpeg")
}

animation <- image_animate(image_join(netplot_list), delay=100)
animation
# Save to file
image_write(animation, path = "images/animation.gif")
image_write(spearman_animation, path = "images/spearman_animation_full.gif")

