library(shiny)
library(shinyWidgets)
library(bslib)
library(shinyBS)
library(bsicons)
library(here)
library(DT)
library(tidyverse)
library(gridExtra)
library(arsenal)
library(htmltools)
library(smacof)
library(splines2)
library(RColorBrewer)
library(ggalluvial)
library(ggdendro)
library(mgcv)
library(igraph)
library(magick)

theme_set(theme_minimal())

set.seed(825)

source(here("Code/Stress.R"))

df <- read.csv("Data/IFEDDemoData.csv")
df <- df %>% select(-ID, -Length, -Weight.for.age, -Height.for.age, 
                    -Weight.for.height, -Length, -Bud.bead.diameter,
                    -Age.at.exam) %>% 
  rename(time=Week)




#### Dynamic MDS with euclidean distance ####
source(here("Code/DynDissmMat.R"))
source(here("Code/Helpers/DynMDSHelpers.R"))
source(here("Code/DynamicMDS.R"))

# dissmilarity matrix
dist_mats <- DynDissimMat(df, method = "euclidean")

# coordinates
system.time({
  coords <- DynamicMDS(dist_mats, 5)
})

# graph
tuniq <- sort(unique(df$time))
fig_list <- list()

for(t in seq_along(tuniq)){
  
  Pt <- dim(dist_mats[[t]])[1]
  graph_t <- make_empty_graph(n = Pt, directed = FALSE)
  V(graph_t)$label <- colnames(dist_mats[[t]])
  coord_t <- coords[[t]]
  
  
  # V(graph_t)$color <- color_vec[rep(1:3, each = 5)] # color by true group
  # if(tvec[t] %in% t_miss){
  #   V(graph_t)$color[1:5] <- NA
  #   coord_t1 <- rbind(matrix(NA, nrow=5, ncol=2), coord_t1)
  # } # for missing variables
  
  # plot dynamic
  pic1 <- paste0("images/CaseStudy/DynMDS_edist_week", tuniq[t], ".jpeg")
  jpeg(filename = pic1, height = 500, width = 500)
  plot(graph_t,
       layout = as.matrix(coord_t),
       # vertex.frame.color=ifelse(is.na(V(graph_t)$color), "grey", NA),
       vertex.frame.color=NA,
       vertex.label.cex=1,
       vertex.size = 20, 
       # vertex.color = V(graph_t)$color, 
       margin = 0, main = paste0("Dynamic, week ", tuniq[t]))
  dev.off()
  fig_list[[t]] <- image_read(pic1)
  
}

# save graphs
DynMDS_edist_all <- image_animate(image_join(fig_list), fps = 1)
DynMDS_edist_noultra <- image_animate(image_join(fig_list[c(1, 2, 4, 6, 8, 12)]), fps = 1)
DynMDS_edist_ultra <- image_animate(image_join(fig_list[c(3, 5, 7, 9:11)]), fps = 1)

image_write(DynMDS_edist_all, path = "images/CaseStudy/DynMDS_edist_all.gif")
image_write(DynMDS_edist_noultra, path = "images/CaseStudy/DynMDS_edist_noultra.gif")
image_write(DynMDS_edist_ultra, path = "images/CaseStudy/DynMDS_edist_ultra.gif")


#### Splines MDS ####

df <- read.csv("Data/IFEDDemoData.csv")
df <- df %>% select(-ID, -Length, -Weight.for.age, -Height.for.age, 
                    -Weight.for.height, -Length, -Bud.bead.diameter,
                    -Week) %>% 
  rename(time=Age.at.exam)

source(here("Code/SplDissmMat.R"))
source(here("Code/Helpers/SplMDSHelpers.R"))
source(here("Code/SplinesMDS.R"))

# dissmilarity matrix
dist_mats <- SplDissimMat(df, method = "euclidean")
dim(dist_mats[[1]])

# coordinates
system.time({
  coords <- SplinesMDS(dist_mats, 7, 12, tuniq)
})

save(coords, file = "CaseStudy/Splcoords.RData")
