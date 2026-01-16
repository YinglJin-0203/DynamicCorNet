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

#### Dynamic MDS with euclidean distance ####

df_org <- read.csv("Data/IFEDDemoData.csv")
df <- df_org %>% select(-ID, -Length, -Weight.for.age, -Height.for.age, 
                    -Weight.for.height, -Length, -Bud.bead.diameter,
                    -Age.at.exam) %>% 
  rename(time=Week)

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
  
  # plot dynamic
  # pic1 <- paste0("images/CaseStudy/DynMDS_edist_week", tuniq[t], ".jpeg")
  pic1 <- paste0("images/CaseStudy/DynMDS_euclidean_week", tuniq[t], ".jpeg")
  jpeg(filename = pic1, height = 500, width = 500)
  plot(graph_t,
       layout = as.matrix(coord_t),
       # vertex.frame.color=ifelse(is.na(V(graph_t)$color), "grey", NA),
       vertex.frame.color=NA,
       vertex.label.cex=1,
       vertex.size = 20, 
       # vertex.color = V(graph_t)$color, 
       margin = 0, main = paste0("Dynamic, Euclidean, week ", tuniq[t]))
  dev.off()
  fig_list[[t]] <- image_read(pic1)
  
}

# save graphs
DynMDS_educlidean_all <- image_animate(image_join(fig_list), fps = 1)
# DynMDS_edist_noultra <- image_animate(image_join(fig_list[c(1, 2, 4, 6, 8, 12)]), fps = 1)
# DynMDS_edist_ultra <- image_animate(image_join(fig_list[c(3, 5, 7, 9:11)]), fps = 1)

image_write(DynMDS_educlidean_all, path = "images/CaseStudy/DynMDS_euclidean_all.gif")
# image_write(DynMDS_edist_noultra, path = "images/CaseStudy/DynMDS_edist_noultra.gif")
# image_write(DynMDS_edist_ultra, path = "images/CaseStudy/DynMDS_edist_ultra.gif")


#### Splines MDS ####

df_org <- read.csv("Data/IFEDDemoData.csv")
df <- df_org %>% select(-ID, -Length, -Weight.for.age, -Height.for.age, 
                    -Weight.for.height, -Length, -Bud.bead.diameter,
                    -Week) %>% 
  rename(time=Age.at.exam)

tuniq <- sort(unique(df$time))

source(here("Code/SplDissmMat.R"))
source(here("Code/Helpers/SplMDSHelpers.R"))
source(here("Code/SplinesMDS.R"))

# dissmilarity matrix
dist_mats <- SplDissimMat(df, method = "euclidean")
warnings()
dim(dist_mats[[1]])

# coordinates
system.time({
  coords <- SplinesMDS(dist_mats, lambda = 5, P = 12, tvec = tuniq)
})

save(coords, file = "CaseStudy/Splcoords_euclidean.RData")

# graph
load("CaseStudy/Splcoords.RData")
## coeeficients
xi1 <- coords$xi1
xi2 <- coords$xi2
init_coords <- coords$init_coord
  
  # c1 <- init_coord[,1] + xi1 %*% t(Xmat)
  # c2 <- init_coord[,2] + xi2 %*% t(Xmat)
  # coords <- lapply(tid, function(x){return(data.frame(c1 = c1[ , x], 
  #                                                     c2 = c2[ , x]))})

tuniq <- sort(unique(df$time))
tmin <- min(tuniq)
tmax <- max(tuniq)
tid <- tmin:tmax
fig_list <- list()
varnames <- colnames(df %>% select(-time))
P <- 12
K <- ncol(xi1)
Xmat <- bSpline(tmin:tmax, 
                df = K,
                degree = 2, derivs = 0)


# find the median of days at each week
tid <- df_org %>% select(Week, Age.at.exam) %>% 
  group_by(Week) %>%
  summarise_at("Age.at.exam", median) %>% select(Age.at.exam) %>% unlist()

for(t in 1:length(tid)){
  
  # number of non-missing variables
  graph_t <- make_empty_graph(n = P, directed = FALSE)
  V(graph_t)$label <- varnames
  # coord_t <- cbind(init_coords[, 1] + xi1 %*% Xmat[t,],
  #                  init_coords[, 2] + xi2 %*% Xmat[t,]
  #                  )
  coord_t <- cbind(init_coords[, 1] + xi1 %*% Xmat[tid[t],],
                   init_coords[, 2] + xi2 %*% Xmat[tid[t],]
  )
  rownames(coord_t) <- varnames
  
  if(tid[t] %in% tuniq){
    NAcols <- df %>% filter(time == tid[t]) %>% select(-time) %>%
      select(where(~ all(is.na(.)))) %>% colnames()
    frame_color <- ifelse(V(graph_t)$label %in% NAcols, "grey", NA)
    fille_color <- ifelse(V(graph_t)$label %in% NAcols, NA, "grey")

  } else {
    frame_color <- "grey"
    fille_color <- NA
  }
    
  # plot dynamic
  pic1 <- paste0("images/CaseStudy/SplMDS_edist_day", tid[t], ".jpeg")
  # pic1 <- paste0("images/CaseStudy/temp.jpeg")
  jpeg(filename = pic1, height = 500, width = 500)
  plot(graph_t,
       layout = coord_t,
       vertex.frame.color=frame_color,
       vertex.color=fille_color,
       vertex.label.cex=1,
       vertex.size = 20, 
       # vertex.color = V(graph_t)$color, 
       margin = 0, main = paste0("Splines, day ", tid[t]))
  dev.off()
  fig_list[[t]] <- image_read(pic1)
  
}

spl_edist_all <- image_animate(image_join(fig_list), fps = 1)
image_write(spl_edist_all, path = "images/CaseStudy/SplMDS_edist_all.gif")
