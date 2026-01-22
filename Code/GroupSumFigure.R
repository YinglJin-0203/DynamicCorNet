
#' Summary plot of group clustering results over time
#'
#' @param group_list Each element is the hierarchical clustering model at a specific time point
#' @param tvec time grid corresponding to the list
#' @param k number of clusters
#' @param plot_method type of plot
#'
#' @returns
#' @export
#'
#' @examples
#' 
GroupSumFigure <- function(group_list, k, tvec, plot_method){
  
  if(plot_method == 1){ # alluvium plot
    group_label_list <- lapply(group_list, function(c){cutree(c, k)}) # group label
    sum_plot <- bind_rows(group_label_list, .id = "time") %>%
      mutate(time=tvec) %>%
      pivot_longer(-time) %>%
      mutate(value = as.factor(value), time = as.numeric(time)) %>%
      ggplot(aes(x=time, stratum = value, fill=value, color=value, alluvium=name))+
      geom_flow()+
      geom_stratum()+
      labs(x="Time", y = "Group", title = "Group flow chart")+
      guides(fill=guide_legend("Group"), color=guide_legend("Group"))+
      theme(legend.position = "bottom", axis.text.y = element_blank())+
      scale_color_brewer(palette = "Accent")+
      scale_fill_brewer(palette = "Accent")
  } else { # tile plot
    group_label_list <- lapply(group_list, function(c){cutree(c, k)}) # group label
    sum_plot <- bind_rows(group_label_list, .id = "time") %>%
      mutate(time=tvec) %>%
      pivot_longer(-time) %>%
      mutate(value = as.factor(value), time = as.numeric(time)) %>%
      ggplot(aes(x=time, y=name, fill=value))+
      geom_tile(alpha=0.7)+
      labs(x="Time", y = " ", fill = "Group", title = "Variable group assignment")+
      theme(legend.position = "bottom")+
      scale_fill_brewer(palette = "Accent")
  } 
  
  return(sum_plot)
  
}

# test 
# df <- read.csv("Data/IFEDDemoData.csv")
# df <- df %>% select(-ID, -Age.at.exam) %>%
#   rename(time=Week)
# t_uniq <- sort(unique(df$time))
# tvec <-seq(min(t_uniq), max(t_uniq), by = min(diff(t_uniq)))
# # diss_list <- DynDissimMat(df)
# diss_list <- SplDissimMat(df)
# # try_lo <- DynamicMDS(diss_list)
# try_lo <- SplinesMDS(diss_list, 7, 12, sort(unique(df$time)))
# group_list <- HclustCoord(try_lo, method = "splines")
# 
# GroupSumFigure(group_list, 4, tvec, plot_method = 2)
