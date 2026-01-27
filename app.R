# This shiny app is built for network visualization 
# of temporal multivariate data

# rm(list=ls())
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

theme_set(theme_minimal())

set.seed(825)


#### Helper functions #### 

source(here("Code/Stress.R"))

source(here("Code/DynDissmMat.R"))
source(here("Code/Helpers/DynMDSHelpers.R"))
source(here("Code/DynamicMDS.R"))

source(here("Code/SplDissmMat.R"))
source(here("Code/Helpers/SplMDSHelpers.R"))
source(here("Code/SplinesMDS.R"))

source(here("Code/HclustCoord.R"))
source(here("Code/GroupSumFigure.R"))

source(here("Code/IntDissMat.R"))

#### User interface ####

# UI includes the following elements
# visualization threshold of correlation
ui <- navbarPage(title = "Temporal network visualization of multidimensional data",
  
  # tab 1: data upload and prespecifications
  tabPanel(title = "Data upload and preview",
           
           sidebarLayout(
             # side bar: upload
             sidebarPanel(
               # upload file
               fileInput(inputId = "df_path", label = "Upload a long format csv file",
                         accept = c(".csv")),
               tagList(
                 icon("info-circle"),
                 em("long format data requires each row to represent an observation at each time point. 
                    Multiple observations at the same time point should be stacked vertically.")
               ),
               br(), br(),
               # specify subject ID and time
               uiOutput("time_var"),
               uiOutput("id_var"),
               width = 3
             ),
               
             
             # main panel: data preview
             mainPanel(h3('Data preview'),
                       dataTableOutput("show_df"),
                       h4("Sample summary"), 
                       htmlOutput("size_info"))
    
           )
           ),
  
  # tab 2: descriptive
  tabPanel(title = "Descriptives",
           tabsetPanel(
              ## subtab 2.1: single variable distribution
                tabPanel(title = "Univariate",
                         # side bar 1
                         sidebarLayout(
                           sidebarPanel(
                           uiOutput("varnames1"),
                           # time type with tooltip
                           selectInput("time_type", "Treat time as", choices = c("Continuous", "Discrete"), selected = "Discrete"),
                           # show detailed missing summary?
                           checkboxInput("show_miss_npct", 
                                         "Show the number and percentage of missing values?",
                                         value=T)
                           ),
  
                         # main panel
                         mainPanel(# summary for selected variable
                                   h3('Single variable summary'),
                                   plotOutput("sum_tb"),
                                   # downloadButton("download_sum", "Download"),
                                   h5(icon("circle-info"), "Notes on Individual Variable Summary"),
                                   htmlOutput("sum_tb_note"),
                                   br(),
                                   h3("Missing pattern of each subject"),
                                   plotOutput("miss_plot"),
                                   h3("N (%) of missing observations"),
                                   htmlOutput("miss_note"),
                                   br(),
                                   dataTableOutput("miss_npct_tb", width = "70%"))
              )),
              ## subtab 2.2: pairwise
              tabPanel(title = "Pairwise",
                       # side bar
                       sidebarLayout(
                        sidebarPanel( # correlation type
                          selectInput("cor_type", label="Type of correlation",
                                      choices = list("pearson", "spearman")),
                          tagList(
                            icon("info-circle"),
                            em("Correlation measures may be unrealiable when the proportion of missing is large!")
                            ),
                          br(), br(),
                          checkboxInput("scaleY", "Scale correlation axis to data?", value = F),
                          uiOutput("varnames2")),
                       
                       # main panel
                         mainPanel(
                           h3("Empirical correlation for a pair of variables"),
                           plotOutput("cor_trend_p"),
                           h5(icon("circle-info"), "Notes on correlation summary plot:"),
                           htmlOutput("cor_trend_note"),
                           br(),
                           h3("Comparision of distribution and temporal trend"),
                           plotOutput("trend_p")
                         ))),
              ## subtab 2.3: overall
              tabPanel(title = "Overall",
                         sidebarLayout(
                           sidebarPanel(
                             selectInput("cor_type2", label="Type of correlation",
                                         choices = list("pearson", "spearman")),
                             tagList(
                               icon("info-circle"),
                               em("Correlation measures may be unrealiable when the proportion of missing is large!")
                             ),
                             br(), br(),
                             uiOutput("time_bar1")),
                           mainPanel(
                             h3("Correlation heatmap"),
                             plotOutput("heatmap")
                           )
                         )
                       )
    )),
  
  # tab 3: network and subnetwork plots
  tabPanel(title = "Temporal visualization",
           sidebarLayout(
             # side bar
             sidebarPanel(
               # choose correlation type
               selectInput("mtype", label="Type of correlation/association", 
                           choices = list("pearson", "spearman", "euclidean")),
               tagList(
                 icon("info-circle"),
                 em("Strong association is indicated by large absolute values of correlation, or small values of euclidean distance.")
               ),
               br(),
               tagList(
                 icon("info-circle"),
                 em("While correlation always has a range of [0, 1], the range of euclidean distance is data-dependent.")
               ),
               br(), br(),
               uiOutput("thres_m"),
               # time bar
               uiOutput("time_bar"),
               # hierarchical grouping
               checkboxInput("hclust", label = "Show groups", 
                             value = FALSE),
               uiOutput("nclust"),
               uiOutput("group_sum"),
               uiOutput("group_type"), 
               # variable list
               uiOutput("varnames3"),
               actionButton("confirm", "Confirm selection")
             ),
             
             # main panel
             mainPanel(h3("Temporal network of correlation"),
                       plotOutput("netp", width = "100%", height = "400px"),
                       h5(icon("circle-info"), "Note on the Temporal Network Graph:"),
                       htmlOutput("mds_note1"),
                       h5(icon("circle-info"), "Note on the Group Label Assignments:"),
                       htmlOutput("mds_note2"),
                       h3("Variable group summary:"),
                       div(align = "center", plotOutput("group_plot", width = "80%", height = "300px")),
                       h5(icon("circle-info"), "Note of Group Summary Plots"),
                       htmlOutput("group_note")
                       )
             )),
  
  # tab 4: Overall structure
  tabPanel(title = "Integrated visualization",
           sidebarLayout(
             sidebarPanel(
               # choose correlation type
               selectInput("mtype2", label="Type of correlation/association", 
                           choices = list("pearson", "spearman", "euclidean")),
               uiOutput("thres_m2"),
               # weight
               checkboxInput("time_wt", label = "Weigh by time interval?", value = T),
               tagList(
                 icon("info-circle"),
                 em("Correlation at each time point is weight by the interval between itself and the next time point. 
                    It assumes the data structure stays unchanged until the next measurement point.")
               ),
               # hierarchical grouping
               checkboxInput("hclust2", label = "Show groups", 
                             value = FALSE),
               uiOutput("nclust2"),
               # variable list
               uiOutput("varnames4"),
               actionButton("confirm2", "Confirm selection")
              ),
             mainPanel(
               h3("Integrated network of correlation"),
               plotOutput("int_net"),
               h3("Integrated variable group summary:"),
               plotOutput("int_tree"),
               h5(icon("circle-info"), "Note on Integrated Network Plot and Group Label Assignment"),
               htmlOutput("int_note")
             )
           ))
  
)
  
  
#### Server ####

server <- function(input, output) {
  options(shiny.maxRequestSize=10*1024^2)
  # tab 1: data upload and preview
  ## data upload
  df <- reactive({
    req(input$df_path)
    read.csv(input$df_path$datapath, check.names = F)
  })
  ## data preview  
  output$show_df <- renderDataTable({df()},
    options = list(scrollX = T, fixedHeader=T)
  )
  output$size_info <- renderPrint({
    req(df(), input$time_var, input$id_var)
    Nsize <- length(unique(df()[, input$id_var]))
    Nfreq <- mean(table(df()[, input$id_var]))
    Trange <- range(df()[, input$time_var])
    cat(paste0("Number of participants: ", round(Nsize, 0), "<br>"))
    cat(paste0("Average number of observations per participant: ", round(Nfreq, 2), "<br>"))
    cat(paste0("Range of time: ", round(Trange[1], 2), " - ", round(Trange[2], 2), "<br>"))
    })
  ## specifying time and ID
  output$time_var <- renderUI({
    req(df())
    selectInput("time_var", label = "Specify time", choices = colnames(df()))
  })
  output$id_var <- renderUI({
    req(df())
    selectInput("id_var", label = "Specify participant ID", choices = colnames(df()))
  })
  
  # tab 2: Descriptives
  ## subtab 2.1: single variable distribution
  output$varnames1 <- renderUI({
    req(df())
    selectInput("select_var1", label = "Variables", choices = colnames(df()), 
                       selected = colnames(df())[5])
  })
  
  ## summary of single variables
  plot_sum <- reactive({
    req(df(), input$select_var1, input$time_var, input$id_var)
    ### summary plot
    if(input$time_type=="Discrete"){
      df_sum <- df()[, c(input$id_var, input$time_var, input$select_var1)] %>%
        rename(time=input$time_var, id=input$id_var, var=input$select_var1) %>%
        group_by(time) %>%
        mutate(med=median(var, na.rm = T)) %>% 
        mutate(Nmiss = sum(is.na(var)), 
               Pctmiss = sum(is.na(var))/length(var))
      t_uniq <- sort(unique(df_sum$time))
      # xlab <- df_sum %>% select(time, Nmiss, Pctmiss) %>% distinct(.) %>%
      #   mutate(Nmiss = paste0(Nmiss, " (", round(100*Pctmiss, 2), "%)"))
      # summary plot
      plot_sum <- df_sum %>% 
        ggplot()+
        geom_boxplot(aes(x=time, y=var, group=time), outlier.size = 0.5, fill = "grey")+
        geom_jitter(aes(x=time, y=var, group=time), size = 0.5)+
        geom_line(data = df_sum %>% filter(!is.na(med)), aes(x=time, y=med))+
        scale_x_continuous(breaks = t_uniq, name = input$time_var)+
                           # sec.axis = sec_axis(~., name = "N (pct) of missing", breaks = xlab$time, label=xlab$Nmiss))+
        # theme(axis.text.x.top = element_text(angle=90))+
        labs(x=input$time_var, y=input$select_var1, 
             title = paste0("Distribution and temporal trend of ", input$select_var1))
    }
    else{
      plot_sum <- df()[, c(input$id_var, input$time_var, input$select_var1)] %>%
        rename(time=input$time_var, id=input$id_var, var=input$select_var1) %>%
        filter(complete.cases(.)) %>%
        ggplot()+
        geom_line(aes(x=time, y=var, group=id), alpha = 0.5, linewidth = 0.5)+
        geom_smooth(aes(x=time, y=var), method = gam, formula = y~s(x))+
        labs(x=input$time_var, y = input$select_var1,
             title = paste0("Distribution and temporal trend of ", input$select_var1))
    }
    plot_sum
    })
  output$sum_tb <- renderPlot({
    plot_sum()
  })
  # missing plot
  output$miss_plot <- renderPlot({
    req(df(), input$select_var1, input$time_var, input$id_var)
    # uni_id <- unique(df[ ,input$id_var])
   df_miss <- df()[, c(input$id_var, input$time_var, input$select_var1)] %>%
      rename(time=input$time_var, id=input$id_var, var=input$select_var1) %>%
      mutate(id = as.factor(id)) %>%
      arrange(time)
   if(input$time_type=="Discrete"){
      df_miss %>%
       pivot_wider(id_cols = "id", names_from = "time", values_from = "var") %>%
       select(-id) %>%
       visdat::vis_miss(.)+
       labs(x=paste0(input$time_var, " (% present)"), y = "ID")
   } else{
     df_miss %>% ggplot()+
       geom_tile(aes(x=time, y=id, fill = is.na(var)))+
       geom_line(aes(x=time, y=id), alpha = 0.3, color = "grey20")+
       scale_fill_manual(name = "", values = c("grey80", "grey20"), labels = c("Present", "Missing"))+
       scale_x_continuous(breaks = seq(0, 300, by = 20))+
       theme(legend.position = "bottom", axis.text.y = element_blank())+
       labs(x=input$time_var, y = "ID")
   }
  })
  ## table for missingness
  output$miss_npct_tb <- renderDataTable({
    req(df(), input$select_var1, input$time_var, input$id_var, input$show_miss_npct)
    df_miss <- df()[, c(input$id_var, input$time_var, input$select_var1)] %>%
      rename(time=input$time_var, id=input$id_var, var=input$select_var1) %>%
      mutate(id = as.factor(id)) %>%
      arrange(time) %>%
      group_by(time) %>%
      summarize(N = sum(is.na(var)), 
                Pct = round(sum(is.na(var))/length(var), 2))
    colnames(df_miss) <- c(input$time_var, "N", "%")
    datatable(df_miss) %>%
      formatStyle("N", target = "row", backgroundColor = styleInterval(20, c(NA,"#ffe6e6")))
  })
  ## notes for the boxplot
  output$sum_tb_note <- renderText({
    if(input$time_type=="Discrete"){
      meg <- HTML("
              <li>Each measurement is considered as a occurrence at a discrete time point.</li>
              <li>Distribution of measurements from different subjects at the same time point is represented with boxplot. </li>
              <li>Trend over time is represented by change of median.</li>
                  ")
    } 
    else{
      meg = HTML("
             <li>Measurements is considered noisy realizations of an underlying smooth process over time.</li>
             <li>Each line represents the trajectory from one subject.</li>
             <li>The blue line is a smoothed trend over all subjects and measurement points.</li>
                 ")
    }
    print(meg)
  })
  ## mising value note
  output$miss_note <- renderText({
    HTML("Correlation measure is sensitive to the missing values and can be unrealiable if the propotion of missing value is large.
         If at any time point the total number of observation is less than 20, 
         we may consider the correlation calculated at this time point too unreliable and unfit for further analysis.
         We recommend treating such correlation as missing at this time point.")
  })
  
  ## subtab 2.2: pairwise correlaion
  output$varnames2 <- renderUI({
    req(df())
    checkboxGroupInput("select_var2", label = "Variables", choices = colnames(df()), 
                       selected = colnames(df())[3:4])
  })
  ### empirical correlations
  output$cor_trend_p <- renderPlot({
    req(df(), input$select_var2, input$time_var, input$id_var, input$time_type)
    validate(need(length(input$select_var2)==2, "Please select a pair of variables."))
    # sub data
    df_pair <- df()[, c(input$time_var, input$id_var, input$select_var2)] %>%
      rename(time=input$time_var, id = input$id_var)
    t_uniq <- unique(df_pair$time)
    N <- length(unique(df()[ , input$id_var]))
    # correlation and proportion of complete pairs
    var1 <- input$select_var2[1]
    var2 <- input$select_var2[2]
    df_cor <- df_pair %>% 
      group_by(time) %>%
      summarize(cor = cor(.data[[var1]], .data[[var2]], 
                          use = "pairwise.complete.obs", method = input$cor_type),
                Npair = sum(complete.cases(.data[[var1]], .data[[var2]]))/N)
    # plot
    p2 <- df_cor %>% 
      filter(complete.cases(.)) %>%
      ggplot()+
      geom_point(aes(x=time, y=cor, size = Npair))+
      geom_line(aes(x=time, y=cor))+
      labs(title = "Empirical correlation", x = input$time_var, y = " ", size = "Proportion of complete pairs")+
      theme(legend.position = "bottom")+
      guides(color = guide_legend(order = 1), 
             size = guide_legend(order = 1))
    if(!input$scaleY){p2 <- p2 + ylim(-1, 1)}# scale correlation axis
    if(input$time_type == "Discrete"){p2 <- p2 + scale_x_continuous(breaks = t_uniq)}
    # display
    p2
  })
  ###### note
  output$cor_trend_note <- renderText({
    HTML("
    <li>The node size indicates the number of complete pairs used to calculated correlation at each time point.</li>
    <li>Correlation/association measures are very sensitive to sample size. If the number of complete pais is very small (i.e < 20), 
    the calculated measures are less realiable and will affect downstream analysis.
         User may consider removing these time points from the dataset </li>")
  })
  
  #### comparision of distribution and trend
  output$trend_p <- renderPlot({
    req(df(), input$select_var2, input$time_var, input$id_var, input$time_type)
    validate(need(length(input$select_var2)==2, "Please select a pair of variables."))
    df_pair <- df()[, c(input$time_var, input$id_var, input$select_var2)] %>%
      rename(time=input$time_var, id = input$id_var)
    t_uniq <- unique(df_pair$time)
    # trend plot
    if(input$time_type == "Discrete"){
      p1 <- df_pair %>%
        pivot_longer(input$select_var2) %>%
        ggplot(aes(x = time, y=value, fill=name, colour = name, group=interaction(time, name)))+
        geom_boxplot(position = "dodge2", alpha = 0.7)+
        geom_jitter(size = 0.5)+
        scale_fill_brewer(palette = "Set2")+
        scale_color_brewer(palette = "Set2")+
        scale_x_continuous(breaks = t_uniq)+
        labs(title = "Variable distribution", x = input$time_var, y = " ")+
        theme(legend.position = "bottom")
    }
    else {
      p1 <- df_pair %>%
        pivot_longer(input$select_var2) %>%
        filter(complete.cases(.)) %>%
        ggplot(aes(x = time, y=value, colour = name, group=interaction(id, name)))+
        geom_line(alpha = 0.7, linewidth=0.5)+
        scale_fill_brewer(palette = "Set2")+
        scale_color_brewer(palette = "Set2")+
        # scale_x_continuous(breaks = t_brk)+
        labs(title = "Variable trend", x = input$time_var, y = " ")+
        theme(legend.position = "bottom")
    }
    # display
    p1
  })
  
  # subtab 2.3: overall correlation heatmap
  output$time_bar1 <- renderUI({
    req(df(), input$time_var, input$id_var, input$time_type)
    # time bar
    tvec <- sort(unique((df()[ ,input$time_var])))
    if(input$time_type == "Discrete"){
      sliderTextInput("time_bar1", "Time of heatmap", choices = tvec, selected = tvec[1],
                      grid = TRUE)
    }
    else {
      sliderTextInput("time_bar1", "Time of heatmap", choices = tvec, selected = tvec[1],
                      grid = FALSE)
    }
    
  })
  output$heatmap <- renderPlot({
    req(df(), input$time_var, input$id_var, input$cor_type2)
    ## correlation matrix
    cormat <- cor(subset(df() %>% rename(time=input$time_var, id=input$id_var) %>%
                           filter(time==input$time_bar1), 
                         select = -c(id, time)), method = input$cor_type2, use = "pairwise.complete.obs")
    ## heatmap
    data.frame(cormat) %>% rownames_to_column("var1") %>%
      pivot_longer(-var1) %>%
      ggplot(aes(x=var1, y=name, fill = value))+
      geom_tile()+
      scale_fill_continuous_diverging(palette = 'Blue-Red 3', mid = 0, breaks = c(-1, 0, 1), limits = c(-1, 1),
                                      rev = TRUE)+
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle=90),
            aspect.ratio = 1)+
      labs(x="", y="", fill="Correlation", title = "Correlation matrix")
  }, height = 500, width = 500)
  
  # tab 3: temporal network and group labels
  ## variable list
  output$varnames3 <- renderUI({
    req(df(), input$time_var, input$id_var)
    all_vars <- colnames(df() %>% select(!c(input$time_var, input$id_var)))
    checkboxGroupInput("select_var3", label = "Variables", choices = all_vars, select = all_vars)
  })
  ## threshold
  output$thres_m <- renderUI({
    req(input$mtype)
    if(input$mtype == "euclidean"){
      diss_range <- round(range(unlist(diss_mat()), na.rm = T))
      sliderInput("thres_m", label = "Show distance below",
                  min = diss_range[1], max = diss_range[2], value = diss_range[1])
    } else {
      sliderInput("thres_m", label = "Show correlation above",
                  min =0, max = 1, value = 1)
    }
  })
  ## data set to analyze
  confirmed <- reactiveVal(NULL)
  observeEvent(input$confirm, {confirmed(input$select_var3)})
  df_net <- reactive({
    req(df(), input$id_var, input$time_var, confirmed())
    df()[, c(input$time_var, confirmed())] %>%
      rename(time = input$time_var) %>%
      filter(!if_all(confirmed(), is.na)) # remove empty columns
  })
  ## time axis
  output$time_bar <- renderUI({
    req(df_net(), input$time_type, input$time_var)
    # time bar: by the original time
    tvec <- sort(unique(df_net()[, "time"]))
   if(input$time_type=="Discrete"){
      time_bar <- sliderTextInput("time_bar", label = input$time_var, choices = tvec, 
                                  selected = tvec[1], grid = TRUE)
    } else {
      time_bar <- sliderInput("time_bar", label = input$time_var, 
                              min = min(tvec, na.rm = T), max = max(tvec, na.rm = T),
                              step = min(diff(tvec), na.rm = T), value = min(tvec, na.rm = T),
                              ticks = FALSE)
    }
    time_bar
  })
  ### time grids
  
  ## grouping results
  output$nclust <- renderUI({
    req(input$hclust)
    numericInput("nclust", label = "Number of groups", value = 1)
  })
  output$group_sum <- renderUI({
    req(input$hclust, input$nclust)
    checkboxInput("group_sum", "Show summary of grouping over time")
    
  })
  output$group_type <- renderUI({
    req(input$group_sum)
    radioButtons("group_plot", label = "Group summary: ", choices = list("By temporal flow" = 1, 
                                                                         "By variable" = 2,
                                                                         "Hierarchical tree" = 3))
  })
  
  # main panel outputs
  ## calculate dissimilarity matrix at each time point
  diss_mat <- reactive({
    req(df_net(), input$time_type, input$mtype)
    if(input$time_type == "Discrete"){
           DynDissimMat(df_net(), method = input$mtype)
      } else {
           SplDissimMat(df_net(), method = input$mtype)
  }
  })
  ## calculate adjacency matrix from dissimilarity matrix
  adj_mat <- reactive({
    req(diss_mat(), input$mtype, input$thres_m)
    # adjacency matrix
    if(input$mtype=="euclidean"){
      max_diss <- max(unlist(diss_mat()), na.rm = T)
      min_diss <- min(unlist(diss_mat()), na.rm = T)
      adj_mat <- lapply(diss_mat(),
                     function(diss_t){(max_diss-diss_t)/(max_diss-min_diss)})
      thres <- (max_diss-input$thres_m)/(max_diss-min_diss)
      adj_mat <- lapply(adj_mat, 
                        function(adj_t){
                          adj_t[adj_t <= thres] <- 0
                          return(adj_t)})
    } else {
      adj_mat <- lapply(diss_mat(), 
                        function(diss_t){
                          adj_t <- 1-diss_t
                          adj_t[adj_t <= input$thres_m] <- 0
                          return(adj_t)
                        })
    }
    adj_mat
  })
  ## calculate coordinates
  coord_list <- reactive({
    req(diss_mat(), df_net(), input$time_type)
    t_uniq <- sort(unique(df_net()[, "time"])) # original time scale
    # t_id <- seq_along(t_uniq) # time index
    if(input$time_type == "Discrete"){
           DynamicMDS(diss_mat(), lambda = 5)
    } else {
        SplinesMDS(diss_mat(), lambda = 7, P = ncol(df_net())-1, tvec = t_uniq)
      # in this case, coord_list includes initial coordinates and cofficients
  }
  })
  ## grouping
  group_list <- reactive({
    req(input$hclust, coord_list())
      if(input$time_type == "Discrete"){
        group_list <- HclustCoord(coord_list(), "dynamic")
      } else {
        group_list <- HclustCoord(coord_list(), "splines")
      }
    group_list
  })
  ## plot
  output$netp <- renderPlot({

      withProgress(
        value=0, message = "Processing", detail="This may take a while...",
        {
          req(df_net(),  coord_list(), input$time_bar, adj_mat())
          t_uniq <- sort(unique(df_net()[, "time"]))
    
          # Dynamic MDS 
          if(input$time_type == "Discrete"){
            tid <- which(input$time_bar == t_uniq)
            # diss_t <- diss_mat()[[tid]]
            coord_t <- as.matrix(coord_list()[[tid]])
            adj_t <- adj_mat()[[tid]]
            nodes_t <- colnames(adj_t)
            vars_t <- colnames(adj_t)
            # graph
            net_t <- graph_from_adjacency_matrix(adj_t, weighted = T, 
                                                 mode = "undirected", 
                                                 diag = FALSE)
            # E(net_t)$label <- round(E(net_t)$weight, 2)
          } else {
            # Splines MDS
            tgrid <-seq(min(t_uniq), max(t_uniq), by = min(diff(t_uniq)))
            tid <- which(tgrid == input$time_bar)
            xi1 <- coord_list()$xi1
            xi2 <- coord_list()$xi2
            init_coords <- coord_list()$init_coord
            Xmat <- coord_list()$Xmat
            # splines 
            coord_t <- cbind(init_coords[, 1] + xi1 %*% Xmat[tid,],
                             init_coords[, 2] + xi2 %*% Xmat[tid,])
            # all nodes on the graph
            nodes_t <- colnames(df_net() %>% select(!time))
            # non-empty variables
            vars_t <- colnames(df_net() %>% filter(time==input$time_bar) %>%
                                 select(!time) %>%
                                 select(where(~!all(is.na(.)))))
            # graph and edges
            if(input$time_bar %in% t_uniq){
              # adjacency matrix
              adj_t <- adj_mat()[input$time_bar == t_uniq][[1]]
              net_t <- graph_from_adjacency_matrix(adj_t, weighted = T, 
                                                   mode = "undirected", 
                                                   diag = FALSE)
            } else {
              net_t <- make_empty_graph(n = length(nodes_t), directed = FALSE)
              V(net_t)$name <- nodes_t
              }
            }
          
          # plot properties
          V(net_t)$color <- rgb(0.2, 0.4, 0.8, alpha=0.4)
          E(net_t)$width = 7*E(net_t)$weight
          
          # if did clustering, color by cluster
          if(input$hclust){
            hclust_t <- group_list()[[tid]]
            group_t <- cutree(hclust_t, input$nclust)
            node_col_t <- brewer.pal(input$nclust, "Accent")[group_t]
            V(net_t)$color <- node_col_t
          }
          
          # plot
          par(
            mar = c(5, 4, 4, 2) + 0.1,  # extra space at bottom
            xpd = NA                   # allow drawing outside plot region
          )
          plot(net_t,
           layout = coord_t,
           vertex.frame.color=V(net_t)$color,
           vertex.color = ifelse(V(net_t)$name %in% vars_t, V(net_t)$color, NA),
           vertex.label.cex=1,
           vertex.size = 20,
           margin = 0)
          # legend
          if(input$time_type == "Continuous"){
            legend(
              "bottom",
              inset = -0.15, 
              legend = c("Present variables", "Unmeasured variables"),
              pch = 21,
              pt.bg = c(rgb(0.2, 0.4, 0.8, alpha=0.4), "white"),
              pt.cex = 2,
              col = rgb(0.2, 0.4, 0.8, alpha=0.4),
              horiz = TRUE,
              ncol = 2,
              bty = "n"
            )
          }
          
    incProgress(1)}
    )
    })
  # note
  output$mds_note1 <- renderPrint({
    mds_type <- ifelse(input$time_type=="Discrete", "Dynamic", "Splines")
    HTML(paste0("This networkplot is generated by ", mds_type, " multidimentional scaling based on ", input$time_type, " time and ", input$cor_type, " correlation.
            To change the type of time variable or correlation, please move back to the previous tab."))
  })
  output$mds_note2 <- renderPrint({
    HTML("Group labels are estiamted based on a hierarchical clustering model fit on 2D coordinates. ")
  })
    ## group label plot
  output$group_plot <- renderPlot({
    req(input$hclust, input$nclust, input$group_sum, group_list(), df_net(), input$group_plot, input$time_bar)
    
    # time scale
    t_uniq <- sort(unique(df_net()[,"time"]))
    if(input$time_type == "Discrete"){
      tgrid <- t_uniq
    } else {
      tgrid <-seq(min(t_uniq), max(t_uniq), by = min(diff(t_uniq)))
    }
    # plot
    if(input$group_plot != 3){
      sum_plot <- GroupSumFigure(group_list(), input$nclust, tgrid, plot_method = input$group_plot)
      sum_plot <- sum_plot+
        geom_vline(xintercept = input$time_bar)
        
    } else{
       sum_plot <- ggdendrogram(group_list()[tgrid==input$time_bar][[1]],
                    rotate = T, size = 2)+
            labs(title = paste0("Hierarchial group at time ", input$time_bar))

    }
    sum_plot
  })
  ## note
  output$group_note <- renderText({
    req(input$hclust, input$nclust, input$group_sum, input$group_plot)
    if(input$group_plot==1){
      HTML("
        <li>This plot visualizes the change of group structure of variables over time using an <a href='https://corybrunson.github.io/ggalluvial/' target='_blank'>Alluvial plot</a></li>
        <li>Band with different colors represents different groups, and the width of band represents the size of group.</li>
        <li>Band flow between different groups across time represents variables that moved from one group to the other, and the width of flow represents number of variables that made the switch.</li>
         ")
    } else if(input$group_plot==2){
      HTML("
        <li>This plot visualizes the change of group label assignment for each variable.</li>
        <li>Each row represents a single variable, and the color represents the group it is assigned to at specific time points.</li>
        <li>Change of color indicates change of group label assignments.</li>
         ")
    } else{
      HTML("
        <li>This plot visualizes the hierachical structure of variables at a specific time.</li>
         ")}
    })
  
 
  # tab 4: integrated correlation and grouping results
  ## sidebar variable list
  output$varnames4 <- renderUI({
    req(df(), input$time_var, input$id_var)
    all_vars <- colnames(df() %>% select(!c(input$time_var, input$id_var)))
    checkboxGroupInput("select_var4", label = "Variables", 
                       choices = all_vars, selected = all_vars)
  })
  confirmed2 <- reactiveVal(NULL)
  observeEvent(input$confirm2, {confirmed2(input$select_var4)})
  # dataset for analysis
  df_net2 <- reactive({
    req(df(), input$id_var, input$time_var, confirmed2())
    df()[, c(input$time_var, confirmed2())] %>%
      rename(time = input$time_var) %>%
      filter(!if_all(confirmed2(), is.na)) # remove empty columns
  })
  # integrated dissimilarity matrix
  int_diss <- reactive({
    req(df_net2(), input$mtype2)
    int_diss <- IntDissmMat(df_net2(), method = input$mtype2, weight = input$time_wt)
    int_diss
  })
  # integrated adjacency matrix from dissimilarity matrix
 int_adj <- reactive({
    req(int_diss(), input$mtype2, input$thres_m2)
    # adjacency matrix
    if(input$mtype2=="euclidean"){
      max_diss <- max(int_diss(), na.rm = T)
      min_diss <- min(int_diss(), na.rm = T)
      int_adj <- (max_diss-int_diss())/(max_diss-min_diss)
      thres <- (max_diss-input$thres_m2)/(max_diss-min_diss)
      int_adj[int_adj <= thres] <- 0
    } else {
      int_adj <- 1-int_diss()
      int_adj[int_adj <= input$thres_m2] <- 0
    }
    int_adj
  })
  # correlation/association threshold
  output$thres_m2 <- renderUI({
    req(input$mtype2, int_diss())
    if(input$mtype2 == "euclidean"){
      diss_range <- round(range(int_diss(), na.rm = T))
      sliderInput("thres_m2", label = "Show distance below",
                  min = diss_range[1], max = diss_range[2], value = diss_range[1])
    } else {
      sliderInput("thres_m2", label = "Show correlation above",
                  min =0, max = 1, value = 1)
    }
  })
 # clustering
  output$nclust2 <- renderUI({
    req(input$hclust2)
    numericInput("nclust2", label = "Number of groups", value = 1)
  })
  # ## values for sanity checks: uniform or all missing across all time points
  # sanity_check2 <- reactive({
  #   req(df(), input$id_var, input$time_var, input$select_var4, confirmed2())
  #   unique_val <- df()[ , c(input$time_var, confirmed2())] %>%
  #     rename(time = input$time_var) %>% group_by(time) %>%
  #     summarize_all(~{length(unique(.))}) %>% ungroup() %>% 
  #     summarise(across(everything(), ~all(.x==1, na.tm = T))) %>% 
  #     select(where(isTRUE))
  #   colnames(unique_val)
  # })
  ## coordinates
  int_coords <- reactive({
    req(int_diss())
    mds(int_diss())$conf
  })
  ## hierarchical clustering result
  int_hclust <- reactive({
    req(int_diss())
    hclust(as.dist(int_diss()))
  })
  ## network plot
  output$int_net <- renderPlot({
    req(int_adj(), input$thres_m2, int_coords())
    int_net <- graph_from_adjacency_matrix(int_adj(), mode = "undirected", weighted = T, diag=F)
    # edges
    E(int_net)$width <- E(int_net)$weight*7
    # node properties
    V(int_net)$color <- rgb(0.2, 0.4, 0.8, alpha=0.4)
    coords <- int_coords()
    # if choose to visualize grouping results
    if(input$hclust2){
      clust_group <- cutree(int_hclust(), k = input$nclust2)
      vcolor <- brewer.pal(input$nclust2, "Accent")[clust_group]
      names(vcolor) <- names(clust_group)
      V(int_net)$color <- vcolor[V(int_net)$name]
    }
    plot(int_net, layout=int_coords(),
           vertex.color =  V(int_net)$color, 
           vertex.frame.color=  V(int_net)$color,
           vertex.label.cex=1,
           vertex.size = 20)
  })
  output$int_tree <- renderPlot({
    req(input$hclust2)
    ggdendrogram(int_hclust(), rotate = F, size = 2)+
            labs(title = "")
    
  })
  ## note
  output$int_note <- renderPrint({
    HTML(paste0("Both plots are generated based on the integerated ", input$mtype2, 
                " correlation matrix. To change the type of correlation, please move back to the second tab."))
  })
  
}



# Run the application 
shinyApp(ui = ui, server = server)
