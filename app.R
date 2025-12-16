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

theme_set(theme_minimal())

set.seed(825)


#### Helper functions #### 

source(here("helpers/Stress.R"))
source(here("helpers/AdjacencyMat.R"))
source(here("helpers/DynNet.R"))
source(here("helpers/DynamicMDS.R"))
source(here("helpers/SplinesMDS.R"))

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
               # uiOutput("bin_yn"),
               # uiOutput("bin_width"),
               uiOutput("id_var")
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
              ## sub-tab 1: single variable distribution
                tabPanel(title = "Univariate",
                         # side bar 1
                         sidebarLayout(
                           sidebarPanel(
                           uiOutput("varnames1"),
                           # time type with tooltip
                           selectInput("time_type", "Treat time as", choices = c("Continuous", "Discrete"), selected = "Discrete")
                           ),
  
                         # main panel
                         mainPanel(# summary for selected variable
                                   h3('Single variable summary'),
                                   plotOutput("sum_tb"),
                                   h4("Note:"), 
                                   htmlOutput("sum_tb_note")
                                             )
              )),
              ## sub-tab 2: pairwise
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
                          uiOutput("varnames2")),
                       
                       # main panel
                         mainPanel(
                           h3("Empirical correlation for a pair of variables"),
                           plotOutput("trendp")
                         ))),
              ## subtab 3: all variables
              tabPanel(title = "Overall",
                         sidebarLayout(
                           sidebarPanel(uiOutput("time_bar1")),
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
               # choose correlation threshold
               sliderInput(inputId = "thres_cor", label = "Show correlation above:", min=0, max=1, value=1, step=0.01,
                           ticks=FALSE),
               # time bar
               uiOutput("time_bar"),
               # hierarchical grouping
               checkboxInput("hclust", label = "Show groups", 
                             value = FALSE),
               uiOutput("nclust"),
               uiOutput("group_type"), 
               # variable list
               uiOutput("varnames3"),
               actionButton("confirm", "Confirm selection")
             ),
             
             # main panel
             mainPanel(h3("Temporal network of correlation"),
                       plotOutput("netp", width = "100%", height = "400px"),
                       h4("Note:"),
                       htmlOutput("mds_note1"),
                       htmlOutput("mds_note2"),
                       h3("Variable group summary:"),
                       div(align = "center", plotOutput("group_plot", width = "80%", height = "300px")),
                       h4("Note:"),
                       htmlOutput("group_note")
                       )
             )),
  
  # tab 4: Overall structure
  tabPanel(title = "Integrated visualization",
           sidebarLayout(
             sidebarPanel(
               sliderInput(inputId = "thres_cor2", 
                           label = "Show correlation above:", min=0, max=1, value=1, step=0.01,
                           ticks=FALSE),
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
               h4("Note:"),
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
  output$sum_tb <- renderPlot({
    req(df(), input$select_var1, input$time_var, input$id_var)
    if(input$time_type=="Discrete"){
      df_sum <- df()[, c(input$id_var, input$time_var, input$select_var1)] %>%
        rename(time=input$time_var, id=input$id_var, var=input$select_var1) %>%
        group_by(time) %>%
        mutate(med=median(var, na.rm = T)) %>% 
        mutate(Nmiss = sum(is.na(var)), 
               Pctmiss = sum(is.na(var))/length(var))
      xlab <- df_sum %>% select(time, Nmiss, Pctmiss) %>% distinct(.) %>%
        mutate(Nmiss = paste0(Nmiss, " (", round(100*Pctmiss, 2), "%)"))
      # summary plot
      plot_sum <- df_sum %>% 
        ggplot()+
        geom_boxplot(aes(x=time, y=var, group=time), outlier.size = 0.5, fill = "grey")+
        geom_jitter(aes(x=time, y=var, group=time), size = 0.5)+
        geom_line(data = df_sum %>% filter(!is.na(med)), aes(x=time, y=med))+
        scale_x_continuous(breaks = xlab$time, name = input$time_var,
                           sec.axis = sec_axis(~., name = "N (pct) of missing", breaks = xlab$time, label=xlab$Nmiss))+
        theme(axis.text.x.top = element_text(angle=90))+
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
  
  ## subtab 2.2: pairwise correlaion
  output$varnames2 <- renderUI({
    req(df())
    checkboxGroupInput("select_var2", label = "Variables", choices = colnames(df()), 
                       selected = colnames(df())[3:4])
  })
  
  output$trendp <- renderPlot({
    req(df(), input$select_var2, input$time_var, input$id_var, input$time_type)
    validate(need(length(input$select_var2)==2, "Please select a pair of variables."))
    
    df_pair <- df()[, c(input$time_var, input$id_var, input$select_var2)] %>%
      rename(time=input$time_var, id = input$id_var)
    t_uniq <- unique(df_pair$time)
    
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
        # correlation trend
        p2 <- df_pair %>% group_by(time) %>%
          group_modify(~{data.frame(cor = cor(.x[, input$select_var2], method = input$cor_type,
                                              use = "pairwise.complete.obs")[1, 2])})  %>%
          ungroup() %>% filter(complete.cases(.)) %>% ggplot()+
          geom_point(aes(x=time, y=cor))+
          geom_line(aes(x=time, y=cor))+
          labs(title = "Empirical correlation", x = input$time_var, y = " ")+
          scale_x_continuous(breaks = t_uniq)}
    else {
      p1 <- df_pair %>%
        pivot_longer(input$select_var2) %>%
        ggplot(aes(x = time, y=value, colour = name, group=interaction(id, name)))+
        geom_line(alpha = 0.7, linewidth=0.5)+
        scale_fill_brewer(palette = "Set2")+
        scale_color_brewer(palette = "Set2")+
        # scale_x_continuous(breaks = t_brk)+
        labs(title = "Variable trend", x = input$time_var, y = " ")+
        theme(legend.position = "bottom")
      # correlation trend
      p2 <- df_pair %>% group_by(time) %>%
        group_modify(~{data.frame(cor = cor(.x[, input$select_var2], method = input$cor_type,
                                            use = "pairwise.complete.obs")[1, 2])})  %>%
        ungroup() %>% filter(complete.cases(.)) %>% ggplot()+
        geom_point(aes(x=time, y=cor))+
        geom_line(aes(x=time, y=cor))+
        labs(title = "Empirical correlation", x = input$time_var, y = " ")
        # scale_x_continuous(breaks = t_brk)}
    }
    pall <- grid.arrange(p2, p1, ncol = 1, heights = c(0.7, 1))
    pall
  },height = 600, width = "auto")
  
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
    req(df(), input$time_var, input$id_var)
    ## correlation matrix
    cormat <- cor(subset(df() %>% rename(time=input$time_var, id=input$id_var) %>%
                           filter(time==input$time_bar1), 
                         select = -c(id, time)), method = input$cor_type, use = "pairwise.complete.obs")
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
    checkboxGroupInput("select_var3", label = "Variables", 
                       choices = colnames(df() %>% select(!c(input$time_var, input$id_var))))
  })
  
  ## data set to analyze
  confirmed <- reactiveVal(NULL)
  observeEvent(input$confirm, {confirmed(input$select_var3)})
  df_net <- reactive({
    req(df(), input$id_var, input$time_var)
    if(is.null(confirmed())){
      df_net <- df() %>% select(!c(input$id_var)) %>% 
        rename(time = input$time_var) %>% 
        filter(!if_all(!time, is.na))
      exclude_time <- df_net %>% group_by(time) %>% summarize(nrow=n()) %>% filter(nrow<=2)
      df_net %>% filter(!time %in% exclude_time$time)
    }
    else{
      df_net <- df()[, c(input$time_var, confirmed())] %>%
        rename(time = input$time_var) %>%
        filter(!if_all(confirmed(), is.na))
      exclude_time <- df_net %>% group_by(time) %>% summarize(nrow=n()) %>% filter(nrow<=2)
      df_net %>% filter(!time %in% exclude_time$time)
    }
  })
  ## time axis
  output$time_bar <- renderUI({
    req(df_net(), input$time_type)
    # time bar: by the original time
    tvec <- sort(unique(df_net()[, "time"]))
   if(input$time_type=="Discrete"){
      time_bar <- sliderTextInput("time_bar", label = input$time_var, choices = tvec, selected = tvec[1],
                      grid = TRUE)
    }
    else{
      time_bar <- sliderTextInput("time_bar", label = input$time_var, choices = tvec, selected = tvec[1],
                                  grid = FALSE)
    }
    time_bar
  }) # what if the time in the data set is not index but actual time (say, 0 to 1)?
  
  
  ## grouping results
  output$nclust <- renderUI({
    req(input$hclust)
    numericInput("nclust", label = "Number of groups", value = 1)
  })
  output$group_type <- renderUI({
    req(input$hclust)
    radioButtons("group_plot", label = "Group summary: ", choices = list("By temporal flow" = 1, 
                                                                         "By variable" = 2,
                                                                         "Hierarchical tree" = 3))
  })

  # calculate adjacency matrix at each time point
  adj_mat <- reactive({
    req(df_net(), input$time_type)
    mds_type <- ifelse(input$time_type=="Discrete", "Dynamic", "Splines")
    # if(is.null(confirmed())){
      GetAdjMat(data= df_net(),
                cor_method = input$cor_type,
                mds_type = mds_type)
  })
  

  # calculate graph
  graph_list <- reactive({
    req(adj_mat())
    graph_dyn_net(adj_mat(), cor_th = input$thres_cor)
  })

  # calculated hierarchical groups
  group_list <- reactive({
    req(adj_mat())
   lapply(adj_mat(),
          function(x){
            dis_mat <- 1-x
            # remove fully unmeasured variable
            var_id <- which(!is.na(diag(dis_mat)))
            dis_mat <- dis_mat[var_id, var_id]
            # remove pairs of variables whose dissimilarity cannot be computed
            # because there is <2 complete pairs
            dis_mat <- dis_mat[complete.cases(dis_mat), complete.cases(dis_mat)]
            dis_mat <- as.dist(dis_mat)
            hclust_fit <- hclust(dis_mat)
            return(hclust_fit)
                          })
  })
  
  ## calculate coordinates
  coord_list <- reactive({
      req(adj_mat(), df_net(), input$time_type)
      mds_type <- ifelse(input$time_type=="Discrete", "Dynamic", "Splines")
      t_uniq <- sort(unique(df_net()[, "time"])) # original time scale
      t_id <- seq_along(t_uniq) # time index
      # asynchronous
      # future({
        if(mds_type=="Splines"){
          SplinesMDS(adj_mat(), lambda = 10, P = dim(adj_mat()[[1]])[1],
                     tvec = t_uniq)
        }
        else{
          DynamicMDS(adj_mat(), 5)
        }
      # })
  })

  # plot
  output$mds_note1 <- renderPrint({
    req(input$time_type)
    mds_type <- ifelse(input$time_type=="Discrete", "Dynamic", "Splines")
    HTML(paste0("This networkplot is generated by ", mds_type, " multidimentional scaling based on ", input$time_type, " time and ", input$cor_type, " correlation. 
            To change the type of time variable or correlation, please move back to the previous tab.")) 
  })
  output$mds_note2 <- renderPrint({
    HTML("
    Observations may be excluded for the following reasons:
        <li>Empty rows or columns</li>
        <li>Less than two observations at any given time. Correlation is not realiable due to insufficiant sample. </li>
         ")
  })
    
  output$netp <- renderPlot({
   
      withProgress(
        value=0, message = "Processing", detail="This may take a while...",
        {
          req(adj_mat(), df_net(),  graph_list(), group_list(), coord_list())
          t_uniq <- sort(unique(df_net()[, "time"])) # original time scale
          t_id <- seq_along(t_uniq) # time index
          
          input_tid <- which(t_uniq==input$time_bar)
          graph_t <- graph_list()[[input_tid]]
          group_t <- group_list()[[input_tid]]
          coord_t <- coord_list()[[input_tid]]
          
          # color 
          if(input$hclust){
            group_t <- cutree(group_t, k = input$nclust)
            group_c <- brewer.pal(input$nclust, "Accent")[group_t]
            names(group_c) <- names(group_t)
            V(graph_t)$color <- group_c[V(graph_t)$name]
          }
          
    
          # plot
          # incProgress(0.5, detail = "Plotting")
          plot(graph_t,
           layout = as.matrix(coord_t),
           vertex.frame.color=ifelse(is.na(V(graph_t)$color), "grey", NA),
           vertex.label.cex=1,
           vertex.size = 20, 
           vertex.color = V(graph_t)$color, 
           margin = 0)
        
        cond <- sapply(adj_mat(), function(x){all(round(abs(x), 10)==1, na.rm = T) & !all(is.na(x))})
        if(any(cond)){warning(paste("Data showed perfect similarity at time", which(cond)))}
        
    incProgress(1)}
    )
    })
    ## group label plot
  output$group_plot <- renderPlot({
    req(input$hclust, input$nclust, group_list(), df_net(), input$group_plot, input$time_bar)
    tvec <- sort(unique((df_net()[ , "time"])))
    tid <- which(tvec==input$time_bar)

    group_label_list <- lapply(group_list(), function(hclust_fit){cutree(hclust_fit, k = input$nclust)})

    if(input$group_plot == 1){
        bind_rows(group_label_list, .id = "time") %>%
          mutate(time=tvec) %>%
          pivot_longer(-time) %>%
          mutate(value = as.factor(value), time = as.numeric(time)) %>%
          ggplot(aes(x=time, stratum = value, fill=value, color=value, alluvium=name))+
          geom_flow()+
          geom_stratum()+
          geom_vline(xintercept = input$time_bar)+
          labs(x="Time", y = "Group", title = "Group flow chart")+
          guides(fill=guide_legend("Group"), color=guide_legend("Group"))+
          theme(legend.position = "bottom", axis.text.y = element_blank())+
          scale_color_brewer(palette = "Accent")+
          scale_fill_brewer(palette = "Accent")+
          scale_x_continuous(breaks = tvec)
    } else if(input$group_plot == 2){
      bind_rows(group_label_list, .id = "time") %>%
        mutate(time=tvec) %>%
        pivot_longer(-time) %>%
        mutate(value = as.factor(value), time = as.numeric(time)) %>%
        ggplot(aes(x=time, y=name, fill=value))+
        geom_tile(alpha=0.7)+
        geom_vline(xintercept = input$time_bar)+
        labs(x="Time", y = " ", fill = "Group", title = "Variable group assignment")+
        theme(legend.position = "bottom")+
        scale_fill_brewer(palette = "Accent")+
        scale_x_continuous(breaks = tvec)
    } else{
       ggdendrogram(group_list()[[tid]], rotate = T, size = 2)+
            labs(title = paste0("Hierarchial group at time ", input$time_bar))
      
    }
  })
  ## note
  output$group_note <- renderText({
    req(input$hclust)
    if(input$group_plot==1){
      HTML("
        <li>This plot visualizes the change of group structure of variables over time using an <a href='https://corybrunson.github.io/ggalluvial/' target='_blank'>Alluvial plot</a></li>
        <li>Band with different colors represents different groups, and the width of band represents the size of group.</li>
        <li>Band flow between different groups across time represents variables that moved from one group to the other, and the width of flow represents number of variables that made the switch.</li>
         ")
    } else if(input$group_plot==2){
      HTML("
        <li>This plot visualizes the change of group assignment for each variable.</li>
        <li>Each row represents a single variable, and the color represents the group it is assigned to at specific time points.</li>
        <li>Change of color indicates change of group assignments.</li>
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
    checkboxGroupInput("select_var4", label = "Variables", 
                       choices = colnames(df() %>% select(!c(input$time_var, input$id_var))))
  })
  confirmed2 <- reactiveVal(NULL)
  observeEvent(input$confirm2, {confirmed2(input$select_var4)})
 ## clustering
  output$nclust2 <- renderUI({
    req(input$hclust2)
    numericInput("nclust2", label = "Number of groups", value = 1)
  })
  # integral adjacency matrix
  int_adj_mat <- reactive({
    req(df(), input$id_var, input$time_var)
    if(is.null(confirmed2())){
      int_adj_mat_list <- GetAdjMat(data= df() %>% select(!c(input$id_var)) %>% rename(time = input$time_var),
                cor_method = input$cor_type,
                mds_type = "Splines")
    }
    else{
      int_adj_mat_list <- GetAdjMat(data=df()[, c(input$time_var, confirmed2())] %>% rename(time = input$time_var),
                cor_method = input$cor_type,
                mds_type = "Splines")
    }
    # integrate
    AveAdj <- apply(simplify2array(int_adj_mat_list), c(1, 2), mean, na.rm = T)
    AveAdj
  })
  ## hierarchical clustering result
  int_hclust <- reactive({
    req(int_adj_mat())
    validate(need(sum(is.na(int_adj_mat())) == 0, "Correlation cannot be calculated, likely due to empty or uniform columns."))
    AveDis <- 1-int_adj_mat()
    hclust(dist(AveDis))
  })
  ## coordinates
  int_coords <- reactive({
    req(int_adj_mat())
    mds(int_adj_mat())$conf
  })
  ## network plot
  output$int_net <- renderPlot({
    req(int_adj_mat(), int_hclust(), input$thres_cor2, int_coords())
    # checks
    validate(need(sum(is.na(int_adj_mat())) == 0, "Correlation cannot be calculated, likely due to empty or uniform columns."))
    # initialize graph
    int_adj <- int_adj_mat()
    int_adj[which(int_adj < input$thres_cor2)] <- 0
    int_net <- graph_from_adjacency_matrix(int_adj,
                                         mode = "undirected", weighted = T, diag=F)

    # edges
    E(int_net)$width <- E(int_net)$weight*5

    # node properties
    V(int_net)$color <- rgb(0.2, 0.4, 0.8, alpha=0.4)
    coords <- mds(int_adj_mat())$conf

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
    req(input$hclust2)
    HTML(paste0("Both plots are generated based on the integerated ", input$cor_type, " correlation matrix. To change the type of correlation, please move back to the second tab."))
  })
  
}



# Run the application 
shinyApp(ui = ui, server = server)
