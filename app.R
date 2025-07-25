# This shiny app is built for network visualization 
# of temporal multivariate data

# rm(list=ls())
library(shiny)
library(shinyWidgets)
library(bslib)
library(here)
library(DT)
library(tidyverse)
library(gridExtra)
library(arsenal)
library(htmltools)
theme_set(theme_minimal())


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

  # tab 1: data summary table
  tabPanel(title = "Data summary", 
           # side bar 1
           sidebarLayout(
             sidebarPanel(
             # upload file
             fileInput(inputId = "df_path", label = "Upload data",
                       accept = c(".csv")),
             uiOutput("varnames1")),
           
           # main panel
           mainPanel(h3('Data overview'),
                     dataTableOutput("show_df"),
                     # also add summary for selected variable
                     h3('Single variable summary'),
                     dataTableOutput("sum_tb")
                     # plotOutput("miss_plot")
                     )
  )),
  
  # tab 2: network and subnetwork plots
  tabPanel(title = "Dynamic network",
           sidebarLayout(
             # side bar
             sidebarPanel(
               # visualization type
               selectInput("mds_type", label = "Method of visualization",
                           choices = list("Splines", "Dynamic"),
                           selected = "Splines"),
               # correlation type
               selectInput("cor_type", label="Type of correlation",
                           choices = list("pearson", "spearman")),
               # choose correlation threshold
               sliderInput(inputId = "thres_cor", label = "Visualization threshold", min=0, max=1, value=0, step=0.01,
                           ticks=FALSE),
               # time bar
               uiOutput("time_bar"),
               # hierarchical grouping
               checkboxInput("hclust", label = tags$span("Show groupng results", 
                                                          style = 'font_weight: bold; font-size: 18px;'), 
                             value = TRUE),
               numericInput("nclust", label = "Number of groups", value = 3),
               # variable list
               uiOutput("varnames2"),
               actionButton("confirm", "Confirm selection")
             ),
             
             # main panel
             mainPanel(plotOutput("netp", width = "100%", height = "600px"))
             )),
  
  # tab 3: pairwise correlation over time
  tabPanel(title = "Pairwise correlation overtime",
           sidebarLayout(
             # side bar
             sidebarPanel(uiOutput("varnames3")),
             
             # main panel
             mainPanel(
               plotOutput("trendp"),
               tableOutput("sum_tb_temp1"),
               tableOutput("sum_tb_temp2")
                 )
             )
               
           ),
  
  # tab 4: correlation heatmap at one time slice
  tabPanel(title = "Correlation structure at static time points",
           sidebarLayout(
             # side bar
             sidebarPanel(uiOutput("time_bar4")),
             
             # main panel
             mainPanel(
               plotOutput("heatmap")
           )))
)
  
  
#### Server ####

server <- function(input, output) {
  # tab 1
  ## data preview
  df <- reactive({
    req(input$df_path)
    read.csv(input$df_path$datapath)
  })
  
  output$show_df <- renderDataTable({df()},
    options = list(scrollX = T, fixedHeader=T)
  )
  
  ## variable list
  output$varnames1 <- renderUI({
    req(df())
    selectInput("select_var1", label = "Variables", choices = colnames(df()), 
                       selected = colnames(df())[5])
  })
  
  ## summary of variables (across all time)
  
  output$sum_tb <- renderDataTable({
    req(df(), input$select_var1)
    df_sum <- df()[, c("time", input$select_var1)]
    # summary table
    tb_sum <- df_sum %>% 
      group_by(time) %>%
      summarise(Min  = min(.data[[input$select_var1]], na.rm = TRUE),
                Mean = mean(.data[[input$select_var1]], na.rm = TRUE),
                Median = median(.data[[input$select_var1]], na.rm = TRUE),
                Max  = max(.data[[input$select_var1]], na.rm = TRUE),
                SD   = sd(.data[[input$select_var1]], na.rm = TRUE),
                Nmiss = sum(is.na(.data[[input$select_var1]]) | is.infinite(.data[[input$select_var1]])),
                N = length(.data[[input$select_var1]]), 
                .groups = "drop") %>%
      mutate(Npct = 100*Nmiss/N) %>%
      mutate(Nmiss = paste0(Nmiss, " (", round(Npct, 2), "%)")) %>%
      select(-N, -Npct) %>%
      mutate_at(vars(-time, -Nmiss), round, 2)
    # caption
    cap_text <- paste0("Variable: ", input$select_var1)
    datatable(tb_sum, 
              caption=tags$caption(style = 'font-weight: bold; font-size: 15px; color: #2C3E50;',
                                   cap_text),
              options = list(dom="t",
                             columnDefs = list(
                               list(targets = 0, className = 'dt-left'),   # First column
                               list(targets = 1:(ncol(tb_sum) - 1), className = 'dt-center')  # All others
                             ))
              )
    
  })

  
  # tab 2
  ## time axis
  output$time_bar <- renderUI({
    req(df())
    # time bar: by the original time 
    tvec <- sort(unique(df()$time))
    sliderTextInput("time_bar", "Time", choices = tvec, selected = tvec[1],
                    grid = TRUE)
  }) # what if the time in the data set is not index but actual time (say, 0 to 1)?
  
  output$varnames2 <- renderUI({
    req(df())
    checkboxGroupInput("select_var2", label = "Variables", 
                       choices = colnames(subset(df(), select = -c(id, time))))
  })

  # calculate adjacency matrix at each time point
  confirmed <- reactiveVal(NULL)
  observeEvent(input$confirm, {confirmed(input$select_var2)})
  adj_mat <- reactive({
    req(df())
    if(is.null(confirmed())){
      GetAdjMat(data=subset(df(), select = -c(id)), 
                cor_method = input$cor_type,
                mds_type = input$mds_type)
    }
    else{
      GetAdjMat(data=df()[, c("time", confirmed())], cor_method = input$cor_type,
                mds_type = input$mds_type)
    }
  })

  # calculate graph
  graph_list <- reactive({
    req(adj_mat())
    graph_dyn_net(adj_mat(), cor_th = input$thres_cor)
  })

  # calculate coordinates
  coord_list <- reactive({
    req(adj_mat(), df())
    t_uniq <- sort(unique(df()$time)) # original time scale
    t_id <- seq_along(t_uniq) # time index
    
    if(input$mds_type=="Splines"){
      SplinesMDS(adj_mat(), lambda = 10, K = 20, P = dim(adj_mat()[[1]])[1], 
                 tvec = t_uniq)
    }
    else{
      DynamicMDS(adj_mat(), 5)
    }
  })
  
  # calculated hierarchical groups
  group_list <- reactive({
    req(adj_mat(), input$hclust, input$nclust)
   lapply(adj_mat(),
          function(x){
            dis_mat <- 1-x
            var_id <- which(!is.na(diag(dis_mat)))
            dis_mat <- dis_mat[var_id, var_id]
            dis_mat <- as.dist(dis_mat)
            hclust_fit <- hclust(dis_mat)
            return(cutree(hclust_fit, k = input$nclust))
                          })
  })

  # plot
  output$netp <- renderPlot({
    # req(input$time_bar, df(), input$hclust, graph_list(), group_list(), coord_list())
    # # find the location index
    # tvec <- sort(unique(df()$time))
    # input_tid <- which(tvec==input$time_bar)
    # 
    # # graph at this time point
    # graph_t <- graph_list()[[input_tid]]
    # color_t <- group_list()[[input_tid]]
    # coord_t <- coord_list()[[input_tid]]
    # 
    # if(input$hclust){
    #   V(graph_t)$color <- color_t[V(graph_t)]
    # }
    
    # plot
    withProgress(
       value=0, message = "Processing", detail="This may take a while...",
        {
          req(input$time_bar, df(), input$hclust, graph_list(), group_list(), coord_list())
          # find the location index
          tvec <- sort(unique(df()$time))
          input_tid <- which(tvec==input$time_bar)
          
          # graph at this time point
          graph_t <- graph_list()[[input_tid]]
          color_t <- group_list()[[input_tid]]
          coord_t <- coord_list()[[input_tid]]
  
          if(input$hclust){
            V(graph_t)$color <- color_t[V(graph_t)]
          }

          # plot
          # incProgress(0.5, detail = "Plotting")
          plot(graph_t,
           layout = as.matrix(coord_t),
           vertex.frame.color=V(graph_t)$color,
           vertex.label.cex=1,
           vertex.size = 20, 
           vertex.color = V(graph_t)$color, 
           margin = 0)
    
    cond <- sapply(adj_mat(), function(x){all(round(abs(x), 10)==1, na.rm = T) & !all(is.na(x))})
    if(any(cond)){warning(paste("Data showed perfect similarity at time", which(cond)))}
    
    incProgress(1)}
    )
    }, height = 600, width = "auto")
  
  # output$try <- renderPrint({
  #   
  #   lapply(1:length(graph_list()),
  #          function(x){
  #            print(V(graph_list()[[x]])$color)
  #          })
  # })
  
  # tab3
  ## variable list
  output$varnames3 <- renderUI({
    req(df())
    checkboxGroupInput("select_var3", label = "Variables", choices = colnames(df()), 
                       selected = colnames(df())[3:4])
  })
  
  ## plot
  output$trendp <- renderPlot({
    req(df(), input$select_var3)
    df_pair <- df()[, c("id", "time", input$select_var3)]
    # individual trend
    p1 <- df_pair %>%
      pivot_longer(input$select_var3) %>%
      ggplot(aes(x=time, y=value, group=time))+
      geom_boxplot()+geom_jitter(size=0.5)+
      facet_wrap(~name, ncol=1)+
      labs(x="Time", y=" ")+
      scale_x_continuous(breaks = 1: length(unique(df()$time)))
    # correlation trend
    p2 <- df_pair %>% group_by(time) %>%
      group_modify(~{data.frame(cor = cor(.x[, input$select_var3], method = input$cor_type,
                                          use = "pairwise.complete.obs")[1, 2])})  %>%
      ungroup() %>% ggplot(aes(x=time, y=cor))+
      geom_point()+
      geom_line()+
      labs(title = "Empirical correlation", x = "Time", y = " ")+
      scale_x_continuous(breaks = 1: length(unique(df()$time)))
   pall <- grid.arrange(p1, p2, ncol = 1, heights = c(2, 1))
   pall
  })
  ## temporal summary statistics
  # output$sum_tb_temp1 <- renderDataTable({
  #  df()[, c("time", input$select_var3[1])] %>%
  #     group_by(time) %>%
  #     group_modify(~{data.frame(N = length(.x[, input$select_var3[1]]), 
  #                               Nmiss=sum(is.na(.x[, input$select_var3[1]])), 
  #                               # Min = min(.x, na.rm = T),
  #                               # Mean=mean(.x, na.rm=T), Median=median(.x, na.rm = T),
  #                               # Max = max(.x, na.rm =T), SD = sd(.x, na.rm = T)
  #                               )

  # })})
  
  # tab 4
  output$time_bar4 <- renderUI({
    req(df())
    # time bar
    tvec <- sort(unique((df()$time)))
    sliderTextInput("time_bar4", "Time", choices = tvec, selected = tvec[1],
                    grid = TRUE)
  })
  output$heatmap <- renderPlot({
    cormat <- cor(subset(df() %>% filter(time==input$time_bar4), 
                         select = -c(id, time)), method = input$cor_type, use = "pairwise.complete.obs")
    col_id <- colnames(cormat)[!is.na(diag(cormat))] 
    cormat <- cormat[col_id, col_id]
    heatmap(cormat)
  })
   
  

}



# Run the application 
shinyApp(ui = ui, server = server)
