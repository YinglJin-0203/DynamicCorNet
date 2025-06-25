# This shiny app is built for network visualization 
# of temporal multivariate data

rm(list=ls())
library(shiny)
library(bslib)
library(here)
library(DT)
library(tidyverse)
library(gridExtra)


#### Helper functions #### 

source(here("helpers/Stress.R"))
source(here("helpers/AdjacencyMat.R"))
source(here("helpers/DynNet.R"))
source(here("helpers/DynamicMDS.R"))


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
           mainPanel(dataTableOutput("show_df"),
                     # also add summary for selected variable
                     dataTableOutput("sum_tb"),
                     plotOutput("sum_plot")
                     )
  )),
  
  # tab 2: network and subnetwork plots
  tabPanel(title = "Dynamic network",
           sidebarLayout(
             # side bar
             sidebarPanel(
               # correlation type
               selectInput("cor_type", label="Type of correlation",
                           choices = list("pearson", "spearman")),
               # choose correlation threshold
               sliderInput(inputId = "thres_cor", label = "Visualization threshold", min=0, max=1, value=0, step=0.01,
                           ticks=FALSE),
               # time bar
               uiOutput("time_bar"),
               # variable list
               uiOutput("varnames2"),
               actionButton("confirm", "Confirm selection")
             ),
             
             # main panel
             mainPanel(plotOutput("netp"))
           )),
  
  # tab 3: pairwise correlation over time
  tabPanel(title = "Pairwise correlation overtime",
           sidebarLayout(
             # side bar
             sidebarPanel(uiOutput("varnames3")),
             
             # main panel
             mainPanel(
               plotOutput("trendp"),
               fluidRow(
                 column(6, tableOutput("sum_tb_temp1")),
                 column(6, tableOutput("sum_tb_temp1"))
                 ))
               
           )),
  
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
                       selected = colnames(df())[3])
  })
  
  ## summary of variables (across all time)
  
  # output$sum_tb <- renderDataTable({
  #   sum_var <- df()[, input$select_var1]
  #   if(is.numeric(sum_var)){
  #     data.frame(N = length(sum_var), Nmiss=sum(is.na(sum_var)), 
  #                Min = min(sum_var, na.rm = T), Mean=mean(sum_var, na.rm=T),
  #                Median=median(sum_var, na.rm = T), Max = max(sum_var, na.rm =T),
  #                SD = sd(sum_var, na.rm = T))
  #   }
  #   else{
  #     data.frame(N = length(sum_var), Nmiss=sum(is.na(sum_var)))
  #   }
  # }, options = list(dom="t"))
  # 
  
  # tab 2
  ## time axis
  output$time_bar <- renderUI({
    req(df())
    # time bar
    Tmax <- length(unique(df()$time))
    sliderInput("time_bar", "Time (index)", min=1, max = Tmax, value=1, step = 1, ticks = F)
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
      GetAdjMat(data=subset(df(), select = -c(id)), cor_method = input$cor_type)
    }
    else{
      GetAdjMat(data=df()[, c("time", confirmed())], cor_method = input$cor_type)
    }
  })

  # calculate graph
  graph_list <- reactive({
    req(adj_mat())
    graph_dyn_net(adj_mat(), cor_th = input$thres_cor)
  })

  # calculate coordinates
  coord_list <- reactive({
    req(adj_mat())
    DynamicMDS(adj_mat(), 5)
  })

  output$netp <- renderPlot(plot(graph_list()[[input$time_bar]],
                                 layout = as.matrix(coord_list()[[input$time_bar]]),
                                 vertex.frame.color=NA, vertex.label.cex=0.5,
                                 vertex.color = rgb(0.2, 0.4, 0.8, alpha=0.4)))
  # tab3
  ## variable list
  output$varnames3 <- renderUI({
    req(df())
    checkboxGroupInput("select_var3", label = "Variables", choices = colnames(df()), 
                       selected = colnames(df())[3:4])
  })
  ## plot
  output$trendp <- renderPlot({
    req(df())
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
  #   tb1 <- df()[, c("id", "time", input$select_var3[1])] %>%
  #     group_by(time) %>%
  #     group_modify(~{data.frame(N = length(.x), Nmiss=sum(is.na(.x)), Min = min(.x, na.rm = T),
  #                               Mean=mean(.x, na.rm=T), Median=median(.x, na.rm = T),
  #                               Max = max(.x, na.rm =T), SD = sd(.x, na.rm = T))
  # 
  # })})
  
  # tab 4
  output$time_bar4 <- renderUI({
    req(df())
    # time bar
    Tmax <- length(unique(df()$time))
    sliderInput("time_bar4", "Time (index)", min=1, max = Tmax, value=1, step = 1, ticks = F)
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
