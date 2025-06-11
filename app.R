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
ui <- page_fluid(
  
  # side bar input
  # maybe change to multipage layout page_nabvar later
  titlePanel("Temporal network visualization of multidimensional data"),
  sidebarLayout(
    sidebarPanel(
      # upload file
      fileInput(inputId = "df_path", label = "Upload data",
                accept = c(".csv")), # Maybe add more data types? 
      # checkboxInput(inputId = "vars", label="Variables", value = T),
      
      # choose correlation threshold
      sliderInput(inputId = "thres_cor", label = "Visualization threshold", min=0, max=1, value=0, step=0.01,
                  ticks=FALSE),
      
      # time bar
      uiOutput("time_bar"),
      
      # variable list 
      uiOutput("varnames")
      
    ),
  
  
  # display
  mainPanel(
    navset_card_underline(
      nav_panel("Data", dataTableOutput("show_df")),
      nav_panel("Network", plotOutput("netp")),
      nav_panel("Temporal trend", plotOutput("trendp"))
    )
  )
  
  )
  
)


#### Server ####

# Define server logic required to draw a histogram
server <- function(input, output) {
  # display loaded data
  df <- reactive({
    req(input$df_path)
    read.csv(input$df_path$datapath)
  })
  
  output$show_df <- renderDataTable({
    df()
  })
  
  # time axis
  output$time_bar <- renderUI({
    req(df())
    # time bar
    Tmax <- length(unique(df()$time))
    sliderInput("time_bar", "Time (index)", min=1, max = Tmax, value=1, step = 1, ticks = F)
  }) # what if the time in the data set is not index but actual time (say, 0 to 1)?
  
  # variable list
  output$varnames <- renderUI({
    req(df())
    checkboxGroupInput("select_var", label = "Variables", choices = colnames(df()))
  })
  
  # calculate adjacency matrix at each time point
  adj_mat <- reactive({
    req(df())
    GetAdjMat(data=df())
  })
  
  # calculate graph
  graph_list <- reactive({
    req(adj_mat())
    graph_dyn_net(adj_mat(), cor_th = input$thres_cor)
  }) 
  
  # calculate coordinates
  coord_list <- reactive({
    req(adj_mat())
    DynamicMDS(adj_mat(), 10)
  })
  
  output$netp <- renderPlot(plot(graph_list()[[input$time_bar]], 
                                 layout = as.matrix(coord_list()[[input$time_bar]])))
  
  output$trendp <- renderPlot({
    req(df())
    df_pair <- df()[, c("id", "time", input$select_var)]
    # individual trend
    p1 <- df_pair %>%
      pivot_longer(input$select_var) %>%
      ggplot(aes(x=time, y=value, group=time))+
      geom_boxplot()+geom_jitter(size=0.5)+
      facet_wrap(~name, ncol=1)
    # correlation trend
    p2 <- df_pair %>% group_by(time) %>%
      group_modify(~{data.frame(cor = cor(.x[, input$select_var], method = "pearson", 
                                          use = "pairwise.complete.obs")[1, 2])})  %>%
      ungroup() %>% ggplot(aes(x=time, y=cor))+
      geom_point()+
      geom_line()
   pall <- grid.arrange(p1, p2, ncol = 1, heights = c(2, 1))
   pall
  })

}

# Run the application 
shinyApp(ui = ui, server = server)
