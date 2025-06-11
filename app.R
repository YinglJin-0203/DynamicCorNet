# This shiny app is built for network visualization 
# of temporal multivariate data

rm(list=ls())
library(shiny)
library(bslib)
library(here)
library(DT)


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
      uiOutput("time_bar")
      
    ),
  
  
  # display
  mainPanel(
    navset_card_underline(
      nav_panel("Data", dataTableOutput("show_df")),
      nav_panel("Network", plotOutput("netp"))
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
  output$time_bar <- renderUI({
    req(df())
    # trange <- range(df()$time)
    Tmax <- length(unique(df()$time))
    sliderInput("time_bar", "Time (index)", min=1, max = Tmax, value=1, step = 1, ticks = F)
  }) # what if the time in the data set is not index but actual time (say, 0 to 1)?
  
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

}

# Run the application 
shinyApp(ui = ui, server = server)
