# This shiny app is built for network visualization 
# of temporal multivariate data

library(shiny)
library(bslib)
library(DT)

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
                  ticks=FALSE)
      
    ),
  
  
  # display
  mainPanel(
    dataTableOutput("show_df")
  )
  
  )
  
)

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


}

# Run the application 
shinyApp(ui = ui, server = server)
