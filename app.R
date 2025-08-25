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
library(smacof)
library(splines2)
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
               fileInput(inputId = "df_path", label = "Upload data",
                         accept = c(".csv")),
               # specify subject ID and time
               uiOutput("time_var"),
               # uiOutput("bin_yn"),
               # uiOutput("bin_width"),
               uiOutput("id_var")
             ),
               
             
             # main panel: data preview
             mainPanel(h3('Data preview'),
                       dataTableOutput("show_df"),
                       br(), 
                       htmlOutput("size_info"))
    
           )
           ),
  
  # tab 2: data summary table
  tabPanel(title = "Data summary",
           # side bar 1
           sidebarLayout(
             sidebarPanel(
             uiOutput("varnames1"),
             uiOutput("time_type") # how to approach time
             ),

           # main panel
           mainPanel(# summary for selected variable
                     h3('Single variable summary'),
                     plotOutput("sum_tb")
                     )
  )),
  
  # tab 3: network and subnetwork plots
  tabPanel(title = "Dynamic network",
           sidebarLayout(
             # side bar
             sidebarPanel(
               # visualization type
               selectInput("mds_type", label = "Method of visualization",
                           choices = list("Splines", "Dynamic"),
                           selected = "Dynamic"),
               # correlation type
               selectInput("cor_type", label="Type of correlation",
                           choices = list("pearson", "spearman")),
               # choose correlation threshold
               sliderInput(inputId = "thres_cor", label = "Visualization threshold", min=0, max=1, value=1, step=0.01,
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
             mainPanel(plotOutput("netp", width = "100%", height = "600px")
                       # htmlOutput("mesg")
                       )
             )),
  
  # tab 4: pairwise correlation over time
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
  
  # tab 5: correlation heatmap at one time slice
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
  options(shiny.maxRequestSize=10*1024^2)
  # tab 1
  ## data upload
  df <- reactive({
    req(input$df_path)
    read.csv(input$df_path$datapath)
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
  
  # tab 2
  ## variable list
  output$varnames1 <- renderUI({
    req(df())
    selectInput("select_var1", label = "Variables", choices = colnames(df()), 
                       selected = colnames(df())[5])
  })
  # type of time
  output$time_type <- renderUI({
    req(df())
    selectInput("time_type", "Treat time as", choices = c("Continuous", "Discrete"), selected = "Discrete")
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
        geom_line(aes(x=time, y=med))+
        scale_x_continuous(breaks = xlab$time, name = input$time_var,
                           sec.axis = sec_axis(~., name = "N (pct) of missing", breaks = xlab$time, label=xlab$Nmiss))+
        theme(axis.text.x.top = element_text(angle=45))+
        labs(x=input$time_var, y=input$select_var1, title = "Temporal trend")
    }
    else{
      plot_sum <- df()[, c(input$id_var, input$time_var, input$select_var1)] %>%
        rename(time=input$time_var, id=input$id_var, var=input$select_var1) %>%
        ggplot()+
        geom_line(aes(x=time, y=var, group=id), alpha = 0.5, linewidth = 0.5)+
        geom_smooth(aes(x=time, y=var), method = gam, formula = y~s(x))
    }
    
    plot_sum
  })
  
  # tab 3
  ## time axis
  output$time_bar <- renderUI({
    req(df(), input$time_var, input$time_type)
    # time bar: by the original time 
    if(input$time_type=="Discrete"){
      tvec <- sort(unique(df()[, input$time_var]))
      time_bar <- sliderTextInput("time_bar", label = input$time_var, choices = tvec, selected = tvec[1],
                      grid = TRUE)
    }
    else{
      trange <- range(df()[ ,input$time_var], na.rm=T)
      time_bar <- sliderInput("time_bar", label = input$time_var, min=trange[1], max=trange[2], value=trange[1],
                              ticks=FALSE)
    }
    time_bar
  }) # what if the time in the data set is not index but actual time (say, 0 to 1)?
  
  output$varnames2 <- renderUI({
    req(df(), input$time_var, input$id_var)
    checkboxGroupInput("select_var2", label = "Variables", 
                       choices = colnames(df() %>% select(!c(input$time_var, input$id_var))))
  })

  # calculate adjacency matrix at each time point
  confirmed <- reactiveVal(NULL)
  observeEvent(input$confirm, {confirmed(input$select_var2)})
  adj_mat <- reactive({
    req(df())
    if(is.null(confirmed())){
      GetAdjMat(data= df() %>% select(!c(input$id_var)) %>% rename(time = input$time_var), 
                cor_method = input$cor_type,
                mds_type = input$mds_type)
    }
    else{
      GetAdjMat(data=df()[, c(input$time_var, confirmed())] %>% rename(time = input$time_var), 
                cor_method = input$cor_type,
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
    t_uniq <- sort(unique(df()[, input$time_var])) # original time scale
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
    req(adj_mat(), input$nclust)
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
            return(cutree(hclust_fit, k = input$nclust))
                          })
  })

  # plot
  output$netp <- renderPlot({
    # plot
    withProgress(
       value=0, message = "Processing", detail="This may take a while...",
        {
          req(input$time_bar, df(), graph_list(), group_list(), coord_list())
          # find the location index
          tvec <- sort(unique(df()[, input$time_var]))
          input_tid <- which(tvec==input$time_bar)
          
          # graph at this time point
          graph_t <- graph_list()[[input_tid]]
          group_t <- group_list()[[input_tid]]
          coord_t <- coord_list()[[input_tid]]
          
          if(input$hclust){
            V(graph_t)$color <- group_t[V(graph_t)$name]
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
    }, height = 600, width = "auto")
  
  ## messages
  # output$mesg <- ({
  #   req(adj_mat(), input$time_bar)
  #   adj_t <- adj_mat()[[input$time_bar]]
  #   
  #   # fully missing variables
  #   full_mis_var <- which(is.na(diag(adj_t)))
  #   
  #   # partiallly missing 
  #   adj_t <- adj_t[-full_mis_var, -full_mis_var]
  #   adj_t_long <- data.frame(adj_t) %>%
  #     rownames_to_column() %>%
  #     pivot_longer(!rowname) 
  # })
  
  
  # tab 4
  ## variable list
  output$varnames3 <- renderUI({
    req(df())
    checkboxGroupInput("select_var3", label = "Variables", choices = colnames(df()), 
                       selected = colnames(df())[3:4])
  })
  
  ## plot
  output$trendp <- renderPlot({
    req(df(), input$select_var3, input$time_var, input$id_var)
    df_pair <- df()[, c(input$time_var, input$id_var, input$select_var3)] %>%
      rename(time=input$time_var, id = input$id_var)
    t_uniq <- unique(df_pair$time)
    # individual trend
    p1 <- df_pair %>%
      pivot_longer(input$select_var3) %>%
      ggplot(aes(x=time, y=value, group=time))+
      geom_boxplot()+geom_jitter(size=0.5)+
      facet_wrap(~name, ncol=1)+
      labs(x=input$time_var, y=" ")+
      scale_x_continuous(breaks = t_uniq)
    # correlation trend
    p2 <- df_pair %>% group_by(time) %>%
      group_modify(~{data.frame(cor = cor(.x[, input$select_var3], method = input$cor_type,
                                          use = "pairwise.complete.obs")[1, 2])})  %>%
      ungroup() %>% filter(complete.cases(.)) %>% ggplot()+
      geom_point(aes(x=time, y=cor))+
      geom_line(aes(x=time, y=cor))+
      labs(title = "Empirical correlation", x = input$time_var, y = " ")+
      scale_x_continuous(breaks = t_uniq)
   pall <- grid.arrange(p1, p2, ncol = 1, heights = c(2, 1))
   pall
  },height = 600, width = "auto")
  
  # tab 4
  output$time_bar4 <- renderUI({
    req(df(), input$time_var, input$id_var)
    # time bar
    tvec <- sort(unique((df()[ ,input$time_var])))
    sliderTextInput("time_bar4", "Time", choices = tvec, selected = tvec[1],
                    grid = TRUE)
  })
  output$heatmap <- renderPlot({
    cormat <- cor(subset(df() %>% rename(time=input$time_var, id=input$id_var) %>%
                           filter(time==input$time_bar4), 
                         select = -c(id, time)), method = input$cor_type, use = "pairwise.complete.obs")
    col_id <- colnames(cormat)[!is.na(diag(cormat))] 
    cormat <- cormat[col_id, col_id]
    heatmap(cormat, distfun = function(mat){as.dist(1-abs(mat))}, margins = c(10, 10), keep.dendro = FALSE)
  }, height = 600, width = "auto")
   
  

}



# Run the application 
shinyApp(ui = ui, server = server)
