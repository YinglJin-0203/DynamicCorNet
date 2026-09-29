# This shiny app is built for network visualization 
# of temporal multivariate data

rm(list=ls())
library(shiny)
library(shinyWidgets)
library(bslib)
library(shinyBS)
library(bsicons)
library(here)
library(DT)
library(tidyverse)
library(gridExtra)
library(htmltools)
library(smacof)
library(splines2)
library(RColorBrewer)
library(mgcv)
library(igraph)
library(ggplot2)
library(rsconnect)


print(sessionInfo())
print(R.version.string)

theme_set(theme_minimal())

set.seed(62)


#### Helper functions #### 

source("Code/dyn_mds.R")
source("Code/lambda_sweep.R")
source("Code/dMDS_Helpers.R")
source("Code/get_similarity.R")
source("Code/lcurve_corner_dist.R")
source("Code/lcurve_corner_menger.R")



#### User interface ####

# UI includes the following elements
# visualization threshold of correlation
ui <- fluidPage(

  # ---- Responsive CSS ----
  # Fallback breakpoint rules for small screens, layered on top of the
  # column(width=, md=) breakpoints used throughout. This forces sidebars
  # and any remaining fixed-width elements to stack full-width below 768px,
  # and lets plot containers shrink instead of clipping/squishing.
  tags$head(
    tags$style(HTML("
      @media (max-width: 768px) {
        .well { width: 100% !important; }
        .col-sm-3 { width: 100% !important; }
        .col-sm-6 { width: 100% !important; }
        .col-sm-9 { width: 100% !important; }
      }
      /* plots and tables should never overflow their container */
      .shiny-plot-output, .dataTables_wrapper {
        max-width: 100%;
      }
    "))
  ),

  navbarPage(title = "Multivariate Longitudinal Exploratory Data Analysis",
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
               uiOutput("id_var"),
               uiOutput("time_var"),
               width = 3
             ),
               
             
             # main panel: data preview
             mainPanel(h3('Data preview'),
                       div(style = "overflow-x: auto;",
                           DTOutput("show_df")
                       ),
                       br(), 
                       h4("Sample summary"), 
                       div(style = "overflow-x: auto;",
                           dataTableOutput("size_info")
                       ),
                       width = 9
                       )
    
           )
           ),
  
  # tab 2: descriptive
  tabPanel(title = "Descriptives",
           tabsetPanel(
              ## subtab 2.1: single variable distribution
                tabPanel(title = "Univariate analysis",
                         # side bar 1
                         sidebarLayout(
                           sidebarPanel(
                           uiOutput("varnames1"),
                           # summarize type
                           selectInput("sum_type", "Summarize by", choices = c("Time", "Individual"), 
                                       selected = "Time", multiple = FALSE),
                           uiOutput("tab1_note"),
                           width = 3),
  
                         # main panel
                         mainPanel(
                                   fluidRow(
                                     column(width = 12, md = 6,
                                            h4("Summary plot"),
                                            plotOutput("sum_plt", height = "45vh")
                                            ),
                                     column(width = 12, md = 6,
                                            h4("Missing pattern"),
                                            plotOutput("miss_plot", height = "45vh")
                                            )
                                   ),
                                   h4("Summary statistics"),
                                   div(style = "overflow-x: auto;",
                                       dataTableOutput("sum_tb")
                                   ),
                                   width = 9
                         )
              )),
              
              ## subtab 2.2: pairwise
              tabPanel(title = "Bivariate analysis",
                       # side bar
                       sidebarLayout(
                        sidebarPanel( 
                          # summarize type
                          selectInput("sum_type2", "Summarize by", choices = c("Time", "Trajectory"), 
                                      selected = "Time", multiple = FALSE),
                          # correlation type
                          selectInput("cor_type", label="Type of correlation",
                                      choices = list("pearson", "spearman")),
                          tagList(
                            icon("info-circle"),
                            em("Correlation measures may be unrealiable when the proportion of missing is large! 
                            When the number of complete pairs is less then 10, the correlation is removed from visualization. 
                            Details of these time points can be examined in the previous subtab.")
                            ),
                          # variable list
                          br(),br(),
                          # scale
                          checkboxInput("scaleY", "Scale correlation axis to data?", value = F),
                          tagList(
                            icon("info-circle"),
                            em("Scaled axis's range is determined by the observed correlation, which is better for observing the change the correlation.
                               \n
                               Unscaled axis's range is fixed to [-1, 1], which is better for observing the magnitude of correlation.  ")
                          ),
                          br(), br(),
                          uiOutput("varnames2"),
                          br(), br(),
                          radioButtons("show_cor", label = "Show correlation table?",
                                       choices = c("No", "Yes"), selected = "No"),
                          width = 3),
                       
                       # main panel
                         mainPanel(
                           fluidRow(
                             column(width = 12, md = 6,
                                  htmlOutput("tab2_plot_title"), 
                                    plotOutput("trend_p", height = "45vh")
                                    ),
                             column(width = 12, md = 6,
                                   h4("Empirical correlation"),
                                   plotOutput("cor_trend_p", height = "45vh")
                                   )
                           ),
                           conditionalPanel(
                             condition = "input.show_cor == 'Yes'",
                             br(),
                             h4("Correlation values by time"),
                             div(style = "overflow-x: auto;",
                                 dataTableOutput("cor_tb")
                             )
                           ),
                           width = 9
                         ))),
              
              ## subtab 2.3: overall
              tabPanel(title = "Multivariate analysis",
                         sidebarLayout(
                           sidebarPanel(
                             selectInput("cor_type2", label="Type of correlation",
                                         choices = list("pearson", "spearman")),
                             tagList(
                               icon("info-circle"),
                               em("Correlation measures may be unrealiable when the proportion of missing is large!")
                             ),
                             br(), br(),
                             uiOutput("time_bar1"),
                             uiOutput("varnames2_3"),
                             width = 3),
                           mainPanel(
                             h3("Correlation heatmap"),
                             # square aspect ratio preserved via CSS max-width, but
                             # scales down on narrow viewports instead of a fixed 600px
                             div(style = "max-width: 700px; width: 100%; margin: 0 auto;",
                                 plotOutput("heatmap", height = "55vh")
                             ),
                             h5(icon("info-circle"),
                                "Correlation measures are very sensitive to sample size. 
                                If the number of complete pairs is very small (i.e < 10), 
                                the calculated measures are less realiable and will affect downstream analysis.
                                User may consider removing these time points from the dataset.
                                Details of these time points can be examined in the previous subtabs. 
                                "),
                             width = 9
                           )
                         )
                       )
    )),
  
  # tab 3: network and subnetwork plots
  tabPanel(title = "Dynamic visualization",
           sidebarLayout(
             # side bar
             sidebarPanel(
               # choose correlation type
               selectInput("cor_type3", label="Type of correlation/association", 
                           choices = list("pearson", "spearman", "euclidean")),
               br(),
               radioButtons("lambda", label = "Stability constraint",
                           choices = c("None", "Mild", "Heavy"),
                           selected = "Mild"),
               br(),
               uiOutput("thres_m"), # show edge or not
               br(),
               uiOutput("time_bar"), # time bar
               uiOutput("varnames3"), # select variables
               actionButton("confirm", "Confirm selection"),
               width = 3
             ),
             
             # main panel
             mainPanel(h3("Dynamic network plot"),
                       plotOutput("netp", height = "65vh"),
                       htmlOutput("vis_info"),
                       width = 9
             )
             )
           ),
  
  # tab 4: Overall structure
  tabPanel(title = "Aggregated visualization",
           sidebarLayout(
             sidebarPanel(
               # choose correlation type
               selectInput("mtype2", label="Type of correlation/association",
                           choices = list("pearson", "spearman", "euclidean")),
               br(),
               sliderInput("thres_m2", label = "Show average similarity above",
                           min =0, max = 1, value = 0.1),
               # weight
               checkboxInput("time_wt", label = "Weigh by time interval?", value = T),
               tagList(
                 icon("info-circle"),
                 em("Correlation at each time point is weight by the interval between itself and the next time point.
                    It assumes the data structure stays unchanged until the next measurement point.")
               ),
               # variable list
               uiOutput("varnames4"),
               actionButton("confirm2", "Confirm selection"),
               width = 3
              ),
             
             mainPanel(
               h3("Integrated network of correlation"),
               plotOutput("int_net", height = "65vh"),
               div(style = "overflow-x: auto;",
                   dataTableOutput("test")
               ),
               width = 9
             )
           ))

)


)


 #### Server ####
# NOTE: server logic is unchanged from the previous version except where
# noted below (df_cor refactor, cor_tb table, tab1_note bullet formatting,
# cor_trend_p color gradient) already covered in earlier turns. No server
# changes are required for the responsive-layout work in this turn — the
# fixes are entirely UI/CSS-side (column breakpoints, vh-based heights,
# overflow-safe wrappers for tables/plots).

server <- function(input, output) {
  options(shiny.maxRequestSize=10*1024^2)
  # tab 1: data upload and preview
  ## data upload
  df <- reactive({
    req(input$df_path)
    read.csv(input$df_path$datapath, check.names = F)
  })
  ## data preview
  output$show_df <- DT::renderDT(
    datatable(df(), options = list(scrollX = T, fixedHeader=T)) %>%
      formatSignif(columns = colnames(df()), digits = 3),
    server = T
  )
  output$size_info <- renderDataTable({
    req(df(), input$time_var, input$id_var)
    Nsize <- length(unique(df()[, input$id_var]))
    Nrange <- range(table(df()[, input$id_var]))
    Nrange <- paste0("(", round(Nrange[1], 2), ", ", round(Nrange[2], 2), ")")
    Nfreq <- round(mean(table(df()[, input$id_var])), 2)
    Trange <- range(df()[, input$time_var])
    Trange <- paste0("(", round(Trange[1], 2), ", ", round(Trange[2], 2), ")")
    medT <- round(median(df()[, input$time_var]), 2)
    sum_tb <- data.frame(c("Number of participants",
                           "Range of number of observations per participant",
                           "Average number of observations per participant",
                           "Range of follow up time",
                           "Median of follow up time"),
                         c(Nsize,Nrange, Nfreq, Trange, medT))
    datatable(sum_tb, rownames = FALSE, caption = " ", colnames = c(" ", " "),
              options = list(dom = "t"))
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
  
  output$tab1_note <- renderUI({
    req(input$sum_type)
    
    if (input$sum_type == "Time") {
      note_html <- "<ul style='font-style: italic; font-weight: normal; margin: 0; padding-left: 20px;'>
      <li>Summary statistics are sensitive to the proportion of missing values and can be unreliable 
      if the proportion of missing values is large</li>
      <li>Time points with no more than 10 observations are highlighted in red</li>
    </ul>"
    } else {
      note_html <- "<ul style='font-style: italic; font-weight: normal; margin: 0; padding-left: 20px;'>
      <li>Summary statistics are sensitive to the proportion of missing values and can be unreliable 
      if the proportion of missing values is large</li>
      <li>Individuals with no more than 3 observations are highlighted in red</li>
      <li>*slope is derived by fitting linear regression models along each individual trajectory. 
      It reflects an average rate of linear change across time</li>
    </ul>"
    }
    
    tagList(
      div(style = "display: flex; align-items: flex-start; gap: 6px;",
          icon("circle-info", style = "margin-top: 3px;"),
          HTML(note_html)
      )
    )
  })

  ## summary of single variables
  df_uni <- reactive({
    req(df(), input$select_var1, input$time_var, input$id_var)
    df()[, c(input$id_var, input$time_var, input$select_var1)] %>%
      rename(time=input$time_var, id=input$id_var, var=input$select_var1)
  })

  plot_sum <- reactive({
    req(df_uni())
    if(input$sum_type=="Time"){
      df_sum <- df_uni() %>%
        group_by(time) %>%
        mutate(med=median(var, na.rm = T)) %>%
        mutate(Nmiss = sum(is.na(var)),
               Pctmiss = sum(is.na(var))/length(var))
      t_uniq <- sort(unique(df_sum$time))
      plot_sum <- df_sum %>%
        ggplot()+
        geom_boxplot(aes(x=time, y=var, group=time), outlier.size = 0.8, fill = "grey")+
        geom_line(data = df_sum %>% filter(!is.na(med)), aes(x=time, y=med))+
        scale_x_continuous(breaks = t_uniq, name = input$time_var)+
        labs(x=input$time_var, y=input$select_var1)
    }
    else{
      plot_sum <- df_uni() %>%
        filter(complete.cases(.)) %>%
        ggplot()+
        geom_line(aes(x=time, y=var, group=id), alpha = 0.5, linewidth = 0.5)+
        geom_smooth(aes(x=time, y=var), na.rm = T)+
        labs(x=input$time_var, y = input$select_var1)
    }
    plot_sum
    })
  output$sum_plt <- renderPlot({
    plot_sum()
  })
  ### summary table
  output$sum_tb <- renderDataTable({
    req(df_uni())
    N <- length(unique(df_uni()$id))
    if(input$sum_type=="Time"){
      sum_tb <- df_uni() %>%
        group_by(time) %>%
        summarise(
          Mean = mean(var, na.rm = T),
          SD = sd(var, na.rm = T),
          Min = min(var, na.rm = T),
          Max = max(var, na.rm = T),
          Nmiss = N-sum(!is.na(var))
        )
      datatable(sum_tb, rownames = FALSE, options = list(dom="tp"),
                colnames = c(input$time_var, "Mean", "SD", "Min", "Max", "# of missing values")
                ) %>%
        formatRound(columns = c("Mean", "SD", "Min", "Max"), digits = 2) %>%
        formatStyle("Nmiss", target = "row", backgroundColor = styleInterval((N-10), c(NA,"#ffe6e6")))
    } else {
      sum_tb <- df_uni() %>% group_by(id) %>%
        summarise(
          nobs = sum(!is.na(var)),
          Mean = mean(var, na.rm = T),
          SD = sd(var, na.rm = T),
          Min = min(var, na.rm = T),
          Max = max(var, na.rm = T),
          vel = tryCatch({coef(lm(var~time))[2]},
                         error = function(e){NA})
        )
      datatable(sum_tb, rownames = FALSE, options = list(dom="ftp"),
                colnames = c("ID", "# of observations", "Mean", "SD", "Min", "Max", "Slope*")
                ) %>%
        formatStyle("nobs", target = "row", backgroundColor = styleInterval(3, c("#ffe6e6", NA))) %>%
        formatRound(columns = c("Mean", "SD", "Min", "Max", "vel"), digits = 2)
    }
  })

  ### missing plot: time summarization only
  output$miss_plot <- renderPlot({
    req(df(), input$select_var1, input$time_var, input$id_var)
   df_miss <- df()[, c(input$id_var, input$time_var, input$select_var1)] %>%
      rename(time=input$time_var, id=input$id_var, var=input$select_var1) %>%
      mutate(id = as.factor(id)) %>%
      arrange(time)
   wide_miss <- df_miss %>%
     pivot_wider(id_cols = "id", names_from = "time", values_from = "var") %>%
     select(-id)
   n_total   <- length(as.matrix(wide_miss))
   n_missing <- sum(is.na(wide_miss))
   pct_missing <- round(100 * n_missing / n_total, 1)
   pct_present <- round(100 - pct_missing, 1)
   plt_miss <- wide_miss %>%
       visdat::vis_miss(.)+
       labs(x=paste0(input$time_var, " (% missing)"), y = "ID",
            caption = paste0("Total observed: ", pct_present, "%   |   missing: ", pct_missing, "%"))+
       theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0),
             axis.title.x = element_text(margin = margin(t = 10)),
             plot.caption = element_text(hjust = 0.5, size = 10))
   plt_miss <- plt_miss +
     scale_fill_manual(
       values = c(`FALSE` = "grey80", `TRUE` = "grey20"),
       labels = c(`FALSE` = "Observed", `TRUE` = "Missing"),
       name = ""
     )
   plt_miss
  })

  ## subtab 2.2: pairwise correlaion
  output$varnames2 <- renderUI({
    req(df())
    checkboxGroupInput("select_var2", label = "Variables", choices = colnames(df()),
                       selected = colnames(df())[3:4])
  })

  ### data
  df_pair <- reactive({
    req(df(), input$select_var2, input$time_var, input$id_var)
    validate(need(length(input$select_var2)==2, "Please select a pair of variables."))
    df_pair <- df()[, c(input$time_var, input$id_var, input$select_var2)] %>%
      rename(time=input$time_var, id = input$id_var)
    df_pair
  })
  #### comparision of distribution and trend
  output$tab2_plot_title <- renderPrint({
    if(input$sum_type2 == "Time"){
      h4("Distribution comparision")
    } else {
      h4("Trajectory comparison")
    }
    
  })
  
  output$trend_p <- renderPlot({
    t_uniq <- unique(df_pair()$time)
    if(input$sum_type2 == "Time"){
      p1 <- df_pair() %>%
        pivot_longer(input$select_var2) %>%
        ggplot(aes(x = time, y=value, fill=name, colour = name, group=interaction(time, name)))+
        geom_boxplot(position = "dodge2", alpha = 0.7)+
        scale_fill_brewer(palette = "Set2")+
        scale_color_brewer(palette = "Set2")+
        scale_x_continuous(breaks = t_uniq)+
        labs(title = " ", x = input$time_var, y = " ", col = " ", fill = " ")+
        theme(legend.position = "bottom")
    }
    else {
      p1 <- df_pair() %>%
        pivot_longer(input$select_var2) %>%
        filter(complete.cases(.)) %>%
        ggplot(aes(x = time, y=value, colour = name, group=interaction(id, name)))+
        geom_line(alpha = 0.7, linewidth=0.5)+
        scale_fill_brewer(palette = "Set2")+
        scale_color_brewer(palette = "Set2")+
        labs(title = " ", x = input$time_var, y = " ", col = " ")+
        theme(legend.position = "bottom")
    }
    p1
  })

  ### empirical correlations
  # shared reactive: correlation + pair counts by time (used by both plot and table)
  df_cor <- reactive({
    req(df_pair(), input$select_var2, input$cor_type)
    validate(need(length(input$select_var2) == 2, "Please select a pair of variables."))
    id_uniq <- unique(df_pair()$id)
    N <- length(id_uniq)
    df_pair() %>%
      group_by(time) %>%
      summarize(cor = cor(.data[[input$select_var2[1]]], .data[[input$select_var2[2]]],
                          use = "pairwise.complete.obs", method = input$cor_type),
                Npair = sum(complete.cases(.data[[input$select_var2[1]]], .data[[input$select_var2[2]]]))) %>%
      mutate(Npct = Npair/N)
  })

  output$cor_trend_p <- renderPlot({
    req(df_cor())
    t_uniq <- unique(df_pair()$time)
    p2 <- df_cor() %>%
      filter(complete.cases(.)) %>%
      filter(Npair >= 10) %>%
      ggplot()+
      geom_point(aes(x=time, y=cor, color = Npct), size = 3)+
      geom_line(aes(x=time, y=cor), color = "black")+
      labs(title = "", x = input$time_var, y = " ",
           color = "Proportion of complete pairs")+
      scale_color_gradient(low = "white", high = "black",
                            limits = c(0, 1), guide = "colorbar")+
      theme(legend.position = "bottom")+
      guides(color = guide_colorbar(order = 1, barwidth = 10, barheight = 0.6))+
      scale_x_continuous(breaks = t_uniq)
      if(!input$scaleY){p2 <- p2 + ylim(-1, 1)}
      p2
  })

  # correlation values table (subtab 2.2)
  output$cor_tb <- renderDataTable({
    req(df_cor())
    tb <- df_cor() %>%
      rename(Correlation = cor, `# complete pairs` = Npair) %>%
      select(time, Correlation, `# complete pairs`)
    datatable(tb, rownames = FALSE, options = list(dom = "tp"),
              colnames = c(input$time_var, "Correlation", "# of complete pairs")) %>%
      formatRound(columns = "Correlation", digits = 3) %>%
      formatStyle("# of complete pairs", target = "row",
                  backgroundColor = styleInterval(9, c("#ffe6e6", NA)))
  })

  # subtab 2.3: overall correlation heatmap
  output$time_bar1 <- renderUI({
    req(df(), input$time_var, input$id_var)
    tvec <- sort(unique((df()[ ,input$time_var])))
      sliderTextInput("time_bar1", label = input$time_var, choices = tvec,
                      selected = tvec[1],
                      grid = TRUE)
  })
  output$varnames2_3 <- renderUI({
    req(df())
    var_choice <-  colnames(df() %>% dplyr::select(!c(input$time_var, input$id_var)))
    checkboxGroupInput("select_var2_3", label = "Variables",
                        choices = var_choice, var_choice)
  })
  ### data
  df_multi <- reactive({
    req(df(), input$time_var, input$id_var, input$cor_type2, input$select_var2_3)
    df_multi <- df()[, c(input$time_var, input$id_var, input$select_var2_3)] %>%
      rename(time=input$time_var, id = input$id_var)

  })
  ## heatmap
  output$heatmap <- renderPlot({
    cor_mat <- cor(df_multi() %>% filter(time==input$time_bar1) %>%
                    dplyr::select(-id, -time),
                  method = input$cor_type2, use = "pairwise.complete.obs")
    cor_long <- reshape2::melt(cor_mat)
    names(cor_long) <- c("Var1", "Var2", "correlation")

    ggplot(cor_long, aes(x = Var1, y = Var2, fill = correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(correlation, 2)), size = 3) +
      scale_fill_gradient2(
        low     = "#E41A1C",
        mid     = "white",
        high    = "#377EB8",
        midpoint = 0,
        limits  = c(-1, 1),
        name    = "Correlation"
      ) +
      theme_minimal() +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1),
        axis.title   = element_blank()
      ) +
      coord_fixed()


  })


  # tab 3: temporal network and group labels
  ## variable list
  output$varnames3 <- renderUI({
    req(df(), input$time_var, input$id_var)
    all_vars <- colnames(df() %>% dplyr::select(!c(input$time_var, input$id_var)))
    checkboxGroupInput("select_var3", label = "Variables", choices = all_vars, select = all_vars)
  })
  ## threshold
  output$thres_m <- renderUI({
    req(input$cor_type3)
      sliderInput("thres_m", label = "Show similarity above",
                  min = 0, max = 1, value = 0.5)
  })

  ## data set to analyze
  confirmed <- reactiveVal(NULL)
  observeEvent(input$confirm, {confirmed(input$select_var3)})
  df_net <- reactive({
    req(df(), input$id_var, input$time_var, confirmed())
    df_net <- df()[, c(input$time_var, confirmed())] %>%
      rename(time = input$time_var) 
    df_net <- df_net %>% group_by(time) %>% group_modify(~clean_sparse_columns(.x, min_obs = 10))
    df_net
  })
  ## time axis
  output$time_bar <- renderUI({
    req(df_net(), input$time_var)
    tvec <- sort(unique(df_net()[["time"]]))
      time_bar <- sliderTextInput("time_bar", label = input$time_var, choices = tvec,
                                  selected = tvec[1], grid = TRUE)
    time_bar
  })

  # main panel outputs
  obs_cors <- reactive({
    req(df_net(), input$cor_type3)
    obs_cors <- df_net() %>%
      group_by(time) %>%
      group_map(~{get_similarity(.x, use = "pairwise.complete.obs", method = input$cor_type3)})
    obs_cors
  })
  # LOCF
  filled_obs_cor <- reactive({
    req(obs_cors())
    filled_obs_cor <- obs_cors()
    filled_obs_cor[[1]][is.na(filled_obs_cor[[1]])] <- 1e-5
    for (i in 2:length(filled_obs_cor)) {
      na_mask <- is.na(filled_obs_cor[[i]])
      filled_obs_cor[[i]][na_mask] <- filled_obs_cor[[i - 1]][na_mask]
    }
    filled_obs_cor
  })
  
  t_uniq <- reactive({
    req(df_net())
    sort(unique(df_net()$time))
  })

  grid_search <- reactive({
    req(filled_obs_cor(), t_uniq())
    lambdas <- seq(0, 10, length.out = 100)
    sweep_smooth <- lambda_sweep(filled_obs_cor(), lambdas)
    sweep_smooth
  })

  layout <- reactive({
    req(grid_search(), input$lambda, filled_obs_cor())
    menger_lam <- lcurve_corner_menger(grid_search(), plot = FALSE)
    maxdist_lam <- lcurve_corner_dist(grid_search())
    lam_options <- sort(c(0, maxdist_lam$lambda_star, menger_lam$lambda_star))
    ui_options <- c("None", "Mild", "Heavy")
    lam_star <- lam_options[which(ui_options == input$lambda)]
    dmds_fit <- dyn_mds(obs_sim = filled_obs_cor(), lambda = lam_star, d = 2)
    dmds_fit$embeddings
  })

 # plot
  output$netp <- renderPlot({
    withProgress(
      value = 0, message = "Processing", detail = "This may take a while...",
      {
        req(layout(), input$time_bar, obs_cors(), t_uniq(), df_net())

        tid <- which(input$time_bar == t_uniq())
        layout_i <- layout()[[tid]]
        cor_i <- obs_cors()[[tid]]
        df_i <- df_net()[df_net()$time == input$time_bar, ]

        edges <- which(abs(cor_i) > input$thres_m & upper.tri(cor_i) & !is.na(cor_i), arr.ind = T)
        from <- rownames(cor_i)[edges[, 1]]
        to   <- colnames(cor_i)[edges[, 2]]
        wt   <- cor_i[edges]

        edge_df <- if (length(wt) > 0) {
          data.frame(from = from, to = to, weight = wt)
        } else {
          data.frame(from = character(0), to = character(0), weight = numeric(0))
        }

        net_i <- igraph::graph_from_data_frame(
          d        = edge_df,
          directed = F,
          vertices = data.frame(name = rownames(cor_i))
        )

        E(net_i)$width <- abs(E(net_i)$weight) * 8
        E(net_i)$color <- ifelse(E(net_i)$weight > 0, "steelblue", "tomato")

        miss_var_id <- sapply(df_i[, rownames(cor_i)], function(x) all(is.na(x)))
        V(net_i)$color       <- ifelse(miss_var_id, NA, "lightgrey")
        V(net_i)$frame.color <- "lightgrey"

        incProgress(1)

        net_i <<- net_i
        layout_i <<- layout_i
        tid <<- tid
      })

        plot(net_i, layout = layout_i,
             vertex.size        = 20,
             vertex.label.cex   = 1,
             vertex.color       = V(net_i)$color,
             vertex.frame.color = V(net_i)$frame.color,
             edge.curved        = 0.2,
             margin             = c(0, 0, 0.2, 0))
        if(input$cor_type3 %in% c("pearson", "spearman")){
          legend(
            x = "bottom",
            legend = c("Positive", "Negative"),
            col    = c("steelblue", "tomato"),
            lwd    = 3,
            lty    = 1,
            bty    = "n",
            horiz = TRUE, xpd = TRUE, inset = c(0, -0.15)
          )
        }
  }, res = 96)


  # tab 4: integrated correlation and grouping results
  
  ## sidebar variable list
  output$varnames4 <- renderUI({
    req(df(), input$time_var, input$id_var)
    all_vars <- colnames(df() %>% dplyr::select(!c(input$time_var, input$id_var)))
    checkboxGroupInput("select_var4", label = "Variables",
                       choices = all_vars, selected = all_vars)
  })
  confirmed2 <- reactiveVal(NULL)
  observeEvent(input$confirm2, {confirmed2(input$select_var4)})
  df_net2 <- reactive({
    req(df(), input$id_var, input$time_var, confirmed2())
    df()[, c(input$time_var, confirmed2())] %>%
      rename(time = input$time_var)
  })
  
  obs_cors2 <- reactive({
    req(df_net2(), input$mtype2)
    obs_cors2 <- df_net2() %>%
      group_by(time) %>%
      group_map(~{get_similarity(.x, use = "pairwise.complete.obs", method = input$cor_type3)})
   obs_cors2
  })
  
  t_uniq2 <- reactive({
    req(df_net2())
    sort(unique(df_net2()$time))
  })
  
  wt_vec <- reactive({
    req(t_uniq2())
    nt <- length(t_uniq2())
    if(input$time_wt){
      wt_vec <- c(diff(t_uniq2()), 0)
    } else {
      wt_vec <- rep(1, nt)
    }
    wt_vec
  })
  
   int_cor <- reactive({
     req(obs_cors2(), wt_vec())
     wts <- wt_vec()
     cor_arr <- array(unlist(obs_cors2()),
                       dim = c(nrow(obs_cors2()[[1]]),
                               ncol(obs_cors2()[[1]]),
                               length(obs_cors2())))
     int_cor <- apply(cor_arr, c(1, 2),
                       function(x){sum(x*wts, na.rm = T)/sum(wts)})
     colnames(int_cor) <- colnames(obs_cors2()[[1]])
     rownames(int_cor) <- rownames(obs_cors2()[[1]])
     int_cor
   })
   
  int_diss <- reactive({
    req(int_cor())
    1-int_cor()
  })
  
 int_coords <- reactive({
   req(int_diss())
   smacofSym(int_diss(), type = "ordinal")$conf
 })
 
 output$int_net <- renderPlot({
   withProgress(
     value = 0, message = "Processing", detail = "This may take a while...",
     {
       req(int_coords(), int_cor())

       int_layout <- int_coords()
       int_cor <- int_cor()

       edges <- which(abs(int_cor) > input$thres_m2 & upper.tri(int_cor) & !is.na(int_cor), arr.ind = T)
       from <- rownames(int_cor)[edges[, 1]]
       to   <- colnames(int_cor)[edges[, 2]]
       wt   <- int_cor[edges]

       edge_df <- if (length(wt) > 0) {
         data.frame(from = from, to = to, weight = wt)
       } else {
         data.frame(from = character(0), to = character(0), weight = numeric(0))
       }

       int_net <- igraph::graph_from_data_frame(
         d        = edge_df,
         directed = F,
         vertices = data.frame(name = rownames(int_cor))
       )

       E(int_net)$width <- abs(E(int_net)$weight) * 8
       E(int_net)$color <- ifelse(E(int_net)$weight > 0, "steelblue", "tomato")

       incProgress(1)

       int_net <<- int_net
       int_layout <<- int_layout
     })

   plot(int_net, layout = int_layout,
        vertex.size        = 20,
        vertex.label.cex   = 1,
        vertex.color       = "lightgrey",
        vertex.frame.color = "lightgrey",
        edge.curved        = 0.2,
        margin             = c(0, 0, 0.2, 0))
   if(input$cor_type3 %in% c("pearson", "spearman")){
     legend(
       x = "bottom",
       legend = c("Positive", "Negative"),
       col    = c("steelblue", "tomato"),
       lwd    = 3,
       lty    = 1,
       bty    = "n",
       horiz = TRUE, xpd = TRUE, inset = c(0, -0.15)
     )
   }
 }, res = 96)
 
}



# Run the application 
shinyApp(ui = ui, server = server)
