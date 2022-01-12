#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


# libraries
library(tidyverse)
library(autonomics)
library(maplet)
library(magrittr)
library(grid)
library(gridExtra)
library(data.table)
library(shiny)
library("shinydashboard")
# source prototype functions (should appear in Maplet one day)

#setwd("/home/rstudio/home/autonomics/nis2034/beef\ metabolomics")
#setwd("/home/rstudio/host/beef_data")
rm(list=ls())

ui <- dashboardPage(
  dashboardHeader(title = "Beef Metabolomics Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      id = "main_menu_items",
      menuItem("Data Cleanup", tabName = "data_cleanup",
               menuSubItem("Missingness", tabName = "missingness"),
               menuSubItem("Box Plot", tabName = "box_plot"),
               menuSubItem("Normalization", tabName = "normalization")),
      menuItem("Global Statistics",
               tabName = "global_statistics",
               menuSubItem("Sub Menu Item 1", tabName = "sub_1"),
               menuSubItem("Sub Menu Item 2", tabName = "sub_2")),
      menuItem("Data Analysis Pipeline",
               tabName = "analysis_pipeline"
      ),
      menuItem("Statistical Results Presentation",
               tabName = "global_statistics",
               menuSubItem("Sub Menu Item 1", tabName = "sub_1"),
               menuSubItem("Sub Menu Item 2", tabName = "sub_2")
      )
    )),
  
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "missingness",
        fluidRow(
          checkboxInput("fil_buff_samples", "Filter Buffer Samples", FALSE),
          h3("Missigness before removing buffer"),
          
          
          
          verbatimTextOutput("filter_buff"),
          
          uiOutput("conditional1"),
          h3(uiOutput("conditionalBuffFilt")),
          uiOutput("conditional2"),
          uiOutput("conditional3")
          
        )
      ),
      tabItem(
        tabName = "box_plot",
        fluidRow(
          h1("Sample Box Plot "),
          selectInput("color", "Select the color variable: ",c("Buffer"="buff", "tissue" = "tiss"),selected = NULL,multiple = FALSE,
                      selectize = TRUE,width = '200px',size = NULL),
          
          verbatimTextOutput("col")
        )
      ),
      tabItem(
        tabName = "normalization",
        fluidRow( 
          #numericInput("feat_max_1", label = "Feature Maximum ", value = 0.2, width = '100px'),
          h1(" Data Normalization "),
          selectInput("norm_1", "Normalization: ",c("NA" = "na", "Quotient Normalization"="quot", "External Sample Annotation" = "external"),selected = NULL,multiple = FALSE,
                      selectize = TRUE,width = '200px',size = NULL),
          uiOutput("norm_sub_menu_1")
          
          
        )
      ),
      tabItem(
        tabName = "analysis_pipeline",
        fluidRow( 
          
          h1("Data Analysis Pipeline "),
          
          radioButtons("fil_buff_samples_pip", "Do you want to filter out buffer samples:",
                       c("Yes" = "yes",
                         "No" = "no")),style = "margin-left: 5px;",
          
          
          textInput("filt_perc", "Enter the % of missingness for pre-filter: ", value = 0.5, width = '250px', placeholder = NULL),
          
         
          selectInput("norm", "Normalization: ",c("NA" = "na", "Quotient Normalization"="quot", "External Sample Annotation" = "external"),selected = NULL,multiple = FALSE,
                      selectize = TRUE,width = '200px',size = NULL),
          uiOutput("norm_sub"),
          
          selectInput("trans", "Data Transformation: ",c("NA" = "na", "log2"="log2", "Exponentiate" = "exp", "Sample Relative" = "rel", "Scale data" = "scale"),selected = NULL,multiple = FALSE,
                      selectize = TRUE,width = '200px',size = NULL),
          
          selectInput("impute", "Impute Data Using: ",c("NA" = "na", "knn"="knn", "Minimum value" = "min"),selected = NULL,multiple = FALSE,
                      selectize = TRUE,width = '200px',size = NULL)
        )
      ),
      tabItem(
        tabName = "sub_1",
        fluidRow(
          h1("Sub Menu Page 1")
        )
      ),
      tabItem(
        tabName = "sub_2",
        fluidRow(
          h1("Sub Menu Page 2")
        )
      )
    )
  )
)
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  check_filtr <- FALSE 

  
  f1 <-  reactive({ input$fil_buff_samples})
  
  output$f2 <- renderText(f1)

  
  test1 =""
  test2 =""
  conditional <- reactive({
    if(input$fil_buff_samples == TRUE){
      test1 <- "test proved"
      
      test2 <- "test proved for second var"
      
      
    }
  })
  
  
  
  output$conditional1 <- renderUI({
    if(input$fil_buff_samples == TRUE){
      renderText(conditional())
      
    }
  })
  
  output$conditionalBuffFilt <- renderUI({
    if(input$fil_buff_samples == TRUE){
      renderText(HTML(paste0("<b>","Missigness after removing buffer", "</b>")))
      
    }
  })
  output$conditional2 <- renderUI({
    if(input$fil_buff_samples == TRUE){
      renderText(test1)
      
    }
  })
  
  output$conditional3 <- renderUI({
    if(input$fil_buff_samples == TRUE){
      renderText(test2)
      
    }
  })
  
  f2 <-  reactive({ input$color})
  output$col <- renderUI({
    renderPrint(f2)})
  
  output$col <- reactive({
    if(input$color == "buff")
      "Buffer"
    else
      "tissue"
  })
  
  # output$quot < -reactive({
  #   input$norm == "quot"
  # })
  # 
  # output$external < -reactive({
  #   input$norm == "external"
  # })
  # 
  # })
  ##### Data analysis pipelline ########
  
  output$norm_sub1 <- renderUI({
    selectInput("ref_samples", "Reference Sample:",
                c("80% Meth + 20% H2O"="comb1", "80% Meth + 20% PBS" = "comb2","50% Meth + 50% H2O" = "comb3","50% Meth + 50% PBS" ="comb4"),selected = NULL,multiple = FALSE,
                selectize = TRUE,width = '200px',size = NULL)
    
  })
  

  
  output$norm_sub <- renderUI({
    if(input$norm == "quot") 
    {
      selectInput("ref_samples", "Reference Sample:",
                  c("80% Meth + 20% H2O"="comb1", "80% Meth + 20% PBS" = "comb2","50% Meth + 50% H2O" = "comb3","50% Meth + 50% PBS" ="comb4"),selected = NULL,multiple = FALSE,
                  selectize = TRUE,width = '200px',size = NULL)
      
    } else if (input$norm == "external") 
      
    { selectInput("col_name", "Numeric-value column:",
                  c("Buffer"="col1", "tissue" = "col2"),selected = NULL,multiple = FALSE,
                  selectize = TRUE,width = '200px',size = NULL)
    }
  })
  
  ############## Normalization ########################
  norm_logic_1 <- reactive({
    if(input$norm_1 == "quot"){
      ref_1 = dQuote(input$ref_samples_1)
      ref_sam_1 = paste0("Buffer ==", ref_1)
      feat_max = input$feat_max_1
      #D %<>% mt_remove_buffer() %>% mt_pre_norm_quot(feat_max = 0.2, ref_samples = Buffer=="80% Meth + 20% H2O") %>%
      #   
      #   mt_plots_dilution_factor(in_col="Buffer") %>% mt_plots_sample_boxplot(color=Buffer, title = "After probabilistic quotient normalization", plot_logged = T)
      # r <- metadata(D)$results
      
    } 
    # else
    # {
    #   D %<>% mt_remove_buffer() %>% mt_pre_norm_external(col_name='Buffer') %>%
    #     
    #     mt_plots_dilution_factor(in_col="Buffer") %>% mt_plots_sample_boxplot(color=Buffer, title = "After probabilistic quotient normalization", plot_logged = T)
    #   r <- metadata(D)$results
    # }
  })
 
  output$norm_sub_menu_1 <- renderUI({
    if(input$norm_1 == "quot") 
    {
      feat_max <- norm_logic_1()
      tagList(
        # current issue is the components below are not changing in the window.
        selectInput("ref_samples_1", "Reference Sample:",
                    c("80% Meth + 20% H2O", "80% Meth + 20% PBS" ,"50% Meth + 50% H2O","50% Meth + 50% PBS"),selected = NULL,multiple = FALSE,
                    selectize = TRUE,width = '200px',size = NULL),
        numericInput("feat_max_1", label = "Feature Maximum ", value = 0.2, width = '100px'),
        hr(),
        textInput("ref_sam_1","Feature Maximum ", value = feat_max)
     
   
      )
    } else if (input$norm_1 == "external") 
      
    { selectInput("col_name_1", "Numeric-value column:",
                  c("Buffer"="col1", "tissue" = "col2"),selected = NULL,multiple = FALSE,
                  selectize = TRUE,width = '200px',size = NULL)
    }
  })
  
  
  
  output$box_plot <- renderUI({
    
    r <- box()
    box_p <- r[[2]]$output[1]
    renderPlot({
      box_p
    })  
    
  })
  observeEvent(input$sidebarItemExpanded, {
    if (input$sidebarItemExpanded == "Global Statistics") {
      print("updating tab items")
      updateTabItems(session, "main_menu_items", "missingness")
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
