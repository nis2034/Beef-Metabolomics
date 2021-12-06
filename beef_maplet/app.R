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
setwd("/home/rstudio/host/beef_data")
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
          
          
         
          verbatimTextOutput("filter_buff")
         
          
        )
      ),
      tabItem(
        tabName = "box_plot",
        fluidRow(
          h1("box_plot 2")
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
  
  # observe({
  #   cat("input$main_menu_items ", input$main_menu_items, "\n")
  #   cat("input$data_cleanup ", input$data_cleanup, "\n")
  #   cat("input$global_statistics ", input$global_statistics, "\n")
  #   cat("input$statistical_results_presentation ", input$statistical_results_presentation, "\n")
  #   cat("input$missingness ", input$missingness, "\n")
  #   cat("input$box_plot ", input$box_plot, "\n")
  #   cat("input$normalization ", input$normalization, "\n")
  

 
 # output$filter_buff <- reactive(input$fil_buff_samples)
 
  
  
  observeEvent(input$fil_buff_samples, { 
    if(input$fil_buff_samples == TRUE){
   
    output$filter_buff ="Checked"
    } 
  })
  
  
  
  
 
  # })
  observeEvent(input$sidebarItemExpanded, {
    if (input$sidebarItemExpanded == "Global Statistics") {
      print("updating tab items")
      updateTabItems(session, "main_menu_items", "missingness")
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
