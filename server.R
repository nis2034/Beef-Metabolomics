library(tidyverse)
library(autonomics)
library(maplet)
library(magrittr)
library(grid)
library(gridExtra)
library(data.table)
library(shiny)
library(shinydashboard)
# source prototype functions (should appear in Maplet one day)
source("mt_read_biocrates.R")
source("merge_SE.R")
source("mt_remove_buffer.R")
source("mt_stats_univ_lm_temp.R")
source("mti_format_se_samplewise_temp.R")

source("read_beef_data.R")



server <- function(input, output, session) {
  
  
  
  # observe({
  #   cat("input$main_menu_items ", input$main_menu_items, "\n")
  #   cat("input$data_cleanup ", input$data_cleanup, "\n")
  #   cat("input$global_statistics ", input$global_statistics, "\n")
  #   cat("input$statistical_results_presentation ", input$statistical_results_presentation, "\n")
  #   cat("input$missingness ", input$missingness, "\n")
  #   cat("input$box_plot ", input$box_plot, "\n")
  #   cat("input$normalization ", input$normalization, "\n")
  
  
  D <- read_beef_data()
  D %<>% mt_plots_missingness(feat_max=0.5)
  
  r <- metadata(D)$results
  
  missing1 <- r[[1]]$output[1]
  missing_heatmap1 <- r[[1]]$output[2]
  output$missingPlot1 <- renderPlot({
    missing1
  })
  output$missing_heatmapPlot1 <- renderPlot({
    missing_heatmap1
  })
  
  reactive({
    if(input$fil_buff_samples == TRUE){
  test = "test proved"
    D %<>% mt_remove_buffer() %>% mt_plots_missingness(feat_max=0.5)
    r <- metadata(D)$results
    missing2 <- r[[2]]$output[1]
    missing_heatmap2 <- r[[2]]$output[2]
 
    
  
  
  } 
  })
  
  
  
  output$conditionalBuffFilt1 <- renderUI({
    if(input$fil_buff_samples == TRUE){
      renderText(test)
     
    }
  })
  
  output$conditionalBuffFilt2 <- renderUI({
    if(input$fil_buff_samples == TRUE){
      
      renderPlot({
        missing2
      })  
    }
  })
  output$conditionalBuffFilt3 <- renderUI({
    if(input$fil_buff_samples == TRUE){
      
      renderPlot({
        missing_heatmap2
      })  
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
