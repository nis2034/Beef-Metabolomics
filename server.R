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
 conditional <- reactive({
    if(input$fil_buff_samples == TRUE){
     
      D %<>% mt_remove_buffer() %>% mt_plots_missingness(feat_max=0.5)
      r <- metadata(D)$results
     
    } 
  })
  
  output$conditionalBuffFilt <- renderUI({
    if(input$fil_buff_samples == TRUE){
      renderText( "Missigness after removing buffer")
     
    }
  })
  
  output$conditionalBuffFilt1 <- renderUI({
    if(input$fil_buff_samples == TRUE){
     r <- conditional()
     missing2 <- r[[2]]$output[1]
     renderPlot({
       missing2
     })  
    }
  })
  
  
  output$conditionalBuffFilt2 <- renderUI({
    if(input$fil_buff_samples == TRUE){
      r <- conditional()
      missing_heatmap2 <- r[[2]]$output[2]
      renderPlot({
        missing_heatmap2
      })  
    }
  })
####### Box Plot ################
  box <- reactive({
    if(input$color == "buff"){
      
      D %<>% mt_plots_sample_boxplot(color=Buffer, title = "Buffer Types", plot_logged = T)
      r <- metadata(D)$results
      
    } 
    else
    {
      D %<>% mt_plots_sample_boxplot(color=tissue, title = "Tissue Types", plot_logged = T)
      r <- metadata(D)$results
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
