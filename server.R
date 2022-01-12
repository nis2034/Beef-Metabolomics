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
  
  
  
 ############## Normalization ########################
  
  norm_logic_1 <- reactive({
    if(input$norm_1 == "quot"){
      ref_sam_1 = dQuote(input$ref_samples_1)
      feat_max = input$feat_max_1
      
      D %<>% mt_remove_buffer() %>% mt_pre_norm_quot(feat_max = feat_max, ref_samples = Buffer==ref_sam_1) %>%
      mt_plots_dilution_factor(in_col="Buffer") %>% mt_plots_sample_boxplot(color=Buffer, title = "After probabilistic quotient normalization", plot_logged = T)
      r <- metadata(D)$results
        
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
      r <- norm_logic_1()
      box_p_norm1 <- r[[1]]$output[1]
      tagList(
        selectInput("ref_samples_1", "Reference Sample:",
                    c("80% Meth + 20% H2O", "80% Meth + 20% PBS" ,"50% Meth + 50% H2O","50% Meth + 50% PBS"),selected = NULL,multiple = FALSE,
                    selectize = TRUE,width = '200px',size = NULL),
        textInput("feat_max_1", "Feature Maximum ", value = 0.2, width = '100px', placeholder = NULL),
        renderPlot({
          box_p_norm1
        })  
        
       
        
        
      )
    } else if (input$norm_1 == "external") 
      
    { selectInput("col_name_1", "Numeric-value column:",
                  c("Buffer"="col1", "tissue" = "col2"),selected = NULL,multiple = FALSE,
                  selectize = TRUE,width = '200px',size = NULL)
    }
  })
  
  
  
  output$box_plot_norm <- renderUI({
    
    r <- box()
    box_p <- r[[2]]$output[1]
    renderPlot({
      box_p
    })  
    
  })
  ######## Data Analysis Pipeline ####################
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
  
  
 
  observeEvent(input$sidebarItemExpanded, {
    if (input$sidebarItemExpanded == "Global Statistics") {
      print("updating tab items")
      updateTabItems(session, "main_menu_items", "missingness")
    }
  })
}
