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
source("mt_read_biocrates.R")
source("merge_SE.R")
source("mt_remove_buffer.R")
source("mt_stats_univ_lm_temp.R")
source("mti_format_se_samplewise_temp.R")
source("mt_biocrates_anno_pathways_xls")

source("beef_metabolomics_maplet_fun.R")
getwd()
#setwd("/home/rstudio/home/autonomics/nis2034/beef\ metabolomics")
setwd("/home/rstudio/host/beef_data")
rm(list=ls())

ui <- fluidPage(  titlePanel("Beef Metabolomics Analysis"),
                  tabsetPanel(
                             tabPanel(h4("Data Cleanup"),column(3, navlistPanel(
                               widths = c(12, 12), "",
                               tabPanel("Missingness"),
                               tabPanel("Box Plot"),
                               tabPanel("Normalization")
                             ))),
                             tabPanel(h4("Global Statistics")),
                             tabPanel(h4("Statistical Results Presentation"))
                  )
    
  )



# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
