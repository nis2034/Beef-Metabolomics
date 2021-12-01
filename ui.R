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

ui <- fluidPage(  titlePanel("Beef Metabolomics Analysis"),
                  tabsetPanel(
                    tabPanel(h4("Data Cleanup"),column(3, navlistPanel(
                      widths = c(12, 12), "",
                      tabPanel("Missingness",  mainPanel(
                        plotOutput("missingPlot")
                      )),
                      tabPanel("Box Plot"),
                      tabPanel("Normalization")
                    ))),
                    tabPanel(h4("Global Statistics")),
                    tabPanel(h4("Statistical Results Presentation"))
                  )
                  
)



#r <- metadata(D)$results