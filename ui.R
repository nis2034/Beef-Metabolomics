<<<<<<< HEAD
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

=======
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
          
          #verbatimTextOutput("filt_buff"),
          plotOutput("missingPlot1"),
          plotOutput("missing_heatmapPlot1"),
          uiOutput("conditionalBuffFilt1"),
          uiOutput("conditionalBuffFilt2"),
          uiOutput("conditionalBuffFilt3")
         
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
>>>>>>> 6bb699b3b8981ab651f82787aecc845d03d2c4d0
