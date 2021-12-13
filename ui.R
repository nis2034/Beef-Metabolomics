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
               tabName = "analysis_pipeline" ),
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
          
         
          plotOutput("missingPlot1"),
          plotOutput("missing_heatmapPlot1"),
          #verbatimTextOutput("filter_buff"),
          h3(uiOutput("conditionalBuffFilt")),
          uiOutput("conditionalBuffFilt1"),
          uiOutput("conditionalBuffFilt2")
          
        )
      ),
      tabItem(
        tabName = "box_plot",
          fluidRow(
            h1("Sample Box Plot "),
            selectInput("color", "Select the color variable: ",c("Buffer"="buff", "tissue" = "tiss"),selected = NULL,multiple = FALSE,
                        selectize = TRUE,width = '200px',size = NULL),
            
            uiOutput("box_plot")
          )
      ),
      tabItem(
        tabName = "normalization",
        fluidRow( 
          
          h1(" Data Normalization "),
          selectInput("norm_tab", "Normalization: ",c("NA" = "na", "Quotient Normalization"="quot", "External Sample Annotation" = "external"),selected = NULL,multiple = FALSE,
                      selectize = TRUE,width = '200px',size = NULL)
         
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