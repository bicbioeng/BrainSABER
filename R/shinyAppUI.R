#' Html blocks for building shiny app User Interface
#' 
#' @import shiny
#' @import shinydashboard
#' @import shinycssloaders
#' @import plotly
#' @importFrom ComplexHeatmap Heatmap draw

# Define UI for application that draws a histogram
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Upload Data Set", tabName = "dataSet",
             icon = icon("cloud-upload")),
    menuItem("View Sample Similarity", tabName = "similarity",
             icon = icon("table")),
    menuItem("View Data Set Similarity", tabName = "similarity2",
             icon = icon("chart-bar"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "dataSet",
            fluidRow(
              tabBox(title = "Import Data", width = 12,
                     tabPanel("Gene Expression Data", status = "primary",
                              fluidRow(box(width = 12,
                                           fluidRow("Please select the location of the AIBSARNA R object or download",
                                                    fileInput('AIBSARNAfile', label = "Choose R object")),
                                           fluidRow(
                                             h5("-OR-", .noWS = "inside")
                                           ),
                                           fluidRow(
                                             actionButton('AIBSARNAdownload', label = 'Download'),
                                             withSpinner(uiOutput("spinner"),
                                                         type = 5, size = 0.25,
                                                         proxy.height = "30px")))
                              ),
                              
                              
                              fluidRow(
                                box(width = 12,
                                    "Please select the CSV file containing a numeric gene
                                    expression matrix for your data set.", br(),
                                    "Rows should correspond to genes, and columns to samples.
                                    Use the check boxes below to adjust the data.", br(),
                                    "Note: Files over 5MB may be very slow to load.",
                                    br(), br(),
                                    fileInput('exprsD', label = 'Choose CSV File',
                                              accept = c(
                                                "text/csv",
                                                "text/comma-separated-values,text/plain",
                                                ".csv")
                                    )
                                )
                              ),
                              fluidRow(
                                infoBoxOutput("numMatchedGenes", width = 12)),
                              br(),
                              fluidRow(
                                p(class = 'text-center', actionButton('confirmExprsD', 'Confirm'))
                              ),
                              fluidRow(withSpinner(uiOutput("logo"),
                                                   type = 5, size = 0.25,
                                                   proxy.height = "30px"))
                     )))
    ),
    tabItem(tabName = "similarity",
            fluidRow(
              box(width = 12,
                  uiOutput("dataSetSelector")
              )
            ),
            fluidRow(box(width = 12,
                         DT::dataTableOutput("recom")
            )),
            fluidRow(
              p(class = 'text-center', downloadButton('download', 'Download Data'))
            ),
            fluidRow(box(width = 12,
                         plotlyOutput("heat"))
            ),
            fluidRow(box(width = 12,
                         plotlyOutput("heat2")))
            
    ),
    tabItem(tabName = "similarity2",
            fluidRow(box(width = 12,
                         plotOutput("heatplot"))),
            fluidRow(p(class = 'text-center', downloadButton('downloadkh', 'Download Heatmap'))),
            fluidRow(box(width = 12,
                         plotOutput("heatplot2"))),
            fluidRow(p(class = 'text-center', downloadButton('downloadsh', 'Download Heatmap')))
    )
  )
)


# Put them together into a dashboardPage
shinyAppUI <- function(){
    dashboardPage(
    dashboardHeader(title = "shinyBrainSABER"),
    sidebar,
    body
  )
}
