
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinydashboard)
library(plotly)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Upload Data Set", tabName = "dataSet", 
             icon = icon("cloud-upload")),
    menuItem("Select Gene Identifiers", icon = icon("th"), 
             tabName = "selectIDs"),
    menuItem("Select Samples", tabName = "samples",
             icon = icon("tasks")),
    menuItem("View Heatmaps", tabName = "heatmaps",
             icon = icon("table"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "dataSet",
      fluidRow(
        tabBox(title = "Import Data", width = 12,
          tabPanel("Gene Expression Data", status = "primary",
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
                          ),
                checkboxGroupInput("exprsDcheckGroup", label = NULL, 
                    choices = list("1st row as header" = 1,
                                   "1st column as rownames" = 2,
                                   "Transpose" = 3)
                ),
                br(),
                actionButton("confirmExprsD", "Confirm")
              )
            ),
            fluidRow(
              box(width = 12, title = "Preview",
                tableOutput('previewExprsD')
              )
            )
          ),
          tabPanel("Gene Identification Data", status = "primary",
            fluidRow(
              box(width = 12,
                "Please select the CSV file containing the gene 
                identification information for your data set.", br(),
                "Rows should correspond to the rows in the gene expression
                matrix, and column(s) to identifier type(s).  Use the check
                boxes below to adjust the data.", br(), 
                "Note: Files over 5MB may be very slow to load.", 
                br(), br(),
                fileInput('geneID', label = 'Choose CSV File'),
                checkboxGroupInput("geneIDcheckGroup", label = NULL, 
                  choices = list("1st row as header" = 1,
                                 "1st column as rownames" = 2,
                                 "Transpose" = 3)
                ),
                br(),
                actionButton("confirmGeneID", "Confirm")
              )
            ),
            fluidRow(
              box(width = 12, title = "Preview",
                tableOutput('previewGeneID')
              )
            )
          ),
          tabPanel("Sample Data", status = "primary",
            fluidRow(
              box(width = 12,
                "Please select the CSV file containing the sample 
                information for your data set.", br(),
                "Rows should correspond to columns of the gene expression 
                matrix, and column(s) to catagories of sample information.
                Use the check boxes below to adjust the data if needed, then 
                click the 'Confirm' button to finalize your selection.", br(), 
                "Note: Files over 5MB may be very slow to load.", 
                br(), br(),
                fileInput('sampleD', label = 'Choose CSV File'),
                checkboxGroupInput("sampleDcheckGroup", label = NULL, 
                  choices = list("1st row as header" = 1,
                                 "1st column as rownames" = 2,
                                 "Transpose" = 3)
                ),
                br(),
                actionButton("confirmSampleD", "Confirm")
              )
            ),
            fluidRow(
              box(width = 12, title = "Preview",
                tableOutput('previewSampleD')
              )
            )
          )
        )
      )
    ),
    
    tabItem(tabName = "selectIDs",
      fluidRow(
        box(width = 6,
          uiOutput("dataSetIDselector")
        ),
        box(width = 6,
            uiOutput("AIBSARNAidSelector")
        )
      ),
      fluidRow(
        infoBoxOutput("numDataSetGenes", width = 6),
        infoBoxOutput("numAIBSARNAgenes", width = 6)
      ),
      fluidRow(
        infoBoxOutput("numMatchedGenes", width = 12)
      ),
      fluidRow(
        box(width = 12,
          actionButton("buildScabbard", "Confirm ID selection"))
      )
    ),
    tabItem(tabName = 'samples',
      fluidRow(
        uiOutput("sampleSelector")
      ),
      fluidRow(
        box(width = 12, title = "Preview Sample A",
          uiOutput('previewSampleA')
        )
      ),
      fluidRow(
        box(width = 12, title = "Preview Sample B",
          uiOutput('previewSampleB')
        )
      ),
      fluidRow(
        box(width = 12, title = "Preview Sample C",
          uiOutput('previewSampleC')
        )
      )
    ),
    tabItem(tabName = "heatmaps",
      fluidRow(
        uiOutput('heatmapDisplay')
      )
    )
  )
)


# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "ShinyBrainSABER"),
  sidebar,
  body
)
