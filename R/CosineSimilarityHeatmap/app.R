library(shiny)
library(d3heatmap)

# Define UI for application
ui <-  navbarPage(title = "BrainSABER Cosine Similarity of 8pcw sample Heatmap",
                  tabPanel("MYC Group", d3heatmapOutput("MYC")),
                  tabPanel("SHH Group", d3heatmapOutput("SHH")),
                  tabPanel("TYR Group", d3heatmapOutput("TYR")),
                  tabPanel("At a Glance", d3heatmapOutput("MYCc"),
                           d3heatmapOutput("SHHc"),
                           d3heatmapOutput("TYRc")),
                  tabPanel("Combined", d3heatmapOutput("ATRT"))
)

server <- function(input, output) {
  output$MYC <- renderD3heatmap({
    d3heatmap(t(MYCcosMat), scale = "column",Colv = FALSE,xaxis_font_size = "2pt")
  })
  output$SHH <- renderD3heatmap({
    d3heatmap(t(SHHcosMat), scale = "column", Colv = FALSE,xaxis_font_size = "2pt")
  })
  output$TYR <- renderD3heatmap({
    d3heatmap(t(TYRcosMat), scale = "column", Colv = FALSE,xaxis_font_size = "2pt")
  })
  output$MYCc <- renderD3heatmap({
    d3heatmap(t(MYCcosMat), scale = "column",Colv = FALSE,xaxis_font_size = "0px")
  })
  output$SHHc <- renderD3heatmap({
    d3heatmap(t(SHHcosMat), scale = "column", Colv = FALSE,xaxis_font_size = "0px")
  })
  output$TYRc <- renderD3heatmap({
    d3heatmap(t(TYRcosMat), scale = "column", Colv = FALSE,xaxis_font_size = "0px")
  })
  output$ATRT <- renderD3heatmap({
    d3heatmap(t(ATRTcosMat), scale = "column", Colv = FALSE,xaxis_font_size = "2pt",
              yaxis_font_size = "4pt")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
