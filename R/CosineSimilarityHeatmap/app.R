library(shiny)
library(d3heatmap)

# Define UI for application
ui <-  navbarPage(title = "BrainSABER Cosine Similarity of 8pcw sample Heatmap",
                  tabPanel("MYC Group", d3heatmapOutput("MYC")),
                  tabPanel("SHH Group", d3heatmapOutput("SHH")),
                  tabPanel("TYR Group", d3heatmapOutput("TYR")),
                  tabPanel("Combined", d3heatmapOutput("ATRT"))
)

server <- function(input, output) {
  output$MYC <- renderD3heatmap({
    d3heatmap(MYCexprsAgeLoc, scale = "column",Colv = FALSE,xaxis_font_size = "2pt")
  })
  output$SHH <- renderD3heatmap({
    d3heatmap(SHHexprsAgeLoc, scale = "column", Colv = FALSE,xaxis_font_size = "2pt")
  })
  output$TYR <- renderD3heatmap({
    d3heatmap(TYRexprsAgeLoc, scale = "column", Colv = FALSE,xaxis_font_size = "2pt")
  })
  output$MYCc <- renderD3heatmap({
    d3heatmap(MYCexprsAgeLoc, scale = "column",Colv = FALSE,xaxis_font_size = "0px")
  })
  output$SHHc <- renderD3heatmap({
    d3heatmap(SHHexprsAgeLoc, scale = "column", Colv = FALSE,xaxis_font_size = "0px")
  })
  output$TYRc <- renderD3heatmap({
    d3heatmap(TYRexprsAgeLoc, scale = "column", Colv = FALSE,xaxis_font_size = "0px")
  })
  output$allgenes <- renderD3heatmap({
    d3heatmap(texprs6, scale = "column", Colv = FALSE,xaxis_font_size = "2pt",
              yaxis_font_size = "4pt")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
