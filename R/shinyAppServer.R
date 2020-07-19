#' server side logic for the shiny app
#' 
#' @param input shiny server input, provided automatically
#' @param output shiny server output, provided automatically
#' 
#' @return a shiny server object
#' 
#' @import shiny
#' @import shinydashboard
#' @import shinycssloaders
#' @import plotly
#' @importFrom grDevices dev.off png terrain.colors
#' @importFrom utils write.csv
#' @importFrom ComplexHeatmap Heatmap draw

# Define server logic
shinyAppServer <- function(input, output) {
  #Load or download AIBSARNA object############################################################
  
  #Define values variable to hold reactive data
  values <- reactiveValues()
  #Set up counter variable for button
  values$n <- 0
  values$n1 <- 0
  
  #If AIBSARNA file is uploaded
  observeEvent(input$AIBSARNAfile, {
    e = new.env()
    AIBSARNA <- load(input$AIBSARNAfile$datapath, envir = e)
    values$AIBSARNA <- e[[AIBSARNA]]
  })
  
  output$spinner<-renderText({
    validate(need(input$AIBSARNAdownload > values$n1, "")) #I'm sending an empty string as message.
    values$n1 <- values$n1 + 1
    values$AIBSARNA <- buildAIBSARNA()
    return("Done")
    Sys.sleep
  })
  
  ##Upload Data Set############################################################################
  
  #When dataset path entered assign to exprs in values
  observeEvent(input$exprsD, {
    values$exprs <- read.csv(input$exprsD$datapath,
                             header = TRUE,
                             stringsAsFactors = FALSE)
  })
  
  #Output the number of matched genes
  output$numMatchedGenes <- renderInfoBox({
    if(!is.null(values$AIBSARNA) && !is.null(values$exprs)){
      values$cs <- mapIDs(dataSet = values$exprs, AIBSARNA = values$AIBSARNA)
      infoBox(
        title = '',
        value = paste0(length(values$cs), "/2000"),
        subtitle = "genes mapped",
        icon = icon("check"), color = "green", width = 6
      )} else {
        infoBox(title = "", subtitle = "Please upload AIBSARNA and dataset")
      }
  })
  
  #Confirm dataset button and begin distance calculation
  output$logo<-renderText({
    validate(need(input$confirmExprsD > values$n, "")) #I'm sending an empty string as message.
    values$n <- values$n + 1
    if(!isTruthy(values$AIBSARNA) || !isTruthy(values$exprs)){
      return("Please upload AIBSARNA and dataset")
    } else {
      values$distance <- calcDistance(dataSet = values$cs, AIBSARNA = values$AIBSARNA)
      return("")
    }
  })
  
  ##View Sample Similarity ############################################################
  
  #Sample selection
  output$dataSetSelector <- renderUI({
    validate(need(values$distance, "Please confirm data set upload"))
    dataIDs <- colnames(values$exprs[-1])
    selectInput("dataSetSelected",
                label = "Choose a sample",
                choices = dataIDs)
  })
  
  #Data table
  output$recom <- DT::renderDataTable({
    if(isTruthy(input$dataSetSelected) && isTruthy(values$distance)){
      values$distanceSet <- values$distance[[input$dataSetSelected]]
      m <- values$distanceSet$meta[order(values$distanceSet$meta$Tau, decreasing =TRUE),]
      m <- cbind(1:nrow(m), m)
      colnames(m)[1] <- "Rank"
      
      dt <- DT::datatable(m, rownames =FALSE, selection = "single")
      dt }
  }, server =FALSE)
  
  #Download for chart button
  output$download <- downloadHandler(filename = 'shinyBrainSABERsimilarity.csv',
                                     content = function(file) {
                                       m <- values$distanceSet$meta[order(values$distanceSet$meta$Tau, decreasing =TRUE),]
                                       m <- cbind(1:nrow(m), m)
                                       colnames(m)[1] <- "Rank"
                                       write.csv(m, file, row.names =FALSE)
                                     })
  #Kendall's tau heatmap
  output$heat <- renderPlotly({
    if(isTruthy(input$dataSetSelected) && isTruthy(values$distance)){
      t <- values$distanceSet$correlation
      
      plot_ly(x = colnames(t), y = rownames(t), z = t, type = "heatmap") %>%
        layout(title = "Kendall's Tau Heatmap",
               xaxis = list(categoryorder = "trace"))
    }
  })
  
  #Spearman's Rho heatmap
  output$heat2 <- renderPlotly({
    if(isTruthy(input$dataSetSelected) && isTruthy(values$distance)){
      t1 <- values$distanceSet$scorrelation
      
      plot_ly(x = colnames(t1), y = rownames(t1), z = t1, type = "heatmap") %>%
        layout(title = "Spearman's Rho Heatmap",
               xaxis = list(categoryorder = "trace"))
    }
  })
  
  ##View Data Set Similarity ############################################################
  
  #Kendall's Tau heatmap build
  kendallHeatplot <- reactive({
    n = length(values$distance)
    k = matrix(unlist(sapply(values$distance, function(k){k$meta$Tau})), ncol = n)
    colnames(k) <- names(values$distance)
    
    Heatmap(k, col = terrain.colors(n = 10),
            name = 'Tau', show_row_names =FALSE,
            column_title = "Clustered Kendall's Tau Heatmap")
  })
  
  #Output kendalls heatmap
  output$heatplot <- renderPlot({
    validate(need(values$distance, "Please confirm data set upload"))
    draw(kendallHeatplot())
  })
  
  #Download handler for heatmap
  output$downloadkh <- downloadHandler(filename = "KendallHeatmap.png",
                                       content = function(file){
                                         png(file, type = 'cairo')
                                         draw(kendallHeatplot())
                                         dev.off()
                                       },
                                       contentType = 'image/png')
  
  #Spearman's heatmap build function
  spearmanHeatplot <- reactive({
    n = length(values$distance)
    s = matrix(unlist(sapply(values$distance, function(s){s$meta$Rho})), ncol = n)
    colnames(s) <- names(values$distance)
    
    Heatmap(s, col = terrain.colors(n = 10),
            name = 'Rho', show_row_names =FALSE,
            column_title = "Clustered Spearman's Rho Heatmap")
  })
  
  
  #Draw Spearman's heatmap
  output$heatplot2 <- renderPlot({
    validate(need(values$distance, "Please confirm data set upload"))
    draw(spearmanHeatplot())
  })
  
  #Download handler for Spearman's heatmap
  output$downloadsh <- downloadHandler(filename = "SpearmanHeatmap.png",
                                       content = function(file){
                                         png(file, type = 'cairo')
                                         draw(spearmanHeatplot())
                                         dev.off()
                                       },
                                       contentType = 'image/png')
  
}

