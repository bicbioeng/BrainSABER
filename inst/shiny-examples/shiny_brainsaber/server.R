
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(Biobase)
library(BrainSABER)
library(fastcluster)
library(heatmaply)
library(plotly)
#library(rio)
options(shiny.maxRequestSize=200*1024^2)

shinyServer(function(input, output) {
  getFile <- function(inFile, header, rowNames, transpose){
    if (header && rowNames && transpose) {
      data <- read.csv(inFile$datapath, 
                       header = FALSE, 
                       stringsAsFactors = FALSE)
      data <- data.table::transpose(data)
      colnames(data) <- data[1, ]
      row.names(data) <- data[, 1]
      data <- data[-1, -1]
    } else if (header && transpose) {
      data <- read.csv(inFile$datapath, 
                       header = FALSE, 
                       stringsAsFactors = FALSE)
      data <- data.table::transpose(data)
      colnames(data) <- data[1, ]
      data <- data[-1, ]
      row.names(data) <- as.character(1:length(data[, 1]))
    } else if (rowNames && transpose) {
      data <- read.csv(inFile$datapath, 
                       header = FALSE, 
                       stringsAsFactors = FALSE)
      data <- data.table::transpose(data)
      row.names(data) <- data[, 1]
      data <- data[, -1] 
      colnames(data) <- as.character(1:length(data[1, ]))
    } else if (transpose) {
      data <- data.table::transpose(read.csv(inFile$datapath, 
                                             header = FALSE, 
                                             stringsAsFactors = FALSE))
      row.names(data) <- as.character(1:length(data[, 1]))
      colnames(data) <- as.character(1:length(data[1, ]))
    } else if (rowNames && header){
      data <- read.csv(inFile$datapath, 
                       header = header, 
                       stringsAsFactors = FALSE)
      row.names(data) <- data[, 1]
      data <- data[, -1]
    } else if (rowNames){
      data <- read.csv(inFile$datapath, 
                       header = FALSE, 
                       stringsAsFactors = FALSE)
      row.names(data) <- data[, 1]
      data <- data[, -1]
      colnames(data) <- as.character(1:length(data[1, ]))
    } else if (header) {
      data <- read.csv(inFile$datapath, 
                       header = header, 
                       stringsAsFactors = FALSE)
      row.names(data) <- as.character(1:length(data[, 1]))
    } else {
      data <- read.csv(inFile$datapath, 
                       header = FALSE, 
                       stringsAsFactors = FALSE)
      row.names(data) <- as.character(1:length(data[, 1]))
      colnames(data) <- as.character(1:length(data[1, ]))
    }
    data
  }
  
  values <- reactiveValues()
  getExprsD <- reactive({
    inFile <- input$exprsD
    if (is.null(inFile)) return(NULL)
    header <- 1 %in% input$exprsDcheckGroup
    rowNames <- 2 %in% input$exprsDcheckGroup
    transpose <- 3 %in% input$exprsDcheckGroup
    getFile(inFile, header, rowNames, transpose)
  })
  output$previewExprsD <- renderTable(rownames = TRUE,{
    head(getExprsD())
  })
  observeEvent(input$confirmExprsD, {
    values$exprs <- getExprsD()
  })
  
  getGeneID <- reactive({
    inFile <- input$geneID
    if (is.null(inFile)) return(NULL)
    header <- 1 %in% input$geneIDcheckGroup
    rowNames <- 2 %in% input$geneIDcheckGroup
    transpose <- 3 %in% input$geneIDcheckGroup
    getFile(inFile, header, rowNames, transpose)
  })
  output$previewGeneID <- renderTable(rownames = TRUE,{
    head(getGeneID())
  })
  observeEvent(input$confirmGeneID, {
    values$fd <- getGeneID()
  })
  
  getSampleD <- reactive({
    inFile <- input$sampleD
    if (is.null(inFile)) return(NULL)
    header <- 1 %in% input$sampleDcheckGroup
    rowNames <- 2 %in% input$sampleDcheckGroup
    transpose <- 3 %in% input$sampleDcheckGroup
    getFile(inFile, header, rowNames, transpose)
  })
  output$previewSampleD <- renderTable(rownames = TRUE,{
    head(getSampleD())
  })
  observeEvent(input$confirmSampleD, {
    values$pd <- getSampleD()
  })

  getExpressionSet <- reactive({
    if (!isTruthy(values$exprs) || !isTruthy(values$fd) 
         || !isTruthy(values$pd)){
      return(NULL)
    } else {
    exprs <- as.matrix(values$exprs)
    fd <- as(values$fd, "AnnotatedDataFrame")
    pd <- as(values$pd, "AnnotatedDataFrame")
    #force the colnames, since read.csv with header=TRUE doesn't seem to 
    #want to play nice
    colnames(exprs) <- row.names.data.frame(values$pd)
    ExpressionSet(assayData = exprs, phenoData = pd, 
                           featureData = fd)
    }
  })
  output$dataSetIDselector <- renderUI({
    if (!isTruthy(values$exprs) || !isTruthy(values$fd) 
        || !isTruthy(values$pd)){
      h4( "Please confirm all three components of your data set in the 
          'Upload Data Set' tab before proceeding.")
    } else {
      dataIDs <- colnames(values$fd)
      selectInput("dataSetID", 
                  label = "Choose a gene identifier from your data set",
                  choices = dataIDs)
    }
  })
  
  output$AIBSARNAidSelector<- renderUI({
    AIBSARNAids <- colnames(fData(AIBSARNA::AIBSARNA))
    selectInput("AIBSARNAid", 
                label = "Choose a gene identifier from Brain Atlas",
                choices = AIBSARNAids)
  })
  output$numDataSetGenes <- renderInfoBox({
    infoBox(
      "Genes in data set:", length(getGeneID()[, 1]), 
      icon = icon("arrow-up"), color = "light-blue", width = 6
    )
  })
  output$numMatchedGenes <- renderInfoBox({
    if (length(input$dataSetID) < 1
        || length(input$AIBSARNAid) < 1){
      infoBox("Please select gene identifiers.", width = 12,
              icon = icon("arrow-up"), color = "red")
    } else {
    infoBox(
      "Matched genes:", 
      length(
        getExternalVector(dataSet = getExpressionSet(),
                          index = 1,
                          dataSetId = as.character(input$dataSetID),
                          AIBSARNAid = as.character(input$AIBSARNAid))), 
      icon = icon("check"), color = "green", width = 12
    )
    }
  })
  output$numAIBSARNAgenes <- renderInfoBox({
    infoBox(
      "Genes in Brain Atlas:", 
      length(fData(AIBSARNA::AIBSARNA)[, 1]), 
      icon = icon("arrow-up"), color = "yellow", width = 6
    )
  })
  observeEvent(input$buildScabbard, {
    withProgress(message = 'Assembling Data', value = .1, {
      values$eset <- getExpressionSet()
      incProgress(amount = .9)
    })
    withProgress(message = 'Gathering relevant Brain Atlas genes', value = .1, {
      values$dataSetID <- input$dataSetID
      values$AIBSARNAid <- input$AIBSARNAid
      values$relevantGenes <- getRelevantGenes(values$eset, input$dataSetID,
                                               input$AIBSARNAid)
      incProgress(amount = .9)
    })
    withProgress(message = 'Trimming data genes to match Brain Atlas genes', 
      value = .1, {
        values$eset <- getTrimmedExternalSet(values$eset, input$dataSetID,
                                               input$AIBSARNAid)
        incProgress(amount = .9)
    })
  })
  
  output$sampleSelector <- renderUI({
    if (!isTruthy(values$exprs) || !isTruthy(values$fd) 
        || !isTruthy(values$pd)){
      box(width = 12,
        "Please confirm all three components of your data set in the 
        'Upload Data Set' tab before proceeding.")
    } else {
      sampleIDs <- sampleNames(values$eset)
      box(width = 12,
        "Please select the samples from your data from which you would like 
        to view similarity heatmaps, then click 'Confirm Samples' to generate 
        them.",
        br(), br(),
        selectInput("sampleA", 
          label = "Select Sample A", 
          choices = sampleIDs),
        br(),
        selectInput("sampleB", 
          label = "Select Sample B", 
          choices = sampleIDs),
        br(),
        selectInput("sampleC", 
          label = "Select Sample C", 
          choices = sampleIDs),
        br(),
        actionButton('generateHeatmaps', 'Confirm Samples')
      )
    }
  })
  
  output$previewSampleA <- renderUI({
    if(!isTruthy(input$sampleA)){
      box(width = 12,
          "Please select Sample A")
    } else {
      tableOutput("psA")
    }
  })
  output$psA <- renderTable(rownames = TRUE, {
    pData(values$eset)[input$sampleA, ]
  })
  
  output$previewSampleB <- renderUI({
    if(!isTruthy(input$sampleB)){
      box(width = 12,
          "Please select Sample B")
    } else {
      tableOutput("psB")
    }
  })
  output$psB <- renderTable(rownames = TRUE, {
    pData(values$eset)[input$sampleB, ]
  })
  
  output$previewSampleC <- renderUI({
    if(!isTruthy(input$sampleC)){
      box(width = 12,
          "Please select Sample C")
    } else {
      tableOutput("psC")
    }
  })
  output$psC <- renderTable(rownames = TRUE, {
    pData(values$eset)[input$sampleC, ]
  })
  
  observeEvent(input$generateHeatmaps, {
    withProgress(message = 'Processing Sample Selections', value = .1, {
      if(isTruthy(input$sampleA)){
        withProgress(message = 'Getting similarity vectors for Sample A',
          value = .1, {
            values$sampleAidx <- which(
              rownames(pData(values$eset)) %in% input$sampleA, arr.ind = TRUE)
            values$sampleA <- getExternalVector(values$eset, 
              index = values$sampleAidx, dataSetId = values$dataSetID,
              AIBSARNAid = values$AIBSARNAid)
            incProgress(amount = .3)
            values$sampleAeuc <- getSimScores(values$sampleA, 
              relevantGenes = values$relevantGenes, 
              similarity_method = "euclidean")
            incProgress(amount = .3)
            values$sampleAcos <- getSimScores(values$sampleA,
              relevantGenes = values$relevantGenes,
              similarity_method = "cosine")
            incProgress(amount = .3)
          })
      }
      incProgress(amount = .3)
      if(isTruthy(input$sampleB)){
        withProgress(message = 'Getting similarity vectors for Sample B',
          value = .1, {
            values$sampleBidx <- which(
              rownames(pData(values$eset)) %in% input$sampleB, arr.ind = TRUE)
            values$sampleB <- getExternalVector(values$eset, 
              index = values$sampleBidx, dataSetId = values$dataSetID,
              AIBSARNAid = values$AIBSARNAid)
            incProgress(amount = .3)
            values$sampleBeuc <- getSimScores(values$sampleB, 
              relevantGenes = values$relevantGenes, 
              similarity_method = "euclidean")
            incProgress(amount = .3)
            values$sampleBcos <- getSimScores(values$sampleB,
              relevantGenes = values$relevantGenes,
              similarity_method = "cosine")
            incProgress(amount = .3)
          })
      }
      incProgress(amount = .3)
      if(isTruthy(input$sampleC)){
        withProgress(message = 'Getting similarity vectors for Sample C',
          value = .1, {
            values$sampleCidx <- which(
              rownames(pData(values$eset)) %in% input$sampleC, arr.ind = TRUE)
            values$sampleC <- getExternalVector(values$eset, 
              index = values$sampleCidx, dataSetId = values$dataSetID,
              AIBSARNAid = values$AIBSARNAid)
            incProgress(amount = .3)
            values$sampleCeuc <- getSimScores(values$sampleC, 
              relevantGenes = values$relevantGenes, 
              similarity_method = "euclidean")
            incProgress(amount = .3)
            values$sampleCcos <- getSimScores(values$sampleC,
              relevantGenes = values$relevantGenes,
              similarity_method = "cosine")
            incProgress(amount = .3)
        })
      }
      incProgress(amount = .3)
    })
    withProgress(message = 'Producing Similarity Matrices', value = .1,{
      values$sampleAeucmat <- getSimMatrix(values$sampleAeuc,
                                           values$relevantGenes)
      incProgress(amount = .1)
      values$sampleAcosmat <- getSimMatrix(values$sampleAcos,
                                           values$relevantGenes)
      incProgress(amount = .2)
      values$sampleBeucmat <- getSimMatrix(values$sampleBeuc,
                                           values$relevantGenes)
      incProgress(amount = .1)
      values$sampleBcosmat <- getSimMatrix(values$sampleBcos,
                                           values$relevantGenes)
      incProgress(amount = .2)
      values$sampleCeucmat <- getSimMatrix(values$sampleCeuc,
                                           values$relevantGenes)
      incProgress(amount = .1)
      values$sampleCcosmat <- getSimMatrix(values$sampleCcos,
                                           values$relevantGenes)
      incProgress(amount = .2)
    })
    withProgress(message = 'Generating Heatmaps', value = .1, {
      values$sampleAeucmap <- heatmapr(values$sampleAeucmat, Rowv = FALSE,
                                       Colv = FALSE)
      incProgress(amount = .1)
      values$sampleAcosmap <- heatmapr(values$sampleAcosmat, Rowv = FALSE,
                                       Colv = FALSE)
      incProgress(amount = .2)
      values$sampleBeucmap <- heatmapr(values$sampleBeucmat, Rowv = FALSE,
                                       Colv = FALSE)
      incProgress(amount = .1)
      values$sampleBcosmap <- heatmapr(values$sampleBcosmat, Rowv = FALSE,
                                       Colv = FALSE)
      incProgress(amount = .2)
      values$sampleCeucmap <- heatmapr(values$sampleCeucmat, Rowv = FALSE,
                                       Colv = FALSE)
      incProgress(amount = .1)
      values$sampleCcosmap <- heatmapr(values$sampleCcosmat, Rowv = FALSE,
                                       Colv = FALSE)
      incProgress(amount = .2)
    })
  })
  output$heatmapDisplay <- renderUI({
    if(!isTruthy(values$sampleAeucmap)){
      NULL
    } else {
      tabBox(width = 12,
        tabPanel("Sample A",
          fluidRow(
            box(width = 12, title = "Euclidean Distance Heatmap",
              plotlyOutput("sampleAeucOut")
            )
          ),
          fluidRow(
            box(width = 12, title = "Cosine Similarity Heatmap",
              plotlyOutput("sampleAcosOut")
            )
          )
        ),
        tabPanel("Sample B",
          fluidRow(
            box(width = 12, title = "Euclidean Distance Heatmap",
              plotlyOutput("sampleBeucOut")
            )
          ),
          fluidRow(
            box(width = 12, title = "Cosine Similarity Heatmap",
              plotlyOutput("sampleBcosOut")
            )
          )
        ),
        tabPanel("Sample C",
          fluidRow(
            box(width = 12, title = "Euclidean Distance Heatmap",
              plotlyOutput("sampleCeucOut")
            )
          ),
          fluidRow(
            box(width = 12, title = "Cosine Similarity Heatmap",
              plotlyOutput("sampleCcosOut")
            )
          )
        )
      )
    }
  })
  output$sampleAeucOut <- renderPlotly({
    heatmaply(values$sampleAeucmap)
  })
  output$sampleAcosOut <- renderPlotly({
    heatmaply(values$sampleAcosmap)
  })
  output$sampleBeucOut <- renderPlotly({
    heatmaply(values$sampleBeucmap)
  })
  output$sampleBcosOut <- renderPlotly({
    heatmaply(values$sampleBcosmap)
  })
  output$sampleCeucOut <- renderPlotly({
    heatmaply(values$sampleCeucmap)
  })
  output$sampleCcosOut <- renderPlotly({
    heatmaply(values$sampleCcosmap)
  })
})
