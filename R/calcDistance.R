#' Calculate Correlation Between User Dataset and Allen Brain Institute in Shiny App
#'
#' This function takes a CellScabbard object as returned by
#' \code{mapIDs} and returns list of matrices/data frame for each sample. This list
#' includes the unordered dataframe of correlation values, significance, and metadata;
#' and one matrix each for Kendall's tau and Spearman's rho correlation coefficients across
#' brain structures and ages in the Allen Brain Developmental Transcriptome. This function
#' is to be used in the context of the Shiny application workflow.
#'
#' @param dataSet a CellScabbard object with non-empty relevantGenes slot
#'
#' @return a list of correlation data frames and matrices, one list element for each sample
#' 
#' @import SummarizedExperiment
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @export

calcDistance <- function(dataSet){
  
  mappedexprsData <- assay(dataSet)
  rownames(mappedexprsData) <- rowData(dataSet)$gene
  
  #subset Allen data to remove features not in gene list
  allen <- assay(relevantGenes(dataSet))
  rownames(allen) <- rowData(relevantGenes(dataSet))[AIBSARNAid(dataSet)][[1]]
  allenMappedSubset <- allen[rownames(mappedexprsData),]
  
  #set up cosine distance function
  cosine <- function(v1, v2){
    v1Dotv2 <- sum(v1*v2)
    distv1 <- sqrt(sum(v1^2))
    distv2 <- sqrt(sum(v2^2))
    v1Dotv2/(distv1*distv2)
  }
  
  #average expression of regions that have the same age/region
  allenColsMeta <- colData(relevantGenes(dataSet))@listData
  allenMetaCollapsed <- paste(allenColsMeta$structure_name, allenColsMeta$age, sep = "_")
  allenSummarized <- tibble(allenMetaCollapsed = allenMetaCollapsed,
                            allen = t(allenMappedSubset))
  
  allenNested <- allenSummarized %>%
    group_by(allenMetaCollapsed) %>%
    tidyr::nest() %>%
    mutate(allenMeanData = purrr::map(.x = data, .f = colSums))
  allenSubset <- matrix(unlist(allenNested$allenMeanData), byrow = F, nrow = nrow(allenMappedSubset),
                        dimnames = list(row = rownames(allenMappedSubset),
                                        col = allenNested$allenMetaCollapsed))
  
  
  #ages and regions
  ages <- unique(allenColsMeta$age)
  regions <- unique(allenColsMeta$structure_name)
  
  #run distance and correlation functions
  distances <- apply(as.matrix(mappedexprsData), 2, function(n1){
    
    #####cosine distance
    distance <- apply(allenSubset, 2, function(n2){
      cosine(n1,n2)
    })
    
    #####kendall's tau
    correlation <- apply(allenSubset, 2, function(n3){
      cor <- cor.test(n1, n3, method = "kendall")
      c(cor$estimate, cor$p.value)
    })
    
    #####spearman's rho
    scorrelation <- apply(allenSubset, 2, function(n4){
      scor <- cor.test(n1, n4, method = "spearman")
      c(scor$estimate, scor$p.value)
    })
    
    #####combine results into a matrix and round
    allenDistance <- apply(cbind(distance, t(correlation), t(scorrelation)), 2, function(r){
      sapply(r, signif, digits = 4)
    })
    colnames(allenDistance) <- c('Cos', 'Tau', 'Tau pVal', 'Rho', 'Rho pVal')
    uniqueAllenColsMeta <- matrix(unlist(strsplit(allenNested$allenMetaCollapsed, "_")),
                                  ncol = 2, byrow = T)
    allenDistanceMeta <- cbind.data.frame(allenDistance, as.data.frame(uniqueAllenColsMeta))
    colnames(allenDistanceMeta)[6:7] <- c("structure_name", "age")
    #allenDistanceMeta <- allenDistanceMeta[,c(1:5, 9, 13)]
    
    #Kendall's tau heatmap matrix
    correlationMatrix <- matrix(0, ncol = length(ages), nrow = length(regions), dimnames = list(regions, ages))
    correlationMatrix[as.matrix(allenDistanceMeta[c("structure_name", "age")])] <- allenDistanceMeta[["Tau"]]
    
    #Spearman's rho heatmap matrix
    scorrelationMatrix <- matrix(0, ncol = length(ages), nrow = length(regions), dimnames = list(regions, ages))
    scorrelationMatrix[as.matrix(allenDistanceMeta[c("structure_name", "age")])] <- allenDistanceMeta[["Rho"]]
    
    ######prepare output
    list("correlation" = correlationMatrix,
         "scorrelation" = scorrelationMatrix,
         "meta" = allenDistanceMeta)
  })
  
  names(distances) <- colnames(mappedexprsData)
  distances
}