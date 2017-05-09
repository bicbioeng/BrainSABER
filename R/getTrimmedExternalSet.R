#' Get a trimmed version of an external data set
#'
#' Returns an ExpressionSet that is a subset of \code{dataSet} containing only
#' genes that are present in AIBSARNA, for use in \code{getUNDmatrix}.
#'
#' @param dataSet a Biobase ExpressionSet
#' @param dataSetColName the name of the column of fData(dataSet) to be used as
#'     the sample vector's names, defaults to "gene_symbol"
#' @param AIBSARNAcolName the name of the column of fData(AIBSARNA) that is
#'     comparable to dataSetColName.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return a Biobase ExpressionSet
#' @export
#'
#' @examples
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes, gene_names = "HGNC")
#' myTrimmedGeneSet <- getTrimmedExternalSet(myGeneSet,
#'      dataSetColName = "gene_symbol", AIBSARNAcolName = "gene_symbol")

getTrimmedExternalSet <- function(dataSet, dataSetColName = "gene_symbol",
                           AIBSARNAcolName = c("gene_id", "ensembl_gene_id",
                            "gene_symbol", "entrez_id", "refseq_ids")){
  # get a trimmed vector faciliate things
  v <- getExternalVector(dataSet = dataSet, index = 1, dataSetColName = dataSetColName, AIBSARNAcolName = AIBSARNAcolName)

  vInd <- which(fData(dataSet)[[dataSetColName]] %in% names(v),
                arr.ind = TRUE)
  # subset feature data
  relfd <- fData(dataSet)[vInd, ]
  # get index(es) of any duplicates
  dup <- anyDuplicated(relfd[[dataSetColName]])
  # if there are duplicates, remove them from vInd and regenerate relfd
  if (dup > 0) {
    vInd <- vInd[-dup]
    relfd <- fData(dataSet)[vInd, ]
  }
  # convert to Annotated Data Frame
  relfd <- new("AnnotatedDataFrame", data = relfd)
  # subset exprs
  relexprs <- exprs(dataSet)[vInd, ]
  # put together ExpressionSet
  relevantGenes <- ExpressionSet(assayData = relexprs,
                               phenoData = phenoData(dataSet),
                               featureData = relfd)
  return(relevantGenes)
}
