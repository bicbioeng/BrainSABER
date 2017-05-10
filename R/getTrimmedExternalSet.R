#' Get a trimmed version of an external data set
#'
#' Returns an ExpressionSet that is a subset of \code{dataSet} containing only
#' genes that are present in AIBSARNA, for use in \code{getSimScores} or
#' \code{getUNDmatrix}.
#'
#' @param dataSet a Biobase ExpressionSet compatible data set
#' @param dataSetId the name of the column of gene identifiers in fData(dataSet)
#'     to be used to compare dataSet to AIBSARNA.
#' @param AIBSARNAid the name of the column of fData(AIBSARNA) that is
#'     comparable to dataSetId.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return a Biobase ExpressionSet
#' @export
#' @import Biobase
#' @import AIBSARNA
#'
#' @examples
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNAid = "gene_symbol")
#' myTrimmedGeneSet <- getTrimmedExternalSet(myGeneSet,
#'      dataSetId = "gene_symbol", AIBSARNAid = "gene_symbol")

getTrimmedExternalSet <- function(dataSet, dataSetId = "gene_symbol",
                           AIBSARNAid = c("gene_id", "ensembl_gene_id",
                            "gene_symbol", "entrez_id", "refseq_ids")){
  # get a trimmed vector faciliate things
  v <- getExternalVector(dataSet = dataSet, index = 1, dataSetId = dataSetId,
                         AIBSARNAid = AIBSARNAid)

  vInd <- which(fData(dataSet)[[dataSetId]] %in% names(v),
                arr.ind = TRUE)
  # subset feature data
  relfd <- fData(dataSet)[vInd, ]
  # get index(es) of any duplicates
  dup <- anyDuplicated(relfd[[dataSetId]])
  # if there are duplicates, remove them from vInd and regenerate relfd
  if (dup > 0) {
    vInd <- vInd[-dup]
    relfd <- fData(dataSet)[vInd, ]
  }
  # remove any unused factor levels
  relfd <- as.data.frame(apply(relfd, 2, function(x) {x[drop = TRUE]}))
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
