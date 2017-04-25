#' getExternalVector
#'
#' Get a sample vector from an outside Biobase ExpressionSet, for use in
#' creating subsets of AIBSARNA with \code{getRelevantGenes}
#'
#' @param dataSet a Biobase ExpressionSet
#' @param index the integer index of the sample of dataSet to be used
#' @param dataSetColName the name of the column of fData(dataSet) to be used as
#'     the sample vector's names, defaults to "gene_symbol"
#' @param AIBSARNAcolName the name of the column of fData(AIBSARNA) that is
#'     comparable to dataSetColName.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return
#' @export
#' @import Biobase
#'
#' @examples
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes, gene_names = "HGNC")
#' myGeneSampleVector <- getExternalVector(myGeneSet, index = 1,
#'      dataSetColName = "gene_symbol", AIBSARNAcolName = "gene_symbol")
#'
getExternalVector <- function(dataSet, index = 1, dataSetColName = "gene_symbol",
                              AIBSARNAcolName = c("gene_id", "ensembl_gene_id",
                                                  "gene_symbol", "entrez_id",
                                                  "refseq_ids")){
  v <- exprs(dataSet)[, index]
  names(v) <- as.character(fData(dataSet)[[dataSetColName]])
  # get gene identifiers common to v and AIBSARNA
  vInAIBSARNA <- which(fData(AIBSARNA::AIBSARNA)[[AIBSARNAcolName]] %in%
                         names(v), arr.ind = TRUE)
  vInAIBSARNA <- as.character(
                      fData(AIBSARNA::AIBSARNA)[[AIBSARNAcolName]][vInAIBSARNA])
  # get indices of v that are present in vInAIBSARNA
  genesToKeep <- which(names(v) %in% vInAIBSARNA, arr.ind = TRUE)
  # update v to only include comparable genes (genes in vInAIBSARNA)
  v <- v[genesToKeep]
  # remove any duplicate genes
  genesToKeep <- unique(names(v))
  v <- v[genesToKeep]
  return(v)
}
